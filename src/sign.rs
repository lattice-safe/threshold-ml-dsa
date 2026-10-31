//! Party-side signing logic for the threshold ML-DSA protocol.
//!
//! Implements the 3-round signing flow from ePrint 2026/013:
//!
//! **Round 1 (Commit):** Each party i generates K commitments
//!   w_{i,k} = `A·r_k` + `e_k` using hyperball sampling, and broadcasts
//!   H(tr ‖ id ‖ w_{i,k}).
//!
//! **Round 2 (Reveal):** Each party reveals the full w_{i,k} vectors.
//!   Other parties verify the commitment hash.
//!
//! **Round 3 (Respond):** Each party computes K partial responses
//!   z_{i,k} = (`c·s₁_I` + `y_k`, `c·s₂_I` + `e_k`) and applies the
//!   hyperball rejection test ‖z_{i,k}‖₂ ≤ r.
//!   Responses that pass are nonzero; those that fail are zero vectors.
//!
//! The coordinator then aggregates and tries each of the K slots.

extern crate alloc;
use alloc::vec::Vec;

use crate::error::Error;
use crate::fvec::{sample_hyperball, FVec};
use crate::nonce::NonceGuard;
use crate::params::{CTILDEBYTES, MAX_PARTIES, TRBYTES, ThresholdParams, K, L, N, Q};
use crate::partition;
use crate::poly::{Poly, PolyVecK, PolyVecL};
use crate::rss::ThresholdPrivateKey;
use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use subtle::ConstantTimeEq;
use zeroize::Zeroize;

/// Round 1 state saved by a party for use in Round 3.
///
/// Contains sensitive ephemeral randomness; zeroized on drop.
/// Fields are `pub(crate)` to prevent external duplication of nonce state.
pub struct StRound1 {
    /// K packed commitment vectors (each is K polynomials in normal domain)
    pub(crate) w_packed: Vec<Vec<u8>>,
    /// K randomness samples as floating vectors from `SampleHyperball`.
    ///
    /// Round 1 commitments use rounded copies of these vectors.
    /// Round 3 adds them in floating-point before the final rounding step.
    pub(crate) rand_fvecs: Vec<FVec>,
    /// [`ThresholdPrivateKey::guard_key_id`] of the key that drew the nonce;
    /// `round3` refuses to answer it with another key.
    pub(crate) key_id: [u8; 32],
}

impl Drop for StRound1 {
    fn drop(&mut self) {
        // w_packed contains commitment data — not secret, but wipe defensively.
        for v in &mut self.w_packed {
            v.zeroize();
        }
        // rand_fvecs are zeroized by FVec's own Drop impl when the Vec drops.
    }
}

/// Round 2 state saved by a party for use in Round 3.
///
/// Owns the Round 1 state: `round2` consumes `StRound1` by value, so a party
/// can reveal its commitment for a given nonce at most once, and `round3`
/// consumes this state by value, so a nonce answers at most one challenge.
/// Fields are `pub(crate)` to prevent external state duplication.
#[derive(Zeroize)]
#[zeroize(drop)]
pub struct StRound2 {
    /// Commitment hashes received from all parties in Round 1
    pub(crate) hashes: Vec<[u8; 32]>,
    /// μ = CRH(tr ‖ msg)
    pub(crate) mu: [u8; 64],
    /// Active signer bitmask
    pub(crate) act: u8,
    /// Session ID the Round 1 commitment was bound to
    pub(crate) session_id: [u8; 32],
    /// Round 1 nonce state (zeroized by its own Drop impl)
    #[zeroize(skip)]
    pub(crate) st1: StRound1,
}

/// Proof-of-verification witness for Round-2 reveals.
///
/// This type can only be constructed by calling
/// [`verify_all_round2_reveals`], which checks every peer's reveal
/// against its Round-1 commitment hash. `round3` requires this
/// token, making it impossible to compute responses without first
/// verifying transcript consistency.
///
/// The witness is **cryptographically bound** to the specific transcript
/// it verified: it contains a binding hash over (session_id, act,
/// commitment hashes, and the aggregated commitments (`wfinals`) derived
/// from those reveals. `round3` re-derives and checks these bindings, so a
/// witness obtained from one reveal set cannot be replayed against a
/// different `StRound2`/`wfinals`.
///
/// `round3` consumes this witness **by value**, preventing reuse across
/// multiple calls or sessions.
pub struct VerifiedReveals {
    /// Cryptographic binding of the verified transcript.
    /// H("th-ml-dsa-verified-reveals-v1" ‖ session_id ‖ act ‖ hashes ‖ reveals)
    binding: [u8; 32],
    /// Session ID this witness was created for.
    session_id: [u8; 32],
    /// Number of parallel slots K used when deriving the aggregate.
    k_reps: u16,
    /// Hash of the expected aggregate commitments derived from the
    /// verified reveals.
    wfinals_hash: [u8; 32],
}

/// Verify all Round-2 reveals against their Round-1 commitment hashes.
///
/// Returns a [`VerifiedReveals`] witness token on success, which is
/// required by [`round3`]. This ensures distributed callers cannot
/// skip per-peer verification.
///
/// # Errors
///
/// Returns `Error::InvalidShare` if any reveal fails verification.
pub fn verify_all_round2_reveals(
    tr: &[u8; TRBYTES],
    party_ids: &[u8],
    act: u8,
    msg: &[u8],
    session_id: &[u8; 32],
    reveals: &[Vec<u8>],
    expected_hashes: &[[u8; 32]],
    k_reps: usize,
) -> Result<VerifiedReveals, Error> {
    if party_ids.len() != reveals.len() || reveals.len() != expected_hashes.len() {
        return Err(Error::InvalidParameters);
    }
    if k_reps == 0 || k_reps > u16::MAX as usize {
        return Err(Error::InvalidParameters);
    }
    validate_party_ids_match_act(party_ids, act)?;
    for (i, reveal) in reveals.iter().enumerate() {
        if !verify_round2_reveal(
            tr,
            party_ids[i],
            act,
            msg,
            session_id,
            reveal,
            &expected_hashes[i],
        ) {
            return Err(Error::InvalidShare);
        }
    }

    // Re-derive the aggregate commitments from verified reveal payloads and
    // seal their hash into the witness. round3 will validate coordinator input
    // against this hash to prevent aggregate-substitution by low-level callers.
    let wfinals_hash = compute_wfinals_hash_from_reveals(reveals, k_reps)?;
    let k_reps_u16 = k_reps as u16;

    // Compute a cryptographic binding over the verified transcript so that
    // this witness cannot be replayed with different hashes or aggregate.
    let binding =
        compute_verified_reveals_binding(session_id, act, expected_hashes, k_reps_u16, &wfinals_hash);

    Ok(VerifiedReveals {
        binding,
        session_id: *session_id,
        k_reps: k_reps_u16,
        wfinals_hash,
    })
}

/// Compute μ = CRH(tr ‖ msg).
#[must_use] 
pub fn compute_mu(tr: &[u8; TRBYTES], msg: &[u8]) -> [u8; 64] {
    let mut h = Shake256::default();
    h.update(tr);
    // ML-DSA context encoding: 0x00 ‖ len(ctx) ‖ ctx ‖ msg
    h.update(&[0u8]); // attribute byte
    h.update(&[0u8]); // ctx length = 0 (no context)
    h.update(msg);
    let mut mu = [0u8; 64];
    h.finalize_xof().read(&mut mu);
    mu
}

/// Round 1: generate K hyperball commitments.
///
/// # Security
///
/// `guard` reserves `session_id` before any randomness is drawn and refuses
/// a session ID already used with this key; its counter is mixed into the
/// nonce seed, so a repeated RNG output alone cannot reproduce a nonce.
/// Reusing a nonce under two challenges would reveal this party's partial
/// secret. See [`crate::nonce`] for what a guard cannot protect against.
///
/// Returns:
/// - `commitment_hash`: 32-byte hash H(tag ‖ tr ‖ id ‖ act ‖ session ‖ μ ‖ `w_packed`)
/// - `StRound1`: saved state for Round 3
///
/// This matches the Go `Round1()` function.
pub fn round1<R: RngCore + CryptoRng, G: NonceGuard + ?Sized>(
    sk: &ThresholdPrivateKey,
    params: &ThresholdParams,
    act: u8,
    msg: &[u8],
    session_id: &[u8; 32],
    rng: &mut R,
    guard: &mut G,
) -> Result<([u8; 32], StRound1), Error> {
    params.validate()?;
    if (sk.t, sk.n) != (params.t, params.n) {
        return Err(Error::InvalidParameters);
    }
    if sk.id >= params.n
        || act >> params.n != 0
        || act & (1 << sk.id) == 0
        || act.count_ones() != u32::from(params.t)
    {
        return Err(Error::InvalidParameters);
    }
    // Reserve the session before drawing randomness: a replayed session ID
    // is refused, and the counter makes every nonce seed input distinct.
    let key_id = sk.guard_key_id();
    let counter = guard.reserve(&key_id, session_id)?;
    let mu = compute_mu(&sk.tr, msg);

    // Hedged nonce derivation: bind RNG entropy to the long-term per-party
    // key, transcript context and the guard counter.
    let mut rng_entropy = [0u8; 64];
    rng.fill_bytes(&mut rng_entropy);
    let mut rhop = [0u8; 64];
    let mut rhop_h = Shake256::default();
    rhop_h.update(b"th-ml-dsa-round1-rhop-v2");
    rhop_h.update(&rng_entropy);
    rhop_h.update(&sk.key);
    rhop_h.update(&sk.tr);
    rhop_h.update(&[sk.id, act]);
    rhop_h.update(session_id);
    rhop_h.update(&mu);
    rhop_h.update(&counter.to_le_bytes());
    rhop_h.finalize_xof().read(&mut rhop);
    rng_entropy.zeroize();

    let k_reps = params.k_reps as usize;
    let mut rand_fvecs: Vec<FVec> = Vec::with_capacity(k_reps);
    let mut w_packed_all: Vec<Vec<u8>> = Vec::with_capacity(k_reps);

    // Expand A from ρ
    let a = expand_a(&sk.rho);

    for i in 0..k_reps {
        // Sample in place: moving a filled FVec would leave an unwiped copy
        // of the nonce on the stack. (Capacity is reserved, so no realloc.)
        rand_fvecs.push(FVec::zero());
        let fv = &mut rand_fvecs[i];
        // Sample from hyperball with radius r₁ (randomness ball)
        sample_hyperball(fv, params.r1, params.nu, &rhop, i as u16);

        // Round FVec to integer polynomials (y, e_) ∈ ℤ^L × ℤ^K
        // CRITICAL: These rounded integers are what we use for BOTH
        // w = A·y + e_ (commitment) and z = c·s + y (response).
        let mut y = PolyVecL::zero();
        let mut e_ = PolyVecK::zero();
        fv.round_to_polyvecs(&mut y, &mut e_);

        // Compute w = A·y + e_
        let mut yh = y.clone();
        yh.ntt();
        let mut w = PolyVecK::zero();
        for (j, a_row) in a.iter().enumerate().take(K) {
            // w[j] = Σ_l A[j][l] · yh[l] (pointwise in NTT domain)
            // acc = (A·y)[j]; with the public w it would give e, so wipe it.
            let mut acc = Poly::zero();
            let mut prod = Poly::zero();
            for (l, a_jl) in a_row.iter().enumerate().take(L) {
                Poly::pointwise_montgomery(&mut prod, a_jl, &yh.polys[l]);
                for (a, p) in acc.coeffs.iter_mut().zip(&prod.coeffs) {
                    *a += p;
                }
            }
            acc.reduce();
            acc.invntt_tomont();
            // w[j] = A·y[j] + e_[j]
            for ((wc, a), e) in w.polys[j].coeffs.iter_mut().zip(&acc.coeffs).zip(&e_.polys[j].coeffs) {
                *wc = a + e;
            }
            acc.zeroize();
            prod.zeroize();
        }
        w.reduce();

        // Pack w for transport (23 bits per coeff × K×N coefficients)
        let packed = pack_w_single(&w);
        w_packed_all.push(packed);

        // y is the integer nonce: with an accepted z it gives c·s_I = z − y.
        y.zeroize();
        e_.zeroize();
        yh.zeroize();
    }
    rhop.zeroize();

    let hash = round1_commitment_hash(
        &sk.tr,
        sk.id,
        act,
        &mu,
        session_id,
        w_packed_all.iter().map(alloc::vec::Vec::as_slice),
    );

    Ok((
        hash,
        StRound1 {
            w_packed: w_packed_all,
            rand_fvecs,
            key_id,
        },
    ))
}

/// Round 2: reveal commitments and compute μ.
///
/// Consumes `st_rd1`: a party reveals its Round 1 commitment exactly once, so
/// it cannot be made to answer a transcript whose peer commitments were chosen
/// after its own reveal was seen.
///
/// Verifies that this party's own reveal is bound to its Round 1 commitment
/// and that this commitment appears in its own slot of `msgs_rd1`.
/// NOTE: Per-peer reveal verification is done separately via
/// `verify_round2_reveal` -- callers in distributed settings MUST
/// call it for every received reveal.
/// Returns the packed commitment data and state for Round 3.
#[allow(clippy::too_many_arguments)]
pub fn round2(
    sk: &ThresholdPrivateKey,
    act: u8,
    msg: &[u8],
    session_id: &[u8; 32],
    msgs_rd1: &[[u8; 32]], // commitment hashes from Round 1
    own_rd1_hash: &[u8; 32],
    st_rd1: StRound1,
    params: &ThresholdParams,
) -> Result<(Vec<u8>, StRound2), Error> {
    params.validate()?;
    if (sk.t, sk.n) != (params.t, params.n) {
        return Err(Error::InvalidParameters);
    }
    // The transcript must hold exactly one commitment per active party, and
    // our own slot must hold our own commitment; otherwise the aggregate
    // could exclude our w and leave the challenge under someone else's control.
    let active = bitmask_to_sorted_ids(act, params.n);
    if act >> params.n != 0 || msgs_rd1.len() != active.len() {
        return Err(Error::InvalidParameters);
    }
    let own_idx = active
        .iter()
        .position(|&id| id == sk.id)
        .ok_or(Error::InvalidParameters)?;
    if !bool::from(msgs_rd1[own_idx].ct_eq(own_rd1_hash)) {
        return Err(Error::InvalidShare);
    }

    // Concatenate all packed w data for this party's reveal
    let mut reveal_data = Vec::new();
    for packed in &st_rd1.w_packed {
        reveal_data.extend_from_slice(packed);
    }

    // Verify this reveal is bound to the Round 1 commitment in constant-time.
    let mu = compute_mu(&sk.tr, msg);
    let computed = round1_commitment_hash(
        &sk.tr,
        sk.id,
        act,
        &mu,
        session_id,
        st_rd1.w_packed.iter().map(alloc::vec::Vec::as_slice),
    );
    if !bool::from(computed.ct_eq(own_rd1_hash)) {
        return Err(Error::InvalidShare);
    }

    let st2 = StRound2 {
        hashes: msgs_rd1.to_vec(),
        mu,
        act,
        session_id: *session_id,
        st1: st_rd1,
    };

    Ok((reveal_data, st2))
}

/// Verify a Round 2 reveal against the sender's Round 1 commitment hash.
#[must_use] 
pub fn verify_round2_reveal(
    tr: &[u8; TRBYTES],
    party_id: u8,
    act: u8,
    msg: &[u8],
    session_id: &[u8; 32],
    reveal_data: &[u8],
    expected_hash: &[u8; 32],
) -> bool {
    let mu = compute_mu(tr, msg);
    let computed = round1_commitment_hash(
        tr,
        party_id,
        act,
        &mu,
        session_id,
        core::iter::once(reveal_data),
    );
    bool::from(computed.ct_eq(expected_hash))
}

/// Round 3: compute K parallel responses with hyperball rejection.
///
/// This matches the Go `ComputeResponses()` function.
/// For each of the K parallel commitments, computes:
///   `z_k` = `c·s_I` + (`r_k`, `e_k`) where (`r_k`, `e_k`) is from the hyperball sample
///
/// Then applies the ν-scaled L₂ norm check: ‖`z_k‖₂` ≤ r.
/// Responses that fail rejection are all-zero (coordinator will skip them).
///
/// # State Consumption
///
/// Takes `st_rd1` by value to enforce one-time use of the nonce randomness.
/// Calling this function consumes the Round 1 state, preventing replay of
/// the same hidden randomness with different coordinator-selected challenges.
///
/// # Transcript Verification (Distributed Callers)
///
/// This function does not receive raw per-party reveals, so it cannot
/// directly re-verify each reveal/hash pair. Instead, it enforces aggregate
/// consistency by checking a hash of `wfinals` against the
/// [`VerifiedReveals`] witness produced by [`verify_all_round2_reveals`].
///
/// Distributed integrators MUST obtain a [`VerifiedReveals`] token by
/// calling [`verify_all_round2_reveals`] before calling this function.
/// The SDK does this automatically.
pub fn round3(
    sk: &ThresholdPrivateKey,
    wfinals: &[PolyVecK], // K aggregated commitment vectors
    st_rd2: StRound2,     // consumed: owns the single-use Round 1 nonce state
    params: &ThresholdParams,
    verified: VerifiedReveals,  // consumed: proof that all reveals were checked
) -> Result<Vec<PolyVecL>, Error> {
    use dilithium::{
        poly::Poly as DPoly,
        polyvec::{polyveck_reduce, polyvecl_reduce, PolyVecK as DPolyVecK, PolyVecL as DPolyVecL},
        ML_DSA_44,
    };

    params.validate()?;
    if (sk.t, sk.n) != (params.t, params.n) || st_rd2.act >> params.n != 0 {
        return Err(Error::InvalidParameters);
    }
    // The nonce must be answered by the key that drew it.
    if !bool::from(st_rd2.st1.key_id.ct_eq(&sk.guard_key_id())) {
        return Err(Error::InvalidParameters);
    }
    let k_reps = params.k_reps as usize;
    let mode = ML_DSA_44;

    // Reject short coordinator input instead of silently truncating.
    if wfinals.len() < k_reps {
        return Err(Error::InvalidParameters);
    }
    if verified.k_reps as usize != k_reps {
        return Err(Error::InvalidParameters);
    }
    if !bool::from(verified.session_id.ct_eq(&st_rd2.session_id)) {
        return Err(Error::InvalidShare);
    }

    // pack_w_single is only injective on [0, Q): reject non-canonical
    // coefficients so an aliased aggregate cannot pass the hash check below
    // while decomposing to a different w₁ (and thus a different challenge).
    for w in &wfinals[..k_reps] {
        for p in &w.polys {
            if p.coeffs.iter().any(|&c| !(0..Q).contains(&c)) {
                return Err(Error::InvalidParameters);
            }
        }
    }

    // Ensure the coordinator-supplied aggregates match the reveals that were
    // previously verified in Round 2.
    let provided_wfinals_hash = compute_wfinals_hash(&wfinals[..k_reps]);
    if !bool::from(verified.wfinals_hash.ct_eq(&provided_wfinals_hash)) {
        return Err(Error::InvalidShare);
    }

    // Verify that the witness is bound to this exact transcript.
    // Re-derive the expected binding from (session_id, act, hashes, k_reps,
    // wfinals_hash) and
    // compare against the one sealed in the witness. The commitment hashes
    // already commit to the individual reveals.
    {
        let expected_binding = compute_verified_reveals_binding(
            &verified.session_id,
            st_rd2.act,
            &st_rd2.hashes,
            verified.k_reps,
            &verified.wfinals_hash,
        );
        if !bool::from(verified.binding.ct_eq(&expected_binding)) {
            return Err(Error::InvalidShare);
        }
    }

    // Recover this party's partial secret for the active signer set
    let active = bitmask_to_sorted_ids(st_rd2.act, params.n);
    let partition = partition::rss_recover(&active, params.n, params.t)?;

    // Find which partition slot corresponds to this party
    let party_idx = active
        .iter()
        .position(|&id| id == sk.id)
        .ok_or(Error::InvalidParameters)?;

    // Sum the shares assigned to this party by the partition.
    // NOTE: recover_partial_secret returns values ALREADY in NTT domain
    // because Share stores s1h = NTT(s1), so summing them gives NTT(Σ s1)
    // directly. Do NOT NTT again — that would double-transform.
    let (mut s1h, mut s2h) = recover_partial_secret(sk, &partition[party_idx])?;
    let mut s1h_d = DPolyVecL::default();
    let mut s2h_d = DPolyVecK::default();
    for j in 0..L {
        for c in 0..N {
            s1h_d.vec[j].coeffs[c] = s1h.polys[j].coeffs[c];
        }
    }
    for j in 0..K {
        for c in 0..N {
            s2h_d.vec[j].coeffs[c] = s2h.polys[j].coeffs[c];
        }
    }

    let mut zs: Vec<PolyVecL> = Vec::with_capacity(k_reps);

    for (i, wfinal) in wfinals.iter().enumerate().take(k_reps) {
        // Decompose w — only w₁ is needed for the challenge hash.
        let w1 = decompose_w1_only(wfinal);

        // c̃ = H(μ ‖ w₁_packed)
        let c_tilde = compute_challenge(&st_rd2.mu, &w1);

        // Expand challenge polynomial c
        let mut ch = DPoly::default();
        DPoly::challenge(mode, &mut ch, &c_tilde);
        ch.ntt();

        // Compute z1 = c·s₁ (integer arithmetic in NTT domain, then back)
        let mut z1 = PolyVecL::zero();
        let mut z1_d = DPolyVecL::default();
        for j in 0..L {
            DPoly::pointwise_montgomery(&mut z1_d.vec[j], &ch, &s1h_d.vec[j]);
            z1_d.vec[j].invntt_tomont();
        }
        polyvecl_reduce(mode, &mut z1_d);
        for j in 0..L {
            for c in 0..N {
                z1.polys[j].coeffs[c] = z1_d.vec[j].coeffs[c];
            }
        }

        // Compute z2 = c·s₂.
        let mut z2 = PolyVecK::zero();
        let mut z2_d = DPolyVecK::default();
        for j in 0..K {
            DPoly::pointwise_montgomery(&mut z2_d.vec[j], &ch, &s2h_d.vec[j]);
            z2_d.vec[j].invntt_tomont();
        }
        polyveck_reduce(mode, &mut z2_d);
        for j in 0..K {
            for c in 0..N {
                z2.polys[j].coeffs[c] = z2_d.vec[j].coeffs[c];
            }
        }

        // Build the full floating response vector and add the sampled
        // hyperball randomness from Round 1 before rejection/rounding.
        let mut zf = FVec::from_polyvecs(&z1, &z2);
        // c·s_I is secret; zf (zeroized on drop) now holds the only copy.
        z1.zeroize();
        z2.zeroize();
        for p in z1_d.vec.iter_mut().chain(z2_d.vec.iter_mut()) {
            p.coeffs.zeroize();
        }
        zf.add_assign(&st_rd2.st1.rand_fvecs[i]);
        if zf.excess(params.r, params.nu) {
            // Rejection: emit zero response
            zs.push(PolyVecL::zero());
            continue;
        }

        // Accepted — round zf back to integer polynomial vectors.
        // By protocol design we only transmit the z1 part (PolyVecL); z2 is
        // used in local rejection algebra but is not part of the final FIPS 204
        // signature encoding.
        let mut z_out = PolyVecL::zero();
        let mut z2_out = PolyVecK::zero();
        zf.round_to_polyvecs(&mut z_out, &mut z2_out);
        zs.push(z_out);
        z2_out.zeroize();
    }

    s1h.zeroize();
    s2h.zeroize();
    for p in s1h_d.vec.iter_mut().chain(s2h_d.vec.iter_mut()) {
        p.coeffs.zeroize();
    }

    Ok(zs)
}

// ─── Internal Helpers ────────────────────────────────────────────────

/// Expand A from ρ using dilithium-rs.
fn expand_a(rho: &[u8; 32]) -> Vec<Vec<Poly>> {
    use dilithium::{
        polyvec::{matrix_expand, PolyVecL as DPolyVecL},
        ML_DSA_44,
    };
    let mode = ML_DSA_44;
    let mut mat: [DPolyVecL; K] = core::array::from_fn(|_| DPolyVecL::default());
    matrix_expand(mode, &mut mat, rho);

    // Convert to our Poly format
    let mut a: Vec<Vec<Poly>> = Vec::with_capacity(K);
    for mat_row in mat.iter().take(K) {
        let mut out_row = Vec::with_capacity(L);
        for j in 0..L {
            let mut p = Poly::zero();
            for c in 0..N {
                p.coeffs[c] = mat_row.vec[j].coeffs[c];
            }
            out_row.push(p);
        }
        a.push(out_row);
    }
    a
}

/// Size of a single packed `PolyVecK` (23 bits per coefficient, K×N coefficients).
#[must_use] 
pub fn pack_w_single_size() -> usize {
    let poly_q_size = (N * 23).div_ceil(8);
    K * poly_q_size
}

/// Pack a single `PolyVecK` as 23 bits per coefficient.
#[must_use] 
pub fn pack_w_single(w: &PolyVecK) -> Vec<u8> {
    let poly_q_size = (N * 23).div_ceil(8);
    let total = K * poly_q_size;
    let mut buf = alloc::vec![0u8; total];

    let mut bit_pos = 0usize;
    for poly in &w.polys {
        for &coeff in &poly.coeffs {
            let val = if coeff < 0 {
                (coeff + Q) as u32
            } else {
                coeff as u32
            };
            // Write 23 bits at bit_pos
            let byte_pos = bit_pos / 8;
            let bit_off = bit_pos % 8;
            let v = u64::from(val) << bit_off;
            for b in 0..4 {
                if byte_pos + b < buf.len() {
                    buf[byte_pos + b] |= (v >> (b * 8)) as u8;
                }
            }
            bit_pos += 23;
        }
    }

    buf
}

/// Unpack a single `PolyVecK` from 23-bit packed data.
#[must_use] 
pub fn unpack_w_single(buf: &[u8]) -> PolyVecK {
    let mut w = PolyVecK::zero();
    let mut bit_pos = 0usize;
    for poly in &mut w.polys {
        for coeff in &mut poly.coeffs {
            let byte_pos = bit_pos / 8;
            let bit_off = bit_pos % 8;
            let mut v: u64 = 0;
            for b in 0..4 {
                if byte_pos + b < buf.len() {
                    v |= u64::from(buf[byte_pos + b]) << (b * 8);
                }
            }
            *coeff = ((v >> bit_off) & 0x7FFFFF) as i32;
            bit_pos += 23;
        }
    }
    w
}

/// Convert bitmask to sorted list of party IDs.
fn bitmask_to_sorted_ids(mask: u8, n: u8) -> Vec<u8> {
    let mut ids = Vec::new();
    for i in 0..n {
        if mask & (1 << i) != 0 {
            ids.push(i);
        }
    }
    ids
}

/// Validate that `party_ids` is the canonical sorted party set encoded by `act`.
fn validate_party_ids_match_act(party_ids: &[u8], act: u8) -> Result<(), Error> {
    if party_ids.is_empty() {
        return Err(Error::InvalidParameters);
    }

    let mut reconstructed = 0u8;
    let mut prev: Option<u8> = None;
    for &id in party_ids {
        if id as usize >= MAX_PARTIES {
            return Err(Error::InvalidParameters);
        }
        if let Some(p) = prev {
            if id <= p {
                return Err(Error::InvalidParameters);
            }
        }
        reconstructed |= 1u8 << id;
        prev = Some(id);
    }

    if reconstructed != act {
        return Err(Error::InvalidParameters);
    }

    Ok(())
}

/// Compute a cryptographic binding for a [`VerifiedReveals`] witness.
///
/// H("th-ml-dsa-verified-reveals-v1" ‖ session_id ‖ act ‖ hashes ‖ k_reps ‖ wfinals_hash)
///
/// The commitment hashes commit to reveal data, and `wfinals_hash` commits
/// to the expected aggregate commitments derived from those verified reveals.
fn compute_verified_reveals_binding(
    session_id: &[u8; 32],
    act: u8,
    hashes: &[[u8; 32]],
    k_reps: u16,
    wfinals_hash: &[u8; 32],
) -> [u8; 32] {
    let mut h = Shake256::default();
    h.update(b"th-ml-dsa-verified-reveals-v1");
    h.update(session_id);
    h.update(&[act]);
    for hash in hashes {
        h.update(hash);
    }
    h.update(&k_reps.to_le_bytes());
    h.update(wfinals_hash);
    let mut out = [0u8; 32];
    h.finalize_xof().read(&mut out);
    out
}

/// Compute a hash of aggregate commitment vectors (`wfinals`).
///
/// H("th-ml-dsa-wfinals-v1" ‖ k_reps ‖ pack_w_single(wfinal[0]) ‖ ... )
fn compute_wfinals_hash(wfinals: &[PolyVecK]) -> [u8; 32] {
    let mut h = Shake256::default();
    h.update(b"th-ml-dsa-wfinals-v1");
    h.update(&(wfinals.len() as u16).to_le_bytes());
    for w in wfinals {
        let packed = pack_w_single(w);
        h.update(&packed);
    }
    let mut out = [0u8; 32];
    h.finalize_xof().read(&mut out);
    out
}

/// Re-derive and hash aggregate commitments from verified reveal payloads.
fn compute_wfinals_hash_from_reveals(reveals: &[Vec<u8>], k_reps: usize) -> Result<[u8; 32], Error> {
    let packed_size = pack_w_single_size();
    let expected_len = k_reps
        .checked_mul(packed_size)
        .ok_or(Error::InvalidParameters)?;
    for reveal in reveals {
        if reveal.len() != expected_len {
            return Err(Error::InvalidParameters);
        }
    }

    let mut wfinals = Vec::with_capacity(k_reps);
    for k in 0..k_reps {
        let start = k * packed_size;
        let end = start + packed_size;
        let mut acc = PolyVecK::zero();
        for reveal in reveals {
            let w = unpack_w_single(&reveal[start..end]);
            acc.add_assign(&w);
        }
        acc.reduce();
        wfinals.push(acc);
    }
    Ok(compute_wfinals_hash(&wfinals))
}

fn round1_commitment_hash<'a, I>(
    tr: &[u8; TRBYTES],
    party_id: u8,
    act: u8,
    mu: &[u8; 64],
    session_id: &[u8; 32],
    packed_chunks: I,
) -> [u8; 32]
where
    I: IntoIterator<Item = &'a [u8]>,
{
    let mut h = Shake256::default();
    h.update(b"th-ml-dsa-round1-commit-v1");
    h.update(tr);
    h.update(&[party_id, act]);
    h.update(session_id);
    h.update(mu);
    for chunk in packed_chunks {
        h.update(chunk);
    }
    let mut out = [0u8; 32];
    h.finalize_xof().read(&mut out);
    out
}

/// Recover this party's partial secret from the assigned share bitmasks.
///
/// Fails if the partition assigns a share this party does not hold.
fn recover_partial_secret(
    sk: &ThresholdPrivateKey,
    assigned_masks: &[u8],
) -> Result<(PolyVecL, PolyVecK), Error> {
    // Check first so an error never leaves a partial sum behind.
    if assigned_masks.iter().any(|m| !sk.shares.contains_key(m)) {
        return Err(Error::InvalidShare);
    }
    let mut s1h = PolyVecL::zero();
    let mut s2h = PolyVecK::zero();
    for share in assigned_masks.iter().filter_map(|m| sk.shares.get(m)) {
        s1h.add_assign(&share.s1h);
        s2h.add_assign(&share.s2h);
    }
    s1h.reduce();
    s2h.reduce();
    Ok((s1h, s2h))
}

/// Extract only w1 (`HighBits`) from a `PolyVecK`, skipping the unused w0 copy.
///
/// CRITICAL: This MUST produce bit-identical results to the coordinator's
/// `dilithium::rounding::decompose()`, otherwise the challenge hash will
/// differ between round3 and combine, causing 100% delta rejection.
fn decompose_w1_only(w: &PolyVecK) -> PolyVecK {
    use dilithium::{polyvec::polyveck_decompose, polyvec::PolyVecK as DPolyVecK, ML_DSA_44};
    let mode = ML_DSA_44;

    let mut w_ref = DPolyVecK::default();
    for i in 0..K {
        w_ref.vec[i].coeffs = w.polys[i].coeffs;
    }

    let mut w1_ref = DPolyVecK::default();
    let mut w0_ref = DPolyVecK::default();
    polyveck_decompose(mode, &mut w1_ref, &mut w0_ref, &w_ref);

    let mut w1 = PolyVecK::zero();
    for i in 0..K {
        w1.polys[i].coeffs = w1_ref.vec[i].coeffs;
    }
    w1
}

/// Compute challenge hash c̃ = H(μ ‖ `w₁_packed`).
fn compute_challenge(mu: &[u8; 64], w1: &PolyVecK) -> [u8; CTILDEBYTES] {
    use dilithium::{polyvec::polyveck_pack_w1, polyvec::PolyVecK as DPolyVecK, ML_DSA_44};
    let mode = ML_DSA_44;

    let mut w1_ref = DPolyVecK::default();
    for i in 0..K {
        w1_ref.vec[i].coeffs = w1.polys[i].coeffs;
    }
    let mut w1_packed = alloc::vec![0u8; mode.k() * mode.polyw1_packedbytes()];
    polyveck_pack_w1(mode, &mut w1_packed, &w1_ref);

    let mut h = Shake256::default();
    h.update(mu);
    h.update(&w1_packed);
    let mut c = [0u8; CTILDEBYTES];
    h.finalize_xof().read(&mut c);
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_unpack_w_roundtrip() {
        let mut w = PolyVecK::zero();
        for i in 0..K {
            for j in 0..N {
                w.polys[i].coeffs[j] = ((i * N + j) as i32 * 137) % Q;
            }
        }
        let packed = pack_w_single(&w);
        let recovered = unpack_w_single(&packed);
        for i in 0..K {
            for j in 0..N {
                assert_eq!(
                    recovered.polys[i].coeffs[j],
                    w.polys[i].coeffs[j],
                    "w pack/unpack mismatch at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn test_pack_w_single_size_matches() {
        let size = pack_w_single_size();
        let w = PolyVecK::zero();
        let packed = pack_w_single(&w);
        assert_eq!(packed.len(), size);
    }

    #[test]
    fn test_bitmask_to_sorted_ids_basic() {
        let ids = bitmask_to_sorted_ids(0b0111, 4);
        assert_eq!(ids, vec![0, 1, 2]);
        let ids = bitmask_to_sorted_ids(0b1010, 4);
        assert_eq!(ids, vec![1, 3]);
        let ids = bitmask_to_sorted_ids(0, 6);
        assert!(ids.is_empty());
    }

    #[test]
    fn test_verify_all_reveals_rejects_malformed_party_set() {
        let tr = [1u8; TRBYTES];
        let session_id = [2u8; 32];
        let msg = b"party set validation";
        let reveals = vec![Vec::new(), Vec::new()];
        let hashes = vec![[0u8; 32], [0u8; 32]];

        assert_eq!(
            verify_all_round2_reveals(&tr, &[0, 0], 0b0000_0011, msg, &session_id, &reveals, &hashes, 1)
                .err(),
            Some(Error::InvalidParameters)
        );
        assert_eq!(
            verify_all_round2_reveals(&tr, &[1, 0], 0b0000_0011, msg, &session_id, &reveals, &hashes, 1)
                .err(),
            Some(Error::InvalidParameters)
        );
        assert_eq!(
            verify_all_round2_reveals(&tr, &[0, 8], 129, msg, &session_id, &reveals, &hashes, 1)
                .err(),
            Some(Error::InvalidParameters)
        );
        assert_eq!(
            verify_all_round2_reveals(&tr, &[0, 2], 0b0000_0011, msg, &session_id, &reveals, &hashes, 1)
                .err(),
            Some(Error::InvalidParameters)
        );
    }

    #[test]
    fn test_compute_mu_deterministic() {
        let tr = [42u8; TRBYTES];
        let msg = b"test message";
        let mu1 = compute_mu(&tr, msg);
        let mu2 = compute_mu(&tr, msg);
        assert_eq!(mu1, mu2, "compute_mu must be deterministic");
    }

    #[test]
    fn test_compute_mu_different_messages() {
        let tr = [42u8; TRBYTES];
        let mu1 = compute_mu(&tr, b"message A");
        let mu2 = compute_mu(&tr, b"message B");
        assert_ne!(mu1, mu2, "different messages must produce different mu");
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;
    use crate::coordinator;
    use crate::nonce::InMemoryNonceGuard;
    use crate::params::get_threshold_params;
    use crate::rss::keygen_from_seed;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    const MSG: &[u8] = b"coverage";
    const SID: [u8; 32] = [7u8; 32];
    const ACT: u8 = 0b11;

    /// A (2, 2) session driven to the inputs of party 0's Round 3.
    struct Fixture {
        sks: Vec<ThresholdPrivateKey>,
        params: ThresholdParams,
        hashes: Vec<[u8; 32]>,
        reveals: Vec<Vec<u8>>,
        st2: StRound2,
        wfinals: Vec<PolyVecK>,
        verified: VerifiedReveals,
    }

    fn fixture() -> Fixture {
        let params = get_threshold_params(2, 2).unwrap();
        let (_, sks) = keygen_from_seed(&[1u8; 32], &params).unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        let mut guard = InMemoryNonceGuard::new();
        let mut hashes = Vec::new();
        let mut st1s = Vec::new();
        for sk in &sks {
            let (h, st1) = round1(sk, &params, ACT, MSG, &SID, &mut rng, &mut guard).unwrap();
            hashes.push(h);
            st1s.push(st1);
        }
        let mut reveals = Vec::new();
        let mut st2s = Vec::new();
        for (i, st1) in st1s.into_iter().enumerate() {
            let (r, st2) = round2(&sks[i], ACT, MSG, &SID, &hashes, &hashes[i], st1, &params).unwrap();
            reveals.push(r);
            st2s.push(st2);
        }
        let k = params.k_reps as usize;
        let size = pack_w_single_size();
        let parsed: Vec<Vec<PolyVecK>> = reveals
            .iter()
            .map(|r| (0..k).map(|j| unpack_w_single(&r[j * size..(j + 1) * size])).collect())
            .collect();
        let wfinals = coordinator::aggregate_commitments(&parsed, k).unwrap();
        let verified =
            verify_all_round2_reveals(&sks[0].tr, &[0, 1], ACT, MSG, &SID, &reveals, &hashes, k).unwrap();
        Fixture { sks, params, hashes, reveals, st2: st2s.into_iter().next().unwrap(), wfinals, verified }
    }

    #[test]
    fn happy_path_round3() {
        let f = fixture();
        let zs = round3(&f.sks[0], &f.wfinals, f.st2, &f.params, f.verified).unwrap();
        assert_eq!(zs.len(), f.params.k_reps as usize);
    }

    #[test]
    fn verify_all_rejects_bad_inputs() {
        let f = fixture();
        let k = f.params.k_reps as usize;
        let tr = &f.sks[0].tr;
        // length mismatch
        assert_eq!(
            verify_all_round2_reveals(tr, &[0], ACT, MSG, &SID, &f.reveals, &f.hashes, k).err(),
            Some(Error::InvalidParameters)
        );
        // k_reps out of range
        assert_eq!(
            verify_all_round2_reveals(tr, &[0, 1], ACT, MSG, &SID, &f.reveals, &f.hashes, 0).err(),
            Some(Error::InvalidParameters)
        );
        assert_eq!(
            verify_all_round2_reveals(tr, &[0, 1], ACT, MSG, &SID, &f.reveals, &f.hashes, 1 << 16).err(),
            Some(Error::InvalidParameters)
        );
        // tampered reveal
        let mut bad = f.reveals.clone();
        bad[1][0] ^= 1;
        assert_eq!(
            verify_all_round2_reveals(tr, &[0, 1], ACT, MSG, &SID, &bad, &f.hashes, k).err(),
            Some(Error::InvalidShare)
        );
        // party ids vs act
        assert!(validate_party_ids_match_act(&[], ACT).is_err());
        assert!(validate_party_ids_match_act(&[0, 2], ACT).is_err());
        // reveal of the wrong length
        assert_eq!(compute_wfinals_hash_from_reveals(&[vec![0u8; 3]], k).err(), Some(Error::InvalidParameters));
    }

    #[test]
    fn round1_and_round2_reject_bad_inputs() {
        let params = get_threshold_params(2, 2).unwrap();
        let (_, sks) = keygen_from_seed(&[2u8; 32], &params).unwrap();
        let mut rng = StdRng::seed_from_u64(2);
        let mut guard = InMemoryNonceGuard::new();
        // act does not contain the party
        assert_eq!(
            round1(&sks[0], &params, 0b10, MSG, &SID, &mut rng, &mut guard).err(),
            Some(Error::InvalidParameters)
        );
        // key for other (T, N)
        let other = get_threshold_params(2, 3).unwrap();
        assert_eq!(
            round1(&sks[0], &other, 0b11, MSG, &SID, &mut rng, &mut guard).err(),
            Some(Error::InvalidParameters)
        );

        // round2: own slot and own hash agree, but neither matches our w
        let (h, st1) = round1(&sks[0], &params, ACT, MSG, &SID, &mut rng, &mut guard).unwrap();
        let fake = [9u8; 32];
        assert_eq!(
            round2(&sks[0], ACT, MSG, &SID, &[fake, h], &fake, st1, &params).err(),
            Some(Error::InvalidShare)
        );
    }

    #[test]
    fn round3_rejects_wrong_key_params() {
        let f = fixture();
        let (_, other) = keygen_from_seed(&[3u8; 32], &get_threshold_params(2, 3).unwrap()).unwrap();
        assert_eq!(
            round3(&other[0], &f.wfinals, f.st2, &f.params, f.verified).err(),
            Some(Error::InvalidParameters)
        );
    }

    #[test]
    fn round3_rejects_short_wfinals() {
        let f = fixture();
        assert_eq!(
            round3(&f.sks[0], &f.wfinals[..1], f.st2, &f.params, f.verified).err(),
            Some(Error::InvalidParameters)
        );
    }

    #[test]
    fn round3_rejects_witness_with_wrong_k() {
        let mut f = fixture();
        f.verified.k_reps += 1;
        assert_eq!(
            round3(&f.sks[0], &f.wfinals, f.st2, &f.params, f.verified).err(),
            Some(Error::InvalidParameters)
        );
    }

    #[test]
    fn round3_rejects_witness_from_other_session() {
        let mut f = fixture();
        f.verified.session_id[0] ^= 1;
        assert_eq!(
            round3(&f.sks[0], &f.wfinals, f.st2, &f.params, f.verified).err(),
            Some(Error::InvalidShare)
        );
    }

    #[test]
    fn round3_rejects_substituted_aggregate() {
        let mut f = fixture();
        f.wfinals[0].polys[0].coeffs[0] = (f.wfinals[0].polys[0].coeffs[0] + 1) % Q;
        assert_eq!(
            round3(&f.sks[0], &f.wfinals, f.st2, &f.params, f.verified).err(),
            Some(Error::InvalidShare)
        );
    }

    #[test]
    fn round3_rejects_forged_binding() {
        let mut f = fixture();
        f.verified.binding[0] ^= 1;
        assert_eq!(
            round3(&f.sks[0], &f.wfinals, f.st2, &f.params, f.verified).err(),
            Some(Error::InvalidShare)
        );
    }

    #[test]
    fn combine_rejects_responses_for_another_message() {
        // Real commitments and responses, but combined for a different
        // message: every slot packs and then fails FIPS 204 verification.
        let params = get_threshold_params(2, 2).unwrap();
        let (pk, sks) = keygen_from_seed(&[5u8; 32], &params).unwrap();
        let mut rng = StdRng::seed_from_u64(5);
        let mut guard = InMemoryNonceGuard::new();
        let k = params.k_reps as usize;
        let size = pack_w_single_size();
        let mut hashes = Vec::new();
        let mut st1s = Vec::new();
        for sk in &sks {
            let (h, st1) = round1(sk, &params, ACT, MSG, &SID, &mut rng, &mut guard).unwrap();
            hashes.push(h);
            st1s.push(st1);
        }
        let mut reveals = Vec::new();
        let mut st2s = Vec::new();
        for (i, st1) in st1s.into_iter().enumerate() {
            let (r, st2) = round2(&sks[i], ACT, MSG, &SID, &hashes, &hashes[i], st1, &params).unwrap();
            reveals.push(r);
            st2s.push(st2);
        }
        let parsed: Vec<Vec<PolyVecK>> = reveals
            .iter()
            .map(|r| (0..k).map(|j| unpack_w_single(&r[j * size..(j + 1) * size])).collect())
            .collect();
        let wfinals = coordinator::aggregate_commitments(&parsed, k).unwrap();
        let mut responses = Vec::new();
        for (i, st2) in st2s.into_iter().enumerate() {
            let v = verify_all_round2_reveals(&sks[i].tr, &[0, 1], ACT, MSG, &SID, &reveals, &hashes, k).unwrap();
            responses.push(round3(&sks[i], &wfinals, st2, &params, v).unwrap());
        }
        let zfinals = coordinator::aggregate_responses(&responses, k).unwrap();
        assert_eq!(
            coordinator::combine(&pk, b"another message", &wfinals, &zfinals, &params).err(),
            Some(Error::InsufficientResponses)
        );
    }

    #[test]
    fn pack_w_single_maps_negative_coefficients() {
        let mut w = PolyVecK::zero();
        w.polys[0].coeffs[0] = -1;
        let back = unpack_w_single(&pack_w_single(&w));
        assert_eq!(back.polys[0].coeffs[0], Q - 1);
    }
}

#[cfg(test)]
mod transcript_tests {
    use super::*;
    use crate::nonce::InMemoryNonceGuard;
    use crate::params::get_threshold_params;
    use crate::rss::keygen_from_seed;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn round2_needs_one_hash_per_signer_and_round3_canonical_aggregates() {
        let params = get_threshold_params(2, 2).unwrap();
        let (_, sks) = keygen_from_seed(&[19u8; 32], &params).unwrap();
        let mut rng = StdRng::seed_from_u64(19);
        let mut guard = InMemoryNonceGuard::new();
        let sid = [1u8; 32];
        let (h, st1) = round1(&sks[0], &params, 0b11, b"m", &sid, &mut rng, &mut guard).unwrap();
        assert_eq!(
            round2(&sks[0], 0b11, b"m", &sid, &[h], &h, st1, &params).err(),
            Some(Error::InvalidParameters)
        );

        let mut hashes = Vec::new();
        let mut st1s = Vec::new();
        let sid = [2u8; 32];
        for sk in &sks {
            let (h, st1) = round1(sk, &params, 0b11, b"m", &sid, &mut rng, &mut guard).unwrap();
            hashes.push(h);
            st1s.push(st1);
        }
        let mut reveals = Vec::new();
        let mut st2s = Vec::new();
        for (i, st1) in st1s.into_iter().enumerate() {
            let (r, st2) = round2(&sks[i], 0b11, b"m", &sid, &hashes, &hashes[i], st1, &params).unwrap();
            reveals.push(r);
            st2s.push(st2);
        }
        let k = params.k_reps as usize;
        let verified = verify_all_round2_reveals(&sks[0].tr, &[0, 1], 0b11, b"m", &sid, &reveals, &hashes, k).unwrap();
        let mut wfinals = alloc::vec![PolyVecK::zero(); k];
        wfinals[0].polys[0].coeffs[0] = -1;
        let st2 = st2s.into_iter().next().unwrap();
        assert_eq!(round3(&sks[0], &wfinals, st2, &params, verified).err(), Some(Error::InvalidParameters));
    }
}

#[cfg(test)]
mod binding_tests {
    use super::*;
    use crate::coordinator;
    use crate::nonce::InMemoryNonceGuard;
    use crate::params::get_threshold_params;
    use crate::rss::keygen_from_seed;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn round1_requires_exactly_t_signers_and_round2_matching_key() {
        let params = get_threshold_params(2, 3).unwrap();
        let (_, sks) = keygen_from_seed(&[40u8; 32], &params).unwrap();
        let mut rng = StdRng::seed_from_u64(40);
        let mut guard = InMemoryNonceGuard::new();
        assert_eq!(
            round1(&sks[0], &params, 0b111, b"m", &[1; 32], &mut rng, &mut guard).err(),
            Some(Error::InvalidParameters)
        );
        // No guard counter is spent on a rejected call.
        assert_eq!(guard.counter(&sks[0].guard_key_id()), 0);

        let (h, st1) = round1(&sks[0], &params, 0b011, b"m", &[2; 32], &mut rng, &mut guard).unwrap();
        let other = get_threshold_params(2, 2).unwrap();
        assert_eq!(
            round2(&sks[0], 0b011, b"m", &[2; 32], &[h, h], &h, st1, &other).err(),
            Some(Error::InvalidParameters)
        );
    }

    #[test]
    fn partial_secret_requires_every_assigned_share() {
        let params = get_threshold_params(2, 3).unwrap();
        let (_, sks) = keygen_from_seed(&[43u8; 32], &params).unwrap();
        // Party 0 does not hold subset {1, 2}.
        assert_eq!(recover_partial_secret(&sks[0], &[0b011, 0b110]).err(), Some(Error::InvalidShare));
        assert!(recover_partial_secret(&sks[0], &[0b011, 0b101]).is_ok());
    }

    #[test]
    fn round3_refuses_a_nonce_drawn_by_another_key() {
        // Two independent (2, 2) key sets with party 0 swapped in round 3.
        let params = get_threshold_params(2, 2).unwrap();
        let (_, sks) = keygen_from_seed(&[41u8; 32], &params).unwrap();
        let (_, alien) = keygen_from_seed(&[42u8; 32], &params).unwrap();
        let mut rng = StdRng::seed_from_u64(41);
        let mut guard = InMemoryNonceGuard::new();
        let sid = [3u8; 32];
        let mut hashes = Vec::new();
        let mut st1s = Vec::new();
        for sk in &sks {
            let (h, st1) = round1(sk, &params, 0b11, b"m", &sid, &mut rng, &mut guard).unwrap();
            hashes.push(h);
            st1s.push(st1);
        }
        let mut reveals = Vec::new();
        let mut st2s = Vec::new();
        for (i, st1) in st1s.into_iter().enumerate() {
            let (r, st2) = round2(&sks[i], 0b11, b"m", &sid, &hashes, &hashes[i], st1, &params).unwrap();
            reveals.push(r);
            st2s.push(st2);
        }
        let k = params.k_reps as usize;
        let size = pack_w_single_size();
        let parsed: Vec<Vec<PolyVecK>> = reveals
            .iter()
            .map(|r| (0..k).map(|j| unpack_w_single(&r[j * size..(j + 1) * size])).collect())
            .collect();
        let wfinals = coordinator::aggregate_commitments(&parsed, k).unwrap();
        let verified = verify_all_round2_reveals(&sks[0].tr, &[0, 1], 0b11, b"m", &sid, &reveals, &hashes, k).unwrap();
        let st2 = st2s.into_iter().next().unwrap();
        assert_eq!(round3(&alien[0], &wfinals, st2, &params, verified).err(), Some(Error::InvalidParameters));
    }
}
