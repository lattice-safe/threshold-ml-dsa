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
use crate::params::{TRBYTES, ThresholdParams, K, L, N, Q, CTILDEBYTES};
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
pub struct StRound1 {
    /// K packed commitment vectors (each is K polynomials in normal domain)
    pub w_packed: Vec<Vec<u8>>,
    /// K randomness samples as floating vectors from `SampleHyperball`.
    ///
    /// Round 1 commitments use rounded copies of these vectors.
    /// Round 3 adds them in floating-point before the final rounding step.
    pub rand_fvecs: Vec<FVec>,
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
/// Contains commitment hashes and message hash; zeroized on drop.
#[derive(Zeroize)]
#[zeroize(drop)]
pub struct StRound2 {
    /// Commitment hashes received from all parties in Round 1
    pub hashes: Vec<[u8; 32]>,
    /// μ = CRH(tr ‖ msg)
    pub mu: [u8; 64],
    /// Active signer bitmask
    pub act: u8,
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
/// Returns:
/// - `commitment_hash`: 32-byte hash H(tag ‖ tr ‖ id ‖ act ‖ session ‖ μ ‖ `w_packed`)
/// - `StRound1`: saved state for Round 3
///
/// This matches the Go `Round1()` function.
pub fn round1<R: RngCore + CryptoRng>(
    sk: &ThresholdPrivateKey,
    params: &ThresholdParams,
    act: u8,
    msg: &[u8],
    session_id: &[u8; 32],
    rng: &mut R,
) -> Result<([u8; 32], StRound1), Error> {
    let mu = compute_mu(&sk.tr, msg);

    // Hedged nonce derivation:
    // bind RNG entropy to long-term per-party key and transcript context so
    // RNG rollback does not trivially re-use nonces across distinct sessions.
    let mut rng_entropy = [0u8; 64];
    rng.fill_bytes(&mut rng_entropy);
    let mut rhop = [0u8; 64];
    let mut rhop_h = Shake256::default();
    rhop_h.update(b"th-ml-dsa-round1-rhop-v1");
    rhop_h.update(&rng_entropy);
    rhop_h.update(&sk.key);
    rhop_h.update(&sk.tr);
    rhop_h.update(&[sk.id, act]);
    rhop_h.update(session_id);
    rhop_h.update(&mu);
    rhop_h.finalize_xof().read(&mut rhop);
    rng_entropy.zeroize();

    let k_reps = params.k_reps as usize;
    let mut rand_fvecs: Vec<FVec> = Vec::with_capacity(k_reps);
    let mut w_packed_all: Vec<Vec<u8>> = Vec::with_capacity(k_reps);

    // Expand A from ρ
    let a = expand_a(&sk.rho);

    for i in 0..k_reps {
        let mut fv = FVec::zero();
        // Sample from hyperball with radius r₁ (randomness ball)
        sample_hyperball(&mut fv, params.r1, params.nu, &rhop, i as u16);

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
            let mut acc = Poly::zero();
            for (l, a_jl) in a_row.iter().enumerate().take(L) {
                let mut prod = Poly::zero();
                Poly::pointwise_montgomery(&mut prod, a_jl, &yh.polys[l]);
                let acc_copy = acc.clone();
                Poly::add(&mut acc, &acc_copy, &prod);
            }
            acc.reduce();
            acc.invntt_tomont();
            // w[j] = A·y[j] + e_[j]
            let acc_copy = acc.clone();
            Poly::add(&mut w.polys[j], &acc_copy, &e_.polys[j]);
        }
        w.reduce();

        // Pack w for transport (23 bits per coeff × K×N coefficients)
        let packed = pack_w_single(&w);
        w_packed_all.push(packed);

        rand_fvecs.push(fv);
    }

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
        },
    ))
}

/// Round 2: reveal commitments and compute μ.
///
/// Verifies that all received commitment hashes match the party's Round 1 hashes.
/// Returns the packed commitment data and state for Round 3.
#[allow(clippy::too_many_arguments)]
pub fn round2(
    sk: &ThresholdPrivateKey,
    act: u8,
    msg: &[u8],
    session_id: &[u8; 32],
    msgs_rd1: &[[u8; 32]], // commitment hashes from Round 1
    own_rd1_hash: &[u8; 32],
    st_rd1: &StRound1,
    _params: &ThresholdParams,
) -> Result<(Vec<u8>, StRound2), Error> {
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
pub fn round3(
    sk: &ThresholdPrivateKey,
    wfinals: &[PolyVecK], // K aggregated commitment vectors
    st_rd1: &StRound1,    // saved randomness from Round 1
    st_rd2: &StRound2,    // μ and active set from Round 2
    params: &ThresholdParams,
) -> Result<Vec<PolyVecL>, Error> {
    use dilithium::{
        poly::Poly as DPoly,
        polyvec::{polyveck_reduce, polyvecl_reduce, PolyVecK as DPolyVecK, PolyVecL as DPolyVecL},
        ML_DSA_44,
    };

    let k_reps = params.k_reps as usize;
    let mode = ML_DSA_44;

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
    let (s1h, s2h) = recover_partial_secret(sk, &partition[party_idx]);
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
        zf.add_assign(&st_rd1.rand_fvecs[i]);
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
fn recover_partial_secret(sk: &ThresholdPrivateKey, assigned_masks: &[u8]) -> (PolyVecL, PolyVecK) {
    let mut s1h = PolyVecL::zero();
    let mut s2h = PolyVecK::zero();

    for &mask in assigned_masks {
        if let Some(share) = sk.shares.get(&mask) {
            s1h.add_assign(&share.s1h);
            s2h.add_assign(&share.s2h);
        }
    }
    s1h.reduce();
    s2h.reduce();
    (s1h, s2h)
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
