//! Coordinator / Aggregator for threshold ML-DSA signing.
//!
//! The coordinator orchestrates the hardened 4-round distributed flow:
//! 0. pre-commitment hash collection
//! 1. commitment reveal and binding verification
//! 2. Fiat-Shamir challenge derivation
//! 3. partial response aggregation and final signature validation

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use crate::error::Error;
use crate::params::*;
use crate::poly::*;
use crate::rss;
use crate::sign::{verify_commitment, verify_session_binding};
use crate::sign::{Commitment, CommitmentHash, PartialSignature, Party};
use crate::verify;
use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use zeroize::Zeroize;

/// The final FIPS 204-compatible ML-DSA signature.
///
/// Encoded as: c̃ ‖ z ‖ h
/// This can be verified by any unmodified ML-DSA verifier.
#[derive(Clone, Debug)]
pub struct Signature {
    /// The challenge hash c̃ (CTILDEBYTES bytes).
    pub c_tilde: Vec<u8>,
    /// The aggregated response vector z ∈ ℤ_q^L.
    pub z: PolyVecL,
    /// The hint vector h (packed as OMEGA + K bytes).
    pub h: Vec<u8>,
}

impl Signature {
    /// Encode the signature into the standard FIPS 204 byte format.
    ///
    /// Layout: c̃ (32 bytes) ‖ z (L × POLYZ_PACKEDBYTES) ‖ h (OMEGA + K bytes)
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut sig = vec![0u8; SIG_BYTES];

        // Pack c̃
        sig[..CTILDEBYTES].copy_from_slice(&self.c_tilde);

        // Pack z
        let mut offset = CTILDEBYTES;
        for j in 0..L {
            self.z.polys[j].pack_z(&mut sig[offset..offset + POLYZ_PACKEDBYTES]);
            offset += POLYZ_PACKEDBYTES;
        }

        // Pack h (already encoded)
        sig[offset..offset + self.h.len()].copy_from_slice(&self.h);

        sig
    }
}

/// Aggregate partial commitments from all parties.
///
/// Computes w = Σ w_i (component-wise polynomial addition).
/// Verifies each commitment against its pre-committed hash (ADV-1).
///
/// # Security Note
/// This is a low-level building block. It does not by itself enforce
/// transcript consistency across protocol rounds or final signature validity.
/// Prefer [`threshold_sign`] for production use.
pub fn aggregate_commitments(
    commitments: &[Commitment],
    hashes: &[CommitmentHash],
) -> Result<PolyVecK, Error> {
    // ADV-1: Verify all commitment bindings
    if commitments.len() != hashes.len() {
        return Err(Error::InvalidParameters);
    }
    for (comm, hash) in commitments.iter().zip(hashes.iter()) {
        if !verify_commitment(comm, hash) {
            return Err(Error::InvalidShare);
        }
    }

    let mut w = PolyVecK::zero();
    for comm in commitments {
        w.add_assign(&comm.w);
    }
    w.reduce();
    w.caddq();
    Ok(w)
}

/// Compute the Fiat-Shamir challenge hash from the aggregate commitment and message.
///
/// c̃ = H(μ) where μ = SHAKE-256(tr ‖ msg)
/// c̃ = SHAKE-256(μ ‖ w1_packed)
///
/// In the threshold setting, `tr` is derived from the public key.
///
/// # Security Note
/// This is a low-level helper. Callers must ensure the commitment set and
/// signer set are fixed and transcript-consistent.
pub fn compute_challenge(w: &PolyVecK, msg: &[u8], tr: &[u8; TRBYTES]) -> Vec<u8> {
    // Use reference decompose/pack path from dilithium-rs to avoid transcript drift.
    let mode = dilithium::params::ML_DSA_44;
    let mut w_ref = dilithium::polyvec::PolyVecK::default();
    for i in 0..K {
        w_ref.vec[i].coeffs = w.polys[i].coeffs;
    }
    let mut w1_high = dilithium::polyvec::PolyVecK::default();
    let mut w0 = dilithium::polyvec::PolyVecK::default();
    dilithium::polyvec::polyveck_decompose(mode, &mut w1_high, &mut w0, &w_ref);
    let mut w1_packed = vec![0u8; mode.k() * mode.polyw1_packedbytes()];
    dilithium::polyvec::polyveck_pack_w1(mode, &mut w1_packed, &w1_high);

    // FIPS 204 Algorithm 7 (ML-DSA.Sign)
    // μ = SHAKE-256(tr ‖ msg, 64)
    let mut hasher_mu = Shake256::default();
    hasher_mu.update(tr);
    hasher_mu.update(msg);
    let mut reader_mu = hasher_mu.finalize_xof();
    let mut mu = [0u8; 64];
    reader_mu.read(&mut mu);

    // c̃ = SHAKE-256(μ ‖ w1_packed, 32)
    let mut hasher = Shake256::default();
    hasher.update(&mu);
    hasher.update(&w1_packed);
    let mut reader = hasher.finalize_xof();
    let mut c_tilde = vec![0u8; CTILDEBYTES];
    reader.read(&mut c_tilde);

    // SECURITY: zeroize intermediate hash
    mu.zeroize();

    c_tilde
}

/// Build a per-attempt session context used by parties for nonce hedging.
fn build_session_context(
    tr: &[u8; TRBYTES],
    msg: &[u8],
    active_party_ids: &[usize],
    attempt: usize,
) -> [u8; 32] {
    let mut hasher = Shake256::default();
    hasher.update(tr);
    hasher.update(msg);
    hasher.update(&(attempt as u64).to_le_bytes());
    hasher.update(&(active_party_ids.len() as u64).to_le_bytes());
    for &id in active_party_ids {
        hasher.update(&(id as u64).to_le_bytes());
    }
    let mut out = [0u8; 32];
    hasher.finalize_xof().read(&mut out);
    out
}

/// Check that the active signer set intersects every RSS subset.
fn active_set_covers_all_subsets(n: usize, t: usize, active_party_ids: &[usize]) -> bool {
    let mut active = [false; MAX_PARTIES];
    for &id in active_party_ids {
        if id >= MAX_PARTIES {
            return false;
        }
        active[id] = true;
    }

    let subsets = rss::enumerate_subsets(n, t);
    subsets
        .iter()
        .all(|subset| subset.iter().any(|&member| active[member]))
}

/// Enumerate combinations of `k` distinct indices from `[0, n)`.
fn index_combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    fn helper(
        n: usize,
        k: usize,
        start: usize,
        current: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == k {
            out.push(current.clone());
            return;
        }
        let remaining = k - current.len();
        for i in start..=(n - remaining) {
            current.push(i);
            helper(n, k, i + 1, current, out);
            current.pop();
        }
    }

    if k == 0 || k > n {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut current = Vec::with_capacity(k);
    helper(n, k, 0, &mut current, &mut out);
    out
}

/// Aggregate partial responses from parties that passed rejection sampling.
///
/// Computes z = Σ z_i and constructs the final signature.
/// Enforces:
/// - Party deduplication (ADV-2)
/// - Session binding verification (ADV-6)
///
/// # Security Note
/// This is a low-level aggregator and does not perform the final external
/// verification step. Production callers should use [`threshold_sign`], which
/// includes an end-to-end validity gate before returning success.
pub fn aggregate_responses(
    partials: &[PartialSignature],
    w: &PolyVecK,
    c_tilde: &[u8],
    t: usize,
    pk: &[u8],
) -> Result<Signature, Error> {
    if pk.len() != PK_BYTES {
        return Err(Error::InvalidParameters);
    }
    if partials.len() < t {
        return Err(Error::InsufficientResponses);
    }

    // SECURITY (ADV-2): Deduplicate by party_id to prevent Sybil attacks.
    // A corrupt party could submit multiple z_i to bias the aggregate.
    {
        let mut seen_ids = [false; MAX_PARTIES];
        for partial in partials {
            if partial.party_id >= MAX_PARTIES || seen_ids[partial.party_id] {
                return Err(Error::InvalidShare);
            }
            seen_ids[partial.party_id] = true;
        }
    }

    // SECURITY (ADV-6): Verify session binding on all partials.
    for partial in partials {
        if !verify_session_binding(partial, c_tilde) {
            return Err(Error::InvalidShare);
        }
    }

    // Sum all z_i
    let mut z = PolyVecL::zero();
    for partial in partials {
        z.add_assign(&partial.z);
    }
    z.reduce();
    // NOTE: Do NOT call z.caddq() here (ALG-6).
    // z must remain in centered form [≈ -(Q/2), Q/2] for pack_z,
    // which computes GAMMA1 - coeff and assumes coeff ∈ [-(γ₁-1), γ₁].
    // caddq() would map to [0, Q), making GAMMA1 - coeff negative → garbage.

    // Check ‖z‖∞ < γ₁ − β (FIPS 204 requirement on the final signature)
    if z.chknorm(GAMMA1 - BETA) {
        return Err(Error::InvalidSignature);
    }

    // ── FIPS 204 Hint Computation ──
    // Recompute w_approx = A·z - c·t₁·2^d
    // Match reference: do subtraction in NTT domain, single INTT at the end (ALG-3)

    // Extract ρ and t₁ from public key
    let mut rho = [0u8; SEEDBYTES];
    rho.copy_from_slice(&pk[0..SEEDBYTES]);
    let mut t1 = PolyVecK::zero();
    for i in 0..K {
        let start = SEEDBYTES + i * POLYT1_PACKEDBYTES;
        Poly::unpack_t1(&mut t1.polys[i], &pk[start..start + POLYT1_PACKEDBYTES]);
    }

    // A·z in NTT domain
    let a_hat = MatrixA::expand(&rho);
    let mut z_hat = z.clone();
    z_hat.ntt();
    let mut w_ntt = a_hat.mul_vec(&z_hat); // w_ntt = A·z (NTT domain)

    // c·t₁·2^d in NTT domain
    let mut c_poly = Poly::zero();
    Poly::challenge(&mut c_poly, c_tilde);
    c_poly.ntt();

    // ALG-2: shiftl without reduction, matching reference exactly
    for i in 0..K {
        for j in 0..N {
            t1.polys[i].coeffs[j] <<= D;
        }
    }
    t1.ntt();

    // Subtract c·t₁·2^d from Az in NTT domain
    for i in 0..K {
        let mut ct1_i = Poly::zero();
        Poly::pointwise_montgomery(&mut ct1_i, &c_poly, &t1.polys[i]);
        for j in 0..N {
            w_ntt.polys[i].coeffs[j] -= ct1_i.coeffs[j];
        }
    }

    // Single INTT at the end (ALG-3)
    w_ntt.reduce();
    w_ntt.invntt_tomont();
    w_ntt.caddq();
    let w_approx = w_ntt;

    // Generate FIPS Hint Encoding
    let mut h = vec![0u8; OMEGA + K];
    let mut hint_count = 0;

    for i in 0..K {
        for j in 0..N {
            let (w1_expected, _) = decompose(w.polys[i].coeffs[j]);
            let (w1_approx, _) = decompose(w_approx.polys[i].coeffs[j]);
            if w1_expected != w1_approx {
                if hint_count >= OMEGA {
                    return Err(Error::InvalidSignature);
                }
                h[hint_count] = j as u8;
                hint_count += 1;
            }
        }
        h[OMEGA + i] = hint_count as u8;
    }

    Ok(Signature {
        c_tilde: c_tilde.to_vec(),
        z,
        h,
    })
}

/// Produce a threshold signature from qualifying parties using the
/// distributed commit/challenge/response flow.
///
/// # Arguments
/// * `parties` — mutable slice of available signing parties (may be a subset)
/// * `msg` — the message to sign
/// * `tr` — the public key hash (tr = H(pk))
/// * `pk` — the public key bytes
/// * `t` — the threshold
/// * `max_retries` — maximum number of retry attempts on local rejection
/// * `rng` — cryptographically secure RNG
///
/// # Returns
/// A FIPS 204-compatible signature on success.
pub fn threshold_sign<R: RngCore + CryptoRng>(
    parties: &mut [Party],
    msg: &[u8],
    tr: &[u8; TRBYTES],
    pk: &[u8],
    t: usize,
    max_retries: usize,
    rng: &mut R,
) -> Result<Signature, Error> {
    if pk.len() != PK_BYTES
        || parties.is_empty()
        || t < 2
        || t > parties.len()
        || parties.len() > MAX_PARTIES
        || max_retries == 0
    {
        return Err(Error::InvalidParameters);
    }

    // Validate party ID uniqueness and sharing metadata.
    let (share_n, share_t) = parties[0].sharing_params();
    if share_n > MAX_PARTIES || share_t != t || share_t > share_n || share_t < 2 {
        return Err(Error::InvalidParameters);
    }
    let mut seen_ids = [false; MAX_PARTIES];
    for p in parties.iter() {
        if p.id >= share_n || seen_ids[p.id] {
            return Err(Error::InvalidParameters);
        }
        seen_ids[p.id] = true;
        if p.sharing_params() != (share_n, share_t) {
            return Err(Error::InvalidParameters);
        }
        if p.id != p.key_share().party_id
            || p.key_share().n != share_n
            || p.key_share().t != share_t
        {
            return Err(Error::InvalidShare);
        }
    }

    // Build all candidate active signer combinations of size t.
    let mut candidate_sets: Vec<Vec<usize>> = Vec::new();
    for combo in index_combinations(parties.len(), t) {
        let mut ids: Vec<usize> = combo.iter().map(|&idx| parties[idx].id).collect();
        ids.sort_unstable();
        if active_set_covers_all_subsets(share_n, share_t, &ids) {
            candidate_sets.push(combo);
        }
    }
    if candidate_sets.is_empty() {
        return Err(Error::InvalidParameters);
    }

    let mut saw_local_rejection = false;
    let mut saw_invalid_signature = false;

    for attempt in 0..max_retries {
        let active = &candidate_sets[attempt % candidate_sets.len()];
        let mut active_party_ids: Vec<usize> = active.iter().map(|&idx| parties[idx].id).collect();
        active_party_ids.sort_unstable();

        let session_context = build_session_context(tr, msg, &active_party_ids, attempt);

        // Round 0: pre-commit
        let mut hashes = Vec::with_capacity(t);
        let mut precommit_failed = false;
        for &idx in active {
            match parties[idx].precommit(rng, &active_party_ids, &session_context) {
                Ok(hash) => hashes.push(hash),
                Err(_) => {
                    precommit_failed = true;
                    break;
                }
            }
        }
        if precommit_failed || hashes.len() != t {
            continue;
        }

        // Round 1: reveal
        let commitments: Vec<Commitment> = active
            .iter()
            .filter_map(|&idx| parties[idx].reveal())
            .collect();
        if commitments.len() != t {
            continue;
        }

        // Round 2: challenge
        let w = match aggregate_commitments(&commitments, &hashes) {
            Ok(w) => w,
            Err(_) => continue,
        };
        let c_tilde = compute_challenge(&w, msg, tr);

        // Round 3: partial signs from the same active signer set
        let mut partials = Vec::with_capacity(t);
        let mut round_rejected = false;
        for &idx in active {
            match parties[idx].sign(&c_tilde) {
                Ok(partial) => partials.push(partial),
                Err(Error::LocalRejectionAbort) => {
                    saw_local_rejection = true;
                    round_rejected = true;
                }
                Err(e) => return Err(e),
            }
        }
        if round_rejected || partials.len() != t {
            continue;
        }

        match aggregate_responses(&partials, &w, &c_tilde, t, pk) {
            Ok(sig) => {
                let sig_bytes = sig.to_bytes();
                if verify::verify(&sig_bytes, msg, pk) {
                    return Ok(sig);
                }
                saw_invalid_signature = true;
            }
            Err(Error::InvalidSignature) => {
                saw_invalid_signature = true;
            }
            Err(e) => return Err(e),
        }
    }

    if saw_invalid_signature {
        return Err(Error::InvalidSignature);
    }
    if saw_local_rejection {
        return Err(Error::InsufficientResponses);
    }
    Err(Error::InvalidSignature)
}
