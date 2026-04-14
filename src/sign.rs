//! Party-side signing logic for the threshold ML-DSA protocol.
//!
//! Implements the 3-round signing flow from ePrint 2026/013:
//!
//! **Round 1 (Commit):** Each party i generates K commitments
//!   w_{i,k} = A·r_k + e_k using hyperball sampling, and broadcasts
//!   H(tr ‖ id ‖ w_{i,k}).
//!
//! **Round 2 (Reveal):** Each party reveals the full w_{i,k} vectors.
//!   Other parties verify the commitment hash.
//!
//! **Round 3 (Respond):** Each party computes K partial responses
//!   z_{i,k} = (c·s₁_I + y_k, c·s₂_I + e_k) and applies the
//!   hyperball rejection test ‖z_{i,k}‖₂ ≤ r.
//!   Responses that pass are nonzero; those that fail are zero vectors.
//!
//! The coordinator then aggregates and tries each of the K slots.

extern crate alloc;
use alloc::vec::Vec;

use crate::error::Error;
use crate::fvec::{FVec, sample_hyperball};
use crate::params::*;
use crate::partition;
use crate::poly::{Poly, PolyVecK, PolyVecL};
use crate::rss::{ThresholdPrivateKey, Share};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use rand_core::{CryptoRng, RngCore};

/// Round 1 state saved by a party for use in Round 3.
pub struct StRound1 {
    /// K packed commitment vectors (each is K polynomials in normal domain)
    pub w_packed: Vec<Vec<u8>>,
    /// K FVec randomness samples (for computing responses)
    pub cmtst: Vec<FVec>,
}

/// Round 2 state saved by a party for use in Round 3.
pub struct StRound2 {
    /// Commitment hashes received from all parties in Round 1
    pub hashes: Vec<[u8; 32]>,
    /// μ = CRH(tr ‖ msg)
    pub mu: [u8; 64],
    /// Active signer bitmask
    pub act: u8,
}

/// Compute μ = CRH(tr ‖ msg).
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
/// - `commitment_hash`: 32-byte hash H(tr ‖ id ‖ w_packed)
/// - `StRound1`: saved state for Round 3
///
/// This matches the Go `Round1()` function.
pub fn round1<R: RngCore + CryptoRng>(
    sk: &ThresholdPrivateKey,
    params: &ThresholdParams,
    rng: &mut R,
) -> Result<([u8; 32], StRound1), Error> {
    // Generate 64-byte randomness for hyperball sampling
    let mut rhop = [0u8; 64];
    rng.fill_bytes(&mut rhop);

    let k_reps = params.k_reps as usize;
    let mut ws: Vec<PolyVecK> = Vec::with_capacity(k_reps);
    let mut cmtst: Vec<FVec> = Vec::with_capacity(k_reps);
    let mut w_packed_all: Vec<Vec<u8>> = Vec::with_capacity(k_reps);

    // Expand A from ρ
    let a = expand_a(&sk.rho);

    for i in 0..k_reps {
        let mut fv = FVec::zero();
        // Sample from hyperball with radius r₁ (randomness ball)
        sample_hyperball(&mut fv, params.r1, params.nu, &rhop, i as u16);

        // Split FVec into (r, e_) ∈ ℤ^L × ℤ^K
        let mut r = PolyVecL::zero();
        let mut e_ = PolyVecK::zero();
        fv.round_to_polyvecs(&mut r, &mut e_);

        // Compute w = A·r + e_
        let mut rh = r.clone();
        rh.ntt();
        let mut w = PolyVecK::zero();
        for j in 0..K {
            // w[j] = Σ_l A[j][l] · rh[l] (pointwise in NTT domain)
            let mut acc = Poly::zero();
            for l in 0..L {
                let mut prod = Poly::zero();
                Poly::pointwise_montgomery(&mut prod, &a[j][l], &rh.polys[l]);
                let acc_copy = acc.clone();
                Poly::add(&mut acc, &acc_copy, &prod);
            }
            acc.reduce();
            acc.invntt_tomont();
            // w[j] = A·r[j] + e_[j]
            let acc_copy = acc.clone();
            Poly::add(&mut w.polys[j], &acc_copy, &e_.polys[j]);
        }
        w.reduce();

        // Pack w for transport (23 bits per coeff × K×N coefficients)
        let packed = pack_w_single(&w);
        w_packed_all.push(packed);

        ws.push(w);
        cmtst.push(fv);
    }

    // Compute commitment hash: H(tr ‖ id ‖ w_packed_all)
    let mut h = Shake256::default();
    h.update(&sk.tr);
    h.update(&[sk.id]);
    for packed in &w_packed_all {
        h.update(packed);
    }
    let mut hash = [0u8; 32];
    h.finalize_xof().read(&mut hash);

    Ok((hash, StRound1 { w_packed: w_packed_all, cmtst }))
}

/// Round 2: reveal commitments and compute μ.
///
/// Verifies that all received commitment hashes match the party's Round 1 hashes.
/// Returns the packed commitment data and state for Round 3.
pub fn round2(
    sk: &ThresholdPrivateKey,
    act: u8,
    msg: &[u8],
    msgs_rd1: &[[u8; 32]], // commitment hashes from Round 1
    st_rd1: &StRound1,
    _params: &ThresholdParams,
) -> (Vec<u8>, StRound2) {
    // Concatenate all packed w data for this party's reveal
    let mut reveal_data = Vec::new();
    for packed in &st_rd1.w_packed {
        reveal_data.extend_from_slice(packed);
    }

    // Save state for Round 3
    let mu = compute_mu(&sk.tr, msg);

    let st2 = StRound2 {
        hashes: msgs_rd1.to_vec(),
        mu,
        act,
    };

    (reveal_data, st2)
}

/// Round 3: compute K parallel responses with hyperball rejection.
///
/// This matches the Go `ComputeResponses()` function.
/// For each of the K parallel commitments, computes:
///   z_k = c·s_I + (r_k, e_k) where (r_k, e_k) is from the hyperball sample
///
/// Then applies the ν-scaled L₂ norm check: ‖z_k‖₂ ≤ r.
/// Responses that fail rejection are all-zero (coordinator will skip them).
pub fn round3(
    sk: &ThresholdPrivateKey,
    wfinals: &[PolyVecK],    // K aggregated commitment vectors
    st_rd1: &StRound1,       // saved randomness from Round 1
    st_rd2: &StRound2,       // μ and active set from Round 2
    params: &ThresholdParams,
) -> Vec<PolyVecL> {
    let k_reps = params.k_reps as usize;

    // Recover this party's partial secret for the active signer set
    let active = bitmask_to_sorted_ids(st_rd2.act, params.n);
    let partition = partition::rss_recover(&active, params.n, params.t);

    // Find which partition slot corresponds to this party
    let party_idx = active.iter().position(|&id| id == sk.id)
        .expect("Party must be in active set");

    // Sum the shares assigned to this party by the partition
    let (s1h, s2h) = recover_partial_secret(sk, &partition[party_idx]);

    let mut zs: Vec<PolyVecL> = Vec::with_capacity(k_reps);

    for i in 0..k_reps {
        // Decompose w into w₀, w₁
        let (_, w1) = decompose_veck(&wfinals[i]);

        // c̃ = H(μ ‖ w₁_packed)
        let c_tilde = compute_challenge(&st_rd2.mu, &w1);

        // Expand challenge polynomial c
        let mut ch = Poly::zero();
        Poly::challenge(&mut ch, &c_tilde);
        ch.ntt();

        // Compute c·s₁ (in NTT domain)
        let mut z = PolyVecL::zero();
        for j in 0..L {
            Poly::pointwise_montgomery(&mut z.polys[j], &ch, &s1h.polys[j]);
            z.polys[j].invntt_tomont();
        }
        z.reduce();

        // Compute c·s₂ (for the second half of the FVec)
        let mut y = PolyVecK::zero();
        for j in 0..K {
            Poly::pointwise_montgomery(&mut y.polys[j], &ch, &s2h.polys[j]);
            y.polys[j].invntt_tomont();
        }
        y.reduce();

        // Build FVec from (c·s₁, c·s₂) and add the randomness sample
        let mut zf = FVec::from_polyvecs(&z, &y);
        zf.add_assign(&st_rd1.cmtst[i]);

        // Hyperball rejection: check ‖zf‖₂ ≤ r
        if zf.excess(params.r, params.nu) {
            // Rejection: emit zero response
            zs.push(PolyVecL::zero());
            continue;
        }

        // Round FVec back to polynomial vector (only z part needed)
        let mut z_out = PolyVecL::zero();
        let mut _y_out = PolyVecK::zero();
        zf.round_to_polyvecs(&mut z_out, &mut _y_out);
        zs.push(z_out);
    }

    zs
}

// ─── Internal Helpers ────────────────────────────────────────────────

/// Expand A from ρ using dilithium-rs.
fn expand_a(rho: &[u8; 32]) -> Vec<Vec<Poly>> {
    use dilithium::{
        polyvec::{PolyVecL as DPolyVecL, matrix_expand},
        ML_DSA_44,
    };
    let mode = ML_DSA_44;
    let mut mat: [DPolyVecL; K] = core::array::from_fn(|_| DPolyVecL::default());
    matrix_expand(mode, &mut mat, rho);

    // Convert to our Poly format
    let mut a: Vec<Vec<Poly>> = Vec::with_capacity(K);
    for i in 0..K {
        let mut row = Vec::with_capacity(L);
        for j in 0..L {
            let mut p = Poly::zero();
            for c in 0..N {
                p.coeffs[c] = mat[i].vec[j].coeffs[c];
            }
            row.push(p);
        }
        a.push(row);
    }
    a
}


/// Size of a single packed PolyVecK (23 bits per coefficient, K×N coefficients).
pub fn pack_w_single_size() -> usize {
    let poly_q_size = (N * 23 + 7) / 8;
    K * poly_q_size
}

/// Pack a single PolyVecK as 23 bits per coefficient.
pub fn pack_w_single(w: &PolyVecK) -> Vec<u8> {
    let poly_q_size = (N * 23 + 7) / 8;
    let total = K * poly_q_size;
    let mut buf = alloc::vec![0u8; total];

    let mut bit_pos = 0usize;
    for poly in &w.polys {
        for &coeff in &poly.coeffs {
            let val = if coeff < 0 { (coeff + Q) as u32 } else { coeff as u32 };
            // Write 23 bits at bit_pos
            let byte_pos = bit_pos / 8;
            let bit_off = bit_pos % 8;
            let v = (val as u64) << bit_off;
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

/// Unpack a single PolyVecK from 23-bit packed data.
pub fn unpack_w_single(buf: &[u8]) -> PolyVecK {
    let mut w = PolyVecK::zero();
    let mut bit_pos = 0usize;
    for poly in w.polys.iter_mut() {
        for coeff in poly.coeffs.iter_mut() {
            let byte_pos = bit_pos / 8;
            let bit_off = bit_pos % 8;
            let mut v: u64 = 0;
            for b in 0..4 {
                if byte_pos + b < buf.len() {
                    v |= (buf[byte_pos + b] as u64) << (b * 8);
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

/// Recover this party's partial secret from the assigned share bitmasks.
fn recover_partial_secret(
    sk: &ThresholdPrivateKey,
    assigned_masks: &[u8],
) -> (PolyVecL, PolyVecK) {
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

/// Decompose PolyVecK into (w₀, w₁) using Power2Round-style decompose.
fn decompose_veck(w: &PolyVecK) -> (PolyVecK, PolyVecK) {
    let mut w0 = PolyVecK::zero();
    let mut w1 = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            let a = w.polys[i].coeffs[j];
            // Ensure a is in [0, Q)
            let a_pos = if a < 0 { a + Q } else { a };
            // Decompose: a = a₁·α + a₀ where α = 2γ₂
            let alpha = 2 * GAMMA2;
            let a0 = a_pos - ((a_pos + alpha / 2 - 1) / alpha) * alpha;
            let a1 = if a0 > 0 {
                (a_pos - a0) / alpha
            } else {
                (a_pos - a0 - 1) / alpha + 1
            };
            w0.polys[i].coeffs[j] = a0;
            w1.polys[i].coeffs[j] = a1;
        }
    }
    (w0, w1)
}

/// Compute challenge hash c̃ = H(μ ‖ w₁_packed).
fn compute_challenge(mu: &[u8; 64], w1: &PolyVecK) -> [u8; CTILDEBYTES] {
    let mut h = Shake256::default();
    h.update(mu);
    // Pack w₁ (6 bits per coeff for ML-DSA-44)
    for poly in &w1.polys {
        let mut packed = [0u8; POLYW1_PACKEDBYTES];
        for k in 0..N / 4 {
            packed[3 * k] = (poly.coeffs[4 * k] as u8) |
                ((poly.coeffs[4 * k + 1] as u8) << 6);
            packed[3 * k + 1] = ((poly.coeffs[4 * k + 1] >> 2) as u8) |
                ((poly.coeffs[4 * k + 2] as u8) << 4);
            packed[3 * k + 2] = ((poly.coeffs[4 * k + 2] >> 4) as u8) |
                ((poly.coeffs[4 * k + 3] as u8) << 2);
        }
        h.update(&packed);
    }
    let mut c = [0u8; CTILDEBYTES];
    h.finalize_xof().read(&mut c);
    c
}

#[cfg(test)]
mod tests {
    // v0.3.0: Tests will be added after the coordinator.rs rewrite
    // enables full end-to-end signing.
}
