//! Coordinator / Aggregator for threshold ML-DSA signing.
//!
//! Implements the `Combine` function from the Go reference.
//! Tries each of the K parallel commitment/response slots until
//! one produces a valid FIPS 204 signature.
//!
//! The coordinator:
//! 1. Aggregates commitments `w_final` = Σ `w_i`
//! 2. Aggregates responses `z_final` = Σ `z_i`
//! 3. For each k ∈ [0, K): checks if `z_k` produces a valid signature

extern crate alloc;
use alloc::vec::Vec;

use crate::error::Error;
use crate::params::{SIG_BYTES, CTILDEBYTES, L, POLYZ_PACKEDBYTES, PK_BYTES, ThresholdParams, TRBYTES, K, GAMMA1, BETA, N, Q, GAMMA2, OMEGA, POLYW1_PACKEDBYTES};
use crate::poly::{PolyVecL, PolyVecK};
use crate::sign;
use crate::verify;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// The final FIPS 204-compatible ML-DSA signature.
///
/// Encoded as: c̃ ‖ z ‖ h
/// This can be verified by any unmodified ML-DSA verifier.
#[derive(Clone, Debug)]
pub struct Signature {
    /// The challenge hash c̃ (CTILDEBYTES bytes).
    pub c_tilde: Vec<u8>,
    /// The aggregated response vector z ∈ `ℤ_q^L`.
    pub z: PolyVecL,
    /// The hint vector h (packed as OMEGA + K bytes).
    pub h: Vec<u8>,
}

impl Signature {
    /// Encode the signature into the standard FIPS 204 byte format.
    #[must_use] 
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut sig = alloc::vec![0u8; SIG_BYTES];
        sig[..CTILDEBYTES].copy_from_slice(&self.c_tilde);
        let mut offset = CTILDEBYTES;
        for j in 0..L {
            self.z.polys[j].pack_z(&mut sig[offset..offset + POLYZ_PACKEDBYTES]);
            offset += POLYZ_PACKEDBYTES;
        }
        sig[offset..offset + self.h.len()].copy_from_slice(&self.h);
        sig
    }
}

/// Aggregate K commitment vectors from T parties.
///
/// For each k ∈ [0, K): wfinal[k] = `Σ_i` w_{i,k}
#[must_use] 
pub fn aggregate_commitments(
    all_reveals: &[Vec<PolyVecK>], // [T parties][K slots]
    k_reps: usize,
) -> Vec<PolyVecK> {
    let mut wfinals: Vec<PolyVecK> = Vec::with_capacity(k_reps);
    for k in 0..k_reps {
        let mut w = PolyVecK::zero();
        for party_reveals in all_reveals {
            w.add_assign(&party_reveals[k]);
        }
        w.reduce();
        wfinals.push(w);
    }
    wfinals
}

/// Aggregate K response vectors from T parties.
///
/// For each k ∈ [0, K): zfinal[k] = `Σ_i` z_{i,k}
#[must_use] 
pub fn aggregate_responses(
    all_responses: &[Vec<PolyVecL>], // [T parties][K slots]
    k_reps: usize,
) -> Vec<PolyVecL> {
    let mut zfinals: Vec<PolyVecL> = Vec::with_capacity(k_reps);
    for k in 0..k_reps {
        let mut z = PolyVecL::zero();
        for party_responses in all_responses {
            z.add_assign(&party_responses[k]);
        }
        z.reduce();
        zfinals.push(z);
    }
    zfinals
}

/// Combine: try each of the K parallel slots to produce a valid signature.
///
/// This matches the Go `Combine()` function.
/// For each slot k:
/// 1. Check ‖`z_k`‖∞ < γ₁ - β
/// 2. Compute `Az_k` - 2^d·c·t₁ = "approximated w"  
/// 3. Check ‖δ‖∞ = ‖`approx_w` - `w_final`‖∞ < γ₂
/// 4. Compute hints
/// 5. If ω hints ≤ OMEGA, pack and return
pub fn combine(
    pk_bytes: &[u8; PK_BYTES],
    msg: &[u8],
    wfinals: &[PolyVecK],
    zfinals: &[PolyVecL],
    params: &ThresholdParams,
) -> Result<[u8; SIG_BYTES], Error> {
    use dilithium::{
        packing::{pack_sig, unpack_pk},
        poly::Poly as DPoly,
        polyvec::{
            matrix_expand, matrix_pointwise_montgomery, PolyVecK as DPolyVecK,
            PolyVecL as DPolyVecL,
        },
        rounding, ML_DSA_44,
    };

    let mode = ML_DSA_44;

    // Unpack public key
    let mut rho = [0u8; 32];
    let mut t1 = DPolyVecK::default();
    unpack_pk(mode, &mut rho, &mut t1, pk_bytes);

    // Compute tr = CRH(pk)
    let mut h_tr = Shake256::default();
    h_tr.update(pk_bytes);
    let mut tr = [0u8; TRBYTES];
    h_tr.finalize_xof().read(&mut tr);

    // Compute μ = CRH(tr ‖ msg)
    let mu = sign::compute_mu(&tr, msg);

    // Expand A
    let mut mat: [DPolyVecL; K] = core::array::from_fn(|_| DPolyVecL::default());
    matrix_expand(mode, &mut mat, &rho);

    let k_reps = params.k_reps as usize;
    for k in 0..k_reps {
        // Check ‖z_k‖∞ < γ₁ - β
        if zfinals[k].chknorm(GAMMA1 - BETA) {
            continue; // exceeds bound
        }

        // Check z is not all-zero (hyperball rejection)
        let all_zero = zfinals[k]
            .polys
            .iter()
            .all(|p| p.coeffs.iter().all(|&c| c == 0));
        if all_zero {
            continue;
        }

        // Compute c̃ = H(μ ‖ w1_packed)
        let (w0, w1) = decompose_polyveck(&wfinals[k]);
        let c_tilde = compute_c_tilde(&mu, &w1);

        // Convert z to dilithium-rs format and NTT
        let mut zh_d = DPolyVecL::default();
        for i in 0..L {
            for j in 0..N {
                zh_d.vec[i].coeffs[j] = zfinals[k].polys[i].coeffs[j];
            }
            zh_d.vec[i].ntt();
        }

        // Compute w_approx = Az - c·t₁·2^d
        // Mirror the FIPS 204 verifier exactly (all in NTT domain first):
        // 1. w_approx = A * NTT(z) (pointwise in NTT domain)
        let mut w_approx = DPolyVecK::default();
        matrix_pointwise_montgomery(mode, &mut w_approx, &mat, &zh_d);

        // 2. cp = NTT(SampleInBall(c̃))
        let mut cp = DPoly::default();
        DPoly::challenge(mode, &mut cp, &c_tilde);
        cp.ntt();

        // 3. ct1 = cp * NTT(t1 << d) in NTT domain
        // 4. w_approx = w_approx - ct1 (still in NTT domain)
        for i in 0..K {
            let mut t1_shifted = t1.vec[i].clone();
            t1_shifted.shiftl();
            t1_shifted.ntt();
            let mut ct1 = DPoly::default();
            DPoly::pointwise_montgomery(&mut ct1, &cp, &t1_shifted);

            let w_copy = w_approx.vec[i].clone();
            DPoly::sub(&mut w_approx.vec[i], &w_copy, &ct1);
        }

        // 5. reduce → invntt → caddq (single pass, matching verifier)
        for i in 0..K {
            w_approx.vec[i].reduce();
            w_approx.vec[i].invntt_tomont();
            w_approx.vec[i].caddq();
        }
        // Now w_approx = Az - c·t₁·2^d in normal domain, [0, Q)

        // Compute hints from:
        //   f = (A*z - c*t1*2^d) - w
        //   w0' = w0 + f
        //   h = MakeHint(w0', w1)
        // with (w0, w1) = Decompose(w).
        //
        // This follows the reference threshold Combine procedure.
        let mut hint = DPolyVecK::default();
        let mut hint_pop = 0usize;
        let mut slot_valid = true;
        for i in 0..K {
            for j in 0..N {
                let target = w1.polys[i].coeffs[j];
                let approx = w_approx.vec[i].coeffs[j];

                // f = approx - w (normalized to [0, Q))
                let mut f = (approx - wfinals[k].polys[i].coeffs[j]).rem_euclid(Q);

                // Reject this slot if ‖f‖∞ > γ2 (using centered representatives).
                if f > (Q - 1) / 2 {
                    f -= Q;
                }
                if f.unsigned_abs() > GAMMA2 as u32 {
                    slot_valid = false;
                    break;
                }

                // Build z0' = w0 + f in [0, Q), then center it for dilithium-rs
                // `make_hint(mode, a0, a1)` which expects signed low bits a0.
                let f_plus_q = if f < 0 { f + Q } else { f };
                let hint_input = (w0.polys[i].coeffs[j] + f_plus_q).rem_euclid(Q);
                let mut hint_low = hint_input;
                if hint_low > (Q - 1) / 2 {
                    hint_low -= Q;
                }
                let hint_bit = rounding::make_hint(mode, hint_low, target);
                if hint_bit {
                    hint.vec[i].coeffs[j] = 1;
                    hint_pop += 1;
                }
            }
            if !slot_valid || hint_pop > OMEGA {
                break;
            }
        }

        if !slot_valid || hint_pop > OMEGA {
            continue;
        }

        // Success! Pack the signature
        let mut z_d = DPolyVecL::default();
        for i in 0..L {
            for j in 0..N {
                let mut centered = zfinals[k].polys[i].coeffs[j].rem_euclid(Q);
                if centered > (Q - 1) / 2 {
                    centered -= Q;
                }
                z_d.vec[i].coeffs[j] = centered;
            }
        }

        let mut sig = [0u8; SIG_BYTES];
        pack_sig(mode, &mut sig, &c_tilde, &z_d, &hint);

        if verify::verify(&sig, msg, pk_bytes) {
            return Ok(sig);
        }
    }

    Err(Error::InsufficientResponses)
}

/// Decompose a `PolyVecK` into (w₀, w₁) using ML-DSA-44 decompose.
fn decompose_polyveck(w: &PolyVecK) -> (PolyVecK, PolyVecK) {
    use dilithium::{rounding, ML_DSA_44};
    let mode = ML_DSA_44;
    let mut w0 = PolyVecK::zero();
    let mut w1 = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            let a = w.polys[i].coeffs[j].rem_euclid(Q);
            let (r1, r0) = rounding::decompose(mode, a);
            let r0_plus_q = if r0 < 0 { r0 + Q } else { r0 };
            w0.polys[i].coeffs[j] = r0_plus_q;
            w1.polys[i].coeffs[j] = r1;
        }
    }
    (w0, w1)
}

/// Compute c̃ = H(μ ‖ `w₁_packed`).
fn compute_c_tilde(mu: &[u8; 64], w1: &PolyVecK) -> [u8; CTILDEBYTES] {
    use dilithium::{poly::Poly as DPoly, ML_DSA_44};
    let mode = ML_DSA_44;
    let mut h = Shake256::default();
    h.update(mu);
    for i in 0..K {
        let mut packed = [0u8; POLYW1_PACKEDBYTES];
        let mut dp = DPoly::default();
        for j in 0..N {
            dp.coeffs[j] = w1.polys[i].coeffs[j];
        }
        DPoly::polyw1_pack(mode, &mut packed, &dp);
        h.update(&packed);
    }
    let mut c = [0u8; CTILDEBYTES];
    h.finalize_xof().read(&mut c);
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aggregate_commitments_zero() {
        // Two parties, K=2 slots, all zero → result should be all zero
        let w_zero = PolyVecK::zero();
        let party0 = vec![w_zero.clone(), w_zero.clone()];
        let party1 = vec![w_zero.clone(), w_zero.clone()];
        let all = vec![party0, party1];
        let wfinals = aggregate_commitments(&all, 2);
        assert_eq!(wfinals.len(), 2);
        for wf in &wfinals {
            for i in 0..K {
                for j in 0..N {
                    assert_eq!(wf.polys[i].coeffs[j], 0);
                }
            }
        }
    }

    #[test]
    fn test_aggregate_commitments_additive() {
        let mut w_a = PolyVecK::zero();
        let mut w_b = PolyVecK::zero();
        w_a.polys[0].coeffs[0] = 100;
        w_b.polys[0].coeffs[0] = 200;
        let party0 = vec![w_a];
        let party1 = vec![w_b];
        let wfinals = aggregate_commitments(&[party0, party1], 1);
        assert_eq!(wfinals[0].polys[0].coeffs[0], 300);
    }

    #[test]
    fn test_aggregate_responses_zero() {
        let z_zero = PolyVecL::zero();
        let party0 = vec![z_zero.clone()];
        let party1 = vec![z_zero.clone()];
        let all = vec![party0, party1];
        let zfinals = aggregate_responses(&all, 1);
        assert_eq!(zfinals.len(), 1);
        for i in 0..L {
            for j in 0..N {
                assert_eq!(zfinals[0].polys[i].coeffs[j], 0);
            }
        }
    }

    #[test]
    fn test_aggregate_responses_additive() {
        let mut z_a = PolyVecL::zero();
        let mut z_b = PolyVecL::zero();
        z_a.polys[0].coeffs[0] = 1000;
        z_a.polys[1].coeffs[5] = 2000;
        z_b.polys[0].coeffs[0] = 3000;
        z_b.polys[1].coeffs[5] = 4000;
        let all = vec![vec![z_a], vec![z_b]];
        let zfinals = aggregate_responses(&all, 1);
        assert_eq!(zfinals[0].polys[0].coeffs[0], 4000);
        assert_eq!(zfinals[0].polys[1].coeffs[5], 6000);
    }
}
