//! Coordinator / Aggregator for threshold ML-DSA signing.
//!
//! Implements the `Combine` function from the Go reference.
//! Tries each of the K parallel commitment/response slots until
//! one produces a valid FIPS 204 signature.
//!
//! The coordinator:
//! 1. Aggregates commitments w_final = Σ w_i
//! 2. Aggregates responses z_final = Σ z_i
//! 3. For each k ∈ [0, K): checks if z_k produces a valid signature

extern crate alloc;
use alloc::vec::Vec;

use crate::error::Error;
use crate::params::*;
use crate::poly::*;
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
    /// The aggregated response vector z ∈ ℤ_q^L.
    pub z: PolyVecL,
    /// The hint vector h (packed as OMEGA + K bytes).
    pub h: Vec<u8>,
}

impl Signature {
    /// Encode the signature into the standard FIPS 204 byte format.
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
/// For each k ∈ [0, K): wfinal[k] = Σ_i w_{i,k}
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
/// For each k ∈ [0, K): zfinal[k] = Σ_i z_{i,k}
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
/// 1. Check ‖z_k‖∞ < γ₁ - β
/// 2. Compute Az_k - 2^d·c·t₁ = "approximated w"  
/// 3. Check ‖δ‖∞ = ‖approx_w - w_final‖∞ < γ₂
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
        poly::Poly as DPoly,
        polyvec::{
            PolyVecK as DPolyVecK, PolyVecL as DPolyVecL,
            matrix_expand, matrix_pointwise_montgomery,
        },
        packing::{unpack_pk, pack_sig},
        ML_DSA_44,
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

        // Compute c̃ = H(μ ‖ w1_packed)
        let (_, w1) = decompose_polyveck(&wfinals[k]);
        let c_tilde = compute_c_tilde(&mu, &w1);

        // Convert z to dilithium-rs format
        let mut zh_d = DPolyVecL::default();
        for i in 0..L {
            for j in 0..N {
                zh_d.vec[i].coeffs[j] = zfinals[k].polys[i].coeffs[j];
            }
            zh_d.vec[i].ntt();
        }

        // Compute Az in NTT domain
        let mut az = DPolyVecK::default();
        matrix_pointwise_montgomery(mode, &mut az, &mat, &zh_d);
        for i in 0..K {
            az.vec[i].reduce();
            az.vec[i].invntt_tomont();
        }

        // Expand challenge c
        let mut ch = DPoly::default();
        DPoly::challenge(mode, &mut ch, &c_tilde);
        ch.ntt();

        // Compute Az - 2^d · c · t₁
        let mut az2dct1 = DPolyVecK::default();
        for i in 0..K {
            // t1 * 2^d
            az2dct1.vec[i] = t1.vec[i].clone();
            az2dct1.vec[i].shiftl();
            az2dct1.vec[i].ntt();
            let tmp = az2dct1.vec[i].clone();
            DPoly::pointwise_montgomery(&mut az2dct1.vec[i], &tmp, &ch);
        }
        // approx_w = Az - 2^d·c·t₁
        for i in 0..K {
            let tmp = az2dct1.vec[i].clone();
            DPoly::sub(&mut az2dct1.vec[i], &az.vec[i], &tmp);
            az2dct1.vec[i].reduce();
            az2dct1.vec[i].invntt_tomont();
            az2dct1.vec[i].caddq();
        }

        // δ = approx_w - w_final
        let mut delta_exceeds = false;
        let mut f = DPolyVecK::default();
        for i in 0..K {
            for j in 0..N {
                let wf = if wfinals[k].polys[i].coeffs[j] < 0 {
                    wfinals[k].polys[i].coeffs[j] + Q
                } else {
                    wfinals[k].polys[i].coeffs[j]
                };
                f.vec[i].coeffs[j] = az2dct1.vec[i].coeffs[j] - wf;
                // Center
                if f.vec[i].coeffs[j] > Q / 2 {
                    f.vec[i].coeffs[j] -= Q;
                } else if f.vec[i].coeffs[j] < -(Q / 2) {
                    f.vec[i].coeffs[j] += Q;
                }
                if f.vec[i].coeffs[j].unsigned_abs() as i32 >= GAMMA2 {
                    delta_exceeds = true;
                }
            }
        }

        if delta_exceeds {
            continue;
        }

        // Compute w₀ + f, then make hint
        let (w0, w1_dec) = decompose_polyveck(&wfinals[k]);
        let mut w0pf = DPolyVecK::default();
        for i in 0..K {
            for j in 0..N {
                let w0_val = if w0.polys[i].coeffs[j] < 0 {
                    w0.polys[i].coeffs[j] + Q
                } else {
                    w0.polys[i].coeffs[j]
                };
                w0pf.vec[i].coeffs[j] = w0_val + f.vec[i].coeffs[j];
                // Normalize
                if w0pf.vec[i].coeffs[j] < 0 {
                    w0pf.vec[i].coeffs[j] += Q;
                } else if w0pf.vec[i].coeffs[j] >= Q {
                    w0pf.vec[i].coeffs[j] -= Q;
                }
            }
        }

        // Make hint
        let mut hint = DPolyVecK::default();
        let mut w1_d = DPolyVecK::default();
        for i in 0..K {
            for j in 0..N {
                w1_d.vec[i].coeffs[j] = w1_dec.polys[i].coeffs[j];
            }
        }
        let mut hint_pop = 0usize;
        for i in 0..K {
            hint_pop += DPoly::make_hint(mode, &mut hint.vec[i], &w0pf.vec[i], &w1_d.vec[i]);
        }

        if hint_pop > OMEGA {
            continue;
        }

        // Success! Pack the signature
        let mut z_d = DPolyVecL::default();
        for i in 0..L {
            for j in 0..N {
                z_d.vec[i].coeffs[j] = zfinals[k].polys[i].coeffs[j];
            }
        }

        let mut sig = [0u8; SIG_BYTES];
        pack_sig(mode, &mut sig, &c_tilde, &z_d, &hint);

        return Ok(sig);
    }

    Err(Error::InsufficientResponses)
}

/// Decompose a PolyVecK into (w₀, w₁) using ML-DSA-44 decompose.
fn decompose_polyveck(w: &PolyVecK) -> (PolyVecK, PolyVecK) {
    use dilithium::{rounding, ML_DSA_44};
    let mode = ML_DSA_44;
    let mut w0 = PolyVecK::zero();
    let mut w1 = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            let a = w.polys[i].coeffs[j];
            let a_pos = if a < 0 { a + Q } else { a };
            let (r1, r0) = rounding::decompose(mode, a_pos);
            w0.polys[i].coeffs[j] = r0;
            w1.polys[i].coeffs[j] = r1;
        }
    }
    (w0, w1)
}

/// Compute c̃ = H(μ ‖ w₁_packed).
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
    // v0.3.0: End-to-end tests will be added after wiring the full
    // round1 → round2 → round3 → combine pipeline.
}
