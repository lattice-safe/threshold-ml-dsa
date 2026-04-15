//! Float vector type and hyperball sampling.
//!
//! Implements `FVec` and `SampleHyperball` from the reference Go implementation.
//! The float vector represents an element of R_q^{ℓ+k} in centered floating-point
//! form, used for hyperball-based rejection sampling (ePrint 2026/013, §2.7).
//!
//! Security note: this module uses floating-point arithmetic (`f64`, `libm`) and
//! therefore does not provide constant-time guarantees on typical hardware.
//! Deployments with strict local side-channel requirements should isolate party
//! execution and treat timing leakage in this path as part of the threat model.

extern crate alloc;

use crate::params::{FVEC_DIM, K, L, N, Q};
use crate::poly::{PolyVecK, PolyVecL};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

/// A floating-point vector of dimension (K+L)×N.
///
/// First L×N entries correspond to z^(1) (the transmitted part),
/// last K×N entries correspond to z^(2) (the recoverable part).
#[derive(Clone)]
pub struct FVec {
    pub coeffs: [f64; FVEC_DIM],
}

impl FVec {
    /// Create a zero FVec
    pub fn zero() -> Self {
        Self {
            coeffs: [0.0; FVEC_DIM],
        }
    }

    /// Set self = a + b (pointwise)
    pub fn add_assign(&mut self, other: &FVec) {
        for i in 0..FVEC_DIM {
            self.coeffs[i] += other.coeffs[i];
        }
    }

    /// Convert from polynomial vectors (s₁ ∈ R_q^ℓ, s₂ ∈ R_q^k) to FVec.
    /// Coefficients are centered mod Q into [-Q/2, Q/2].
    pub fn from_polyvecs(s1: &PolyVecL, s2: &PolyVecK) -> Self {
        let mut fv = Self::zero();
        for i in 0..(L + K) {
            for j in 0..N {
                let u = if i < L {
                    s1.polys[i].coeffs[j]
                } else {
                    s2.polys[i - L].coeffs[j]
                };
                // Center mod Q: map [0, Q) → [-(Q-1)/2, (Q-1)/2]
                let mut centered = u + Q / 2;
                let t = centered - Q;
                centered = t + ((t >> 31) & Q);
                centered -= Q / 2;
                fv.coeffs[i * N + j] = centered as f64;
            }
        }
        fv
    }

    /// Round FVec back to polynomial vectors (s₁, s₂).
    /// Adds +Q to negative coefficients.
    pub fn round_to_polyvecs(&self, s1: &mut PolyVecL, s2: &mut PolyVecK) {
        for i in 0..(L + K) {
            for j in 0..N {
                let u = libm::round(self.coeffs[i * N + j]) as i32;
                // Map negative to [0, Q)
                let t = u >> 31;
                let mapped = u + (t & Q);
                if i < L {
                    s1.polys[i].coeffs[j] = mapped;
                } else {
                    s2.polys[i - L].coeffs[j] = mapped;
                }
            }
        }
    }

    /// Check if the ν-scaled L₂-norm exceeds radius `r`.
    ///
    /// Computes: Σ_i (x_i / ν)² for i < L·N, plus Σ_i x_i² for i ≥ L·N.
    /// Returns true if the sum exceeds r².
    pub fn excess(&self, r: f64, nu: f64) -> bool {
        let mut sq = 0.0f64;
        for i in 0..(L + K) {
            for j in 0..N {
                let val = self.coeffs[i * N + j];
                if i < L {
                    sq += val * val / (nu * nu);
                } else {
                    sq += val * val;
                }
            }
        }
        sq > r * r
    }
}

/// Sample a point uniformly on the hyperball B_{R, ℓ+k}(radius) using
/// Box-Muller transform for Gaussian direction + radial scaling.
///
/// The first L·N coordinates are scaled by ν (expansion factor).
///
/// This matches the reference Go `SampleHyperball()` in sample.go.
///
/// # Arguments
/// * `out` — output float vector
/// * `radius` — the hyperball radius (r₁ for randomness, r for target)
/// * `nu` — expansion factor (ν=3 for ML-DSA-44)
/// * `rhop` — 64-byte randomness seed
/// * `nonce` — 16-bit nonce
pub fn sample_hyperball(out: &mut FVec, radius: f64, nu: f64, rhop: &[u8; 64], nonce: u16) {
    // Dimension: (K+L)*N = 8*256 = 2048, plus 2 extra for Box-Muller pairing
    let total_samples = FVEC_DIM + 2;

    // Generate random bytes via SHAKE-256
    let mut hasher = Shake256::default();
    hasher.update(b"H"); // Domain separator (matches Go reference)
    hasher.update(rhop);
    hasher.update(&nonce.to_le_bytes());

    let mut buf = alloc::vec![0u8; total_samples * 8];
    let mut reader = hasher.finalize_xof();
    reader.read(&mut buf);

    // Box-Muller: generate normally-distributed samples
    let mut samples = alloc::vec![0.0f64; total_samples];
    let mut sq = 0.0f64;

    let mut idx = 0;
    while idx < total_samples {
        // Read two uint64 values
        let u1 = u64::from_le_bytes([
            buf[idx * 8],
            buf[idx * 8 + 1],
            buf[idx * 8 + 2],
            buf[idx * 8 + 3],
            buf[idx * 8 + 4],
            buf[idx * 8 + 5],
            buf[idx * 8 + 6],
            buf[idx * 8 + 7],
        ]);
        let u2 = u64::from_le_bytes([
            buf[(idx + 1) * 8],
            buf[(idx + 1) * 8 + 1],
            buf[(idx + 1) * 8 + 2],
            buf[(idx + 1) * 8 + 3],
            buf[(idx + 1) * 8 + 4],
            buf[(idx + 1) * 8 + 5],
            buf[(idx + 1) * 8 + 6],
            buf[(idx + 1) * 8 + 7],
        ]);

        // Convert to float64 in [0, 1)
        let f1 = (u1 as f64) / ((1u128 << 64) as f64);
        let f2 = (u2 as f64) / ((1u128 << 64) as f64);

        // Box-Muller transform
        let log_val = if f1 < f64::MIN_POSITIVE {
            f64::MIN_POSITIVE
        } else {
            f1
        };
        let r = libm::sqrt(-2.0 * libm::log(log_val));
        let theta = 2.0 * core::f64::consts::PI * f2;
        let z1 = r * libm::cos(theta);
        let z2 = r * libm::sin(theta);

        // Apply ν scaling to the first L·N coordinates, then accumulate norm
        // in that scaled space so normalization matches the target geometry.
        let mut s1 = z1;
        if idx < N * L {
            s1 *= nu;
        }
        samples[idx] = s1;
        sq += s1 * s1;

        let mut s2 = z2;
        if idx + 1 < N * L {
            s2 *= nu;
        }
        samples[idx + 1] = s2;
        sq += s2 * s2;

        idx += 2;
    }

    // Normalize to the sphere surface and scale to radius
    let factor = radius / libm::sqrt(sq);
    for (i, sample) in samples.iter().enumerate().take(FVEC_DIM) {
        out.coeffs[i] = sample * factor;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sample_hyperball_radius() {
        let rhop = [0u8; 64];
        let mut fv = FVec::zero();
        sample_hyperball(&mut fv, 252778.0, 3.0, &rhop, 0);

        // The ν-scaled norm should be ≤ radius
        // ‖(x₁/ν, x₂)‖₂ ≤ r
        let mut sq = 0.0f64;
        for i in 0..(L + K) {
            for j in 0..N {
                let val = fv.coeffs[i * N + j];
                if i < L {
                    sq += val * val / (3.0 * 3.0);
                } else {
                    sq += val * val;
                }
            }
        }
        let norm = libm::sqrt(sq);
        assert!(
            norm <= 252778.0 * 1.001, // small tolerance for float precision
            "Hyperball sample norm {} exceeds radius 252778",
            norm
        );
    }

    #[test]
    fn test_fvec_excess() {
        let mut fv = FVec::zero();
        // All zeros should not exceed any positive radius
        assert!(!fv.excess(1.0, 3.0));

        // Set a large value
        fv.coeffs[0] = 1e6;
        assert!(fv.excess(1000.0, 3.0));
    }

    #[test]
    fn test_from_and_round_roundtrip() {
        let mut s1 = PolyVecL::zero();
        let mut s2 = PolyVecK::zero();
        // Set some small values
        s1.polys[0].coeffs[0] = 1;
        s1.polys[0].coeffs[1] = Q - 1; // -1 mod Q
        s2.polys[0].coeffs[0] = 5;

        let fv = FVec::from_polyvecs(&s1, &s2);
        assert!((fv.coeffs[0] - 1.0).abs() < 0.01);
        assert!((fv.coeffs[1] - (-1.0)).abs() < 0.01);
        assert!((fv.coeffs[L * N] - 5.0).abs() < 0.01);

        // Round-trip back
        let mut s1_out = PolyVecL::zero();
        let mut s2_out = PolyVecK::zero();
        fv.round_to_polyvecs(&mut s1_out, &mut s2_out);
        assert_eq!(s1_out.polys[0].coeffs[0], 1);
        assert_eq!(s1_out.polys[0].coeffs[1], Q - 1);
        assert_eq!(s2_out.polys[0].coeffs[0], 5);
    }
}
