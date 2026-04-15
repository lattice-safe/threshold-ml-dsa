//! ML-DSA-44 parameters (FIPS 204 §4, Table 1) and threshold-specific constants
//! from ePrint 2026/013, Figures 8-9.

// ─── FIPS 204 ML-DSA-44 Parameters ───────────────────────────────────────

/// The prime modulus q for `ℤ_q`.
/// q = 2²³ − 2¹³ + 1 = 8380417.
pub const Q: i32 = 8380417;

/// Polynomial degree: coefficients in `ℤ_q`\[X\]/(X^N + 1).
pub const N: usize = 256;

/// Dropped bits from public key (power-of-2 rounding parameter).
pub const D: u32 = 13;

/// Number of rows in the matrix A (k × l).
pub const K: usize = 4;

/// Number of columns in the matrix A (k × l).
pub const L: usize = 4;

/// CBD sampling bound η for secret key coefficients s₁, s₂.
pub const ETA: u32 = 2;

/// Number of ±1 coefficients in the challenge polynomial c.
pub const TAU: u32 = 39;

/// Masking vector bound: γ₁ = 2¹⁷.
pub const GAMMA1: i32 = 1 << 17;

/// Decomposition parameter: γ₂ = (q − 1) / 88.
pub const GAMMA2: i32 = (Q - 1) / 88;

/// Norm bound for challenge × secret: β = τ · η.
pub const BETA: i32 = (TAU * ETA) as i32;

/// Maximum number of 1-bits in hint polynomial vector.
pub const OMEGA: usize = 80;

// ─── Sizes (bytes) ────────────────────────────────────────────────────
// FIPS 204 §5, Table 2

/// Seed size (ρ, ρ', K).
pub const SEEDBYTES: usize = 32;

/// CRH output size (μ, tr).
pub const CRHBYTES: usize = 64;

/// tr size in the secret key.
pub const TRBYTES: usize = 64;

/// Packed polynomial t₁ size.
pub const POLYT1_PACKEDBYTES: usize = 320;

/// Packed polynomial t₀ size.
pub const POLYT0_PACKEDBYTES: usize = 416;

/// Packed polynomial with η-bounded coefficients.
pub const POLYETA_PACKEDBYTES: usize = 96; // η=2 → 3 bits/coeff → 256×3/8

/// Packed polynomial z (γ₁-bounded).
pub const POLYZ_PACKEDBYTES: usize = 576; // γ₁=2¹⁷ → 18 bits/coeff

/// Packed w₁ polynomial.
pub const POLYW1_PACKEDBYTES: usize = 192;

/// Public key size: ρ ‖ t₁.
pub const PK_BYTES: usize = SEEDBYTES + K * POLYT1_PACKEDBYTES;

/// Secret key size: ρ ‖ K ‖ tr ‖ s₁ ‖ s₂ ‖ t₀.
pub const SK_BYTES: usize = 3 * SEEDBYTES
    + TRBYTES
    + L * POLYETA_PACKEDBYTES
    + K * POLYETA_PACKEDBYTES
    + K * POLYT0_PACKEDBYTES;

/// Challenge hash (c̃) size.
pub const CTILDEBYTES: usize = 32;

/// Signature size: c̃ ‖ z ‖ h.
pub const SIG_BYTES: usize = CTILDEBYTES + L * POLYZ_PACKEDBYTES + OMEGA + K;

// ─── Threshold-specific parameters (ePrint 2026/013) ─────────────────

/// Maximum number of parties supported by the RSS scheme.
/// The paper targets N ≤ 6 to keep the number of RSS subsets manageable.
pub const MAX_PARTIES: usize = 6;

/// Dimension of the `FVec` float vector: (K+L)×N coefficients.
pub const FVEC_DIM: usize = (K + L) * N;

/// Threshold-specific parameters for a given (T, N) pair.
/// From ePrint 2026/013, Figure 8 (ML-DSA-44).
///
/// - `r`: target hyperball radius for `χ_z`
/// - `r1`: randomness hyperball radius for `χ_r`  
/// - `k_reps`: number of parallel repetitions K
/// - `nu`: expansion factor for first ℓ coordinates (ν=3 for ML-DSA-44)
#[derive(Debug, Clone, Copy)]
pub struct ThresholdParams {
    /// Threshold T — minimum signers needed
    pub t: u8,
    /// Total parties N
    pub n: u8,
    /// Target ball radius r (for `χ_z`)
    pub r: f64,
    /// Randomness ball radius r₁ (for `χ_r`), r₁ ≥ r
    pub r1: f64,
    /// Number of parallel protocol repetitions K
    pub k_reps: u16,
    /// Expansion factor ν for the first ℓ·n coords of the response
    pub nu: f64,
}

/// Lookup the paper-exact parameters for a given (T, N) pair.
///
/// Returns `None` if (T, N) is not in the supported range
/// (2 ≤ T ≤ N ≤ 6).
///
/// Parameters from ePrint 2026/013, Figure 8 (ML-DSA-44).
/// All sets use the same multiplicative factor ν = 3.
#[must_use] 
pub fn get_threshold_params(t: u8, n: u8) -> Option<ThresholdParams> {
    if t < 2 || t > n || n > 6 {
        return None;
    }

    let nu = 3.0_f64;

    let (r, r1, k_reps) = match (t, n) {
        (2, 2) => (252778.0, 252833.0, 2u16),
        (2, 3) => (310060.0, 310138.0, 3),
        (3, 3) => (246490.0, 246546.0, 4),
        (2, 4) => (305919.0, 305997.0, 3),
        (3, 4) => (279235.0, 279314.0, 7),
        (4, 4) => (243463.0, 243519.0, 8),
        (2, 5) => (285363.0, 285459.0, 3),
        (3, 5) => (282800.0, 282912.0, 14),
        (4, 5) => (259427.0, 259526.0, 30),
        (5, 5) => (239924.0, 239981.0, 16),
        (2, 6) => (300265.0, 300362.0, 4),
        (3, 6) => (277014.0, 277139.0, 19),
        (4, 6) => (268705.0, 268831.0, 74),
        (5, 6) => (250590.0, 250686.0, 100),
        (6, 6) => (219245.0, 219301.0, 37),
        _ => return None,
    };

    Some(ThresholdParams {
        t,
        n,
        r,
        r1,
        k_reps,
        nu,
    })
}

/// Compute C(n, k) = n! / (k! (n-k)!)
#[must_use] 
pub const fn binomial(n: u8, k: u8) -> usize {
    if k > n {
        return 0;
    }
    if k == 0 || k == n {
        return 1;
    }
    // Use the smaller of k and n-k for efficiency
    let k = if k > n - k { n - k } else { k };
    let mut result: usize = 1;
    let mut i: u8 = 0;
    while i < k {
        result = result * (n - i) as usize / (i + 1) as usize;
        i += 1;
    }
    result
}

/// Number of RSS subsets for a given (N, T): C(N, N-T+1)
#[must_use] 
pub const fn num_subsets(n: u8, t: u8) -> usize {
    binomial(n, n - t + 1)
}

// ─── Legacy constant (for backward compatibility in tests) ───────────

/// Squared L₂ hyperball rejection bound B² (legacy, superseded by ThresholdParams.r²).
pub const HYPERBALL_BOUND_SQ: u64 = (L as u64) * (N as u64) * (GAMMA1 as u64) * (GAMMA1 as u64);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_param_sets_exist() {
        // All 15 valid (T, N) pairs from Figure 8
        let valid = [
            (2, 2),
            (2, 3),
            (3, 3),
            (2, 4),
            (3, 4),
            (4, 4),
            (2, 5),
            (3, 5),
            (4, 5),
            (5, 5),
            (2, 6),
            (3, 6),
            (4, 6),
            (5, 6),
            (6, 6),
        ];
        for (t, n) in valid {
            let p = get_threshold_params(t, n)
                .unwrap_or_else(|| panic!("Missing params for ({}, {})", t, n));
            assert!(p.r1 >= p.r, "r1 < r for ({}, {})", t, n);
            assert!(p.k_reps >= 2, "K < 2 for ({}, {})", t, n);
            assert!((p.nu - 3.0).abs() < 1e-10, "ν ≠ 3 for ({}, {})", t, n);
        }
    }

    #[test]
    fn test_invalid_params_rejected() {
        assert!(get_threshold_params(1, 2).is_none()); // T < 2
        assert!(get_threshold_params(3, 2).is_none()); // T > N
        assert!(get_threshold_params(2, 7).is_none()); // N > 6
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(3, 2), 3);
        assert_eq!(binomial(4, 2), 6);
        assert_eq!(binomial(6, 3), 20);
        assert_eq!(binomial(6, 5), 6);
    }

    #[test]
    fn test_num_subsets() {
        assert_eq!(num_subsets(3, 2), 3); // C(3,2) = 3
        assert_eq!(num_subsets(4, 3), 6); // C(4,2) = 6
        assert_eq!(num_subsets(6, 2), 6); // C(6,5) = 6
    }
}
