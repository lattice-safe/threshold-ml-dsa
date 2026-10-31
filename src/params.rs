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

/// Secret key size: ρ ‖ K ‖ tr ‖ s₁ ‖ s₂ ‖ t₀ (2560 bytes in FIPS 204).
pub const SK_BYTES: usize = 2 * SEEDBYTES
    + TRBYTES
    + L * POLYETA_PACKEDBYTES
    + K * POLYETA_PACKEDBYTES
    + K * POLYT0_PACKEDBYTES;

/// Challenge hash (c̃) size.
pub const CTILDEBYTES: usize = 32;

/// Signature size: c̃ ‖ z ‖ h.
pub const SIG_BYTES: usize = CTILDEBYTES + L * POLYZ_PACKEDBYTES + OMEGA + K;

// ─── FIPS 204 ML-DSA-65 Parameters ───────────────────────────────────────

/// Parameter constants for ML-DSA-65 (NIST Security Category 3).
pub mod ml_dsa_65 {
    pub use super::{N, Q, SEEDBYTES, TRBYTES};

    pub const K: usize = 6;
    pub const L: usize = 5;
    pub const ETA: u32 = 4;
    pub const TAU: u32 = 49;
    pub const GAMMA1: i32 = 1 << 19;
    pub const GAMMA2: i32 = (Q - 1) / 32;
    pub const BETA: i32 = (TAU * ETA) as i32;
    pub const OMEGA: usize = 55;

    pub const POLYETA_PACKEDBYTES: usize = 128;
    pub const POLYZ_PACKEDBYTES: usize = 640;
    pub const POLYW1_PACKEDBYTES: usize = 128;
    pub const CTILDEBYTES: usize = 48;

    pub const PK_BYTES: usize = 1952;
    pub const SK_BYTES: usize = 4032;
    pub const SIG_BYTES: usize = 3309;
}

// ─── FIPS 204 ML-DSA-87 Parameters ───────────────────────────────────────

/// Parameter constants for ML-DSA-87 (NIST Security Category 5).
pub mod ml_dsa_87 {
    pub use super::{N, Q, SEEDBYTES, TRBYTES};

    pub const K: usize = 8;
    pub const L: usize = 7;
    pub const ETA: u32 = 2;
    pub const TAU: u32 = 60;
    pub const GAMMA1: i32 = 1 << 19;
    pub const GAMMA2: i32 = (Q - 1) / 32;
    pub const BETA: i32 = (TAU * ETA) as i32;
    pub const OMEGA: usize = 75;

    pub const POLYETA_PACKEDBYTES: usize = 96;
    pub const POLYZ_PACKEDBYTES: usize = 640;
    pub const POLYW1_PACKEDBYTES: usize = 128;
    pub const CTILDEBYTES: usize = 64;

    pub const PK_BYTES: usize = 2592;
    pub const SK_BYTES: usize = 4896;
    pub const SIG_BYTES: usize = 4627;
}

// ─── Threshold-specific parameters (ePrint 2026/013) ─────────────────

/// Maximum number of parties supported by the RSS scheme (N ≤ 8).
///
/// ePrint 2026/013 gives parameters for N ≤ 6; N = 7, 8 are derived here
/// with the reference scripts (see [`get_threshold_params`]). Not every
/// (T, 8) is supported.
pub const MAX_PARTIES: usize = 8;

/// Dimension of the `FVec` float vector: (K+L)×N coefficients.
pub const FVEC_DIM: usize = (K + L) * N;

/// Threshold-specific parameters for a given (T, N) pair.
/// From ePrint 2026/013 full version (ML-DSA-44; values as in the Go reference).
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

impl ThresholdParams {
    /// Check that these parameters are exactly the canonical set for (T, N).
    ///
    /// The fields are public, so callers can build arbitrary values; a radius
    /// of NaN or +∞ would disable hyperball rejection and release every
    /// response. Protocol entry points call this before using `self`.
    pub fn validate(&self) -> Result<(), crate::error::Error> {
        let canon = get_threshold_params(self.t, self.n).ok_or(crate::error::Error::InvalidParameters)?;
        if self.k_reps != canon.k_reps
            || self.r.to_bits() != canon.r.to_bits()
            || self.r1.to_bits() != canon.r1.to_bits()
            || self.nu.to_bits() != canon.nu.to_bits()
        {
            return Err(crate::error::Error::InvalidParameters);
        }
        Ok(())
    }
}

/// Lookup the paper parameters for a given (T, N) pair.
///
/// Returns `None` if (T, N) is unsupported: outside 2 ≤ T ≤ N ≤ 8, or one
/// of (5, 8), (6, 8), (7, 8).
///
/// N ≤ 6: ePrint 2026/013 full version, Fig. 9 (values as in the Go
/// reference). Note that Fig. 9's (2,4), (4,5), (5,5) radii correspond to
/// ν = 3.5, 2.5, 2.5 in the reference script, and its K for (4,6), (5,6),
/// (6,6) is below what the script gives (per-attempt success ≈ 0.39–0.46
/// instead of 1/2); these values are kept for compatibility.
///
/// N = 7, 8: derived with the reference script `params/hyperball.sage`
/// (ν = 3, per-attempt success 1/2, Monte Carlo for the final-check
/// probabilities) and the max-flow partitions of `params/recover.py`. K is
/// the value at the −2σ lower bound of the estimated success probability.
/// Per party and attempt, communication is about 32 + 5248·K bytes and
/// signer state about 19.3 KB·K: e.g. (4, 8) ≈ 6.3 MB / 23 MB. (5, 8),
/// (6, 8) and (7, 8) would need K = 14 785, 38 473 and 5 162 (76–198 MB of
/// traffic, 98–728 MB of memory per attempt) and are not supported.
///
/// All sets use the same multiplicative factor ν = 3.
#[must_use] 
pub fn get_threshold_params(t: u8, n: u8) -> Option<ThresholdParams> {
    // The table below lists exactly the supported (T, N).
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
        (2, 7) => (313939.0, 314051.0, 4),
        (3, 7) => (285453.0, 285601.0, 31),
        (4, 7) => (261941.0, 262110.0, 322),
        (5, 7) => (241018.0, 241166.0, 1071),
        (6, 7) => (220521.0, 220633.0, 1008),
        (7, 7) => (192000.0, 192055.0, 91),
        (2, 8) => (313939.0, 314051.0, 4),
        (3, 8) => (291228.0, 291406.0, 58),
        (4, 8) => (267837.0, 268048.0, 1195),
        (8, 8) => (182918.0, 182974.0, 227),
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

/// Number of RSS subsets for a given (N, T): C(N, N-T+1); 0 unless 1 ≤ T ≤ N.
#[must_use] 
pub const fn num_subsets(n: u8, t: u8) -> usize {
    if t == 0 || t > n {
        return 0;
    }
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
        // All 15 valid (T, N) pairs from the parameter table
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
            let p = get_threshold_params(t, n);
            assert!(p.is_some(), "Missing params for ({t}, {n})");
            let p = p.unwrap();
            assert!(p.r1 >= p.r, "r1 < r for ({}, {})", t, n);
            assert!(p.k_reps >= 2, "K < 2 for ({}, {})", t, n);
            assert!((p.nu - 3.0).abs() < 1e-10, "ν ≠ 3 for ({}, {})", t, n);
        }
    }

    #[test]
    fn test_invalid_params_rejected() {
        assert!(get_threshold_params(1, 2).is_none()); // T < 2
        assert!(get_threshold_params(3, 2).is_none()); // T > N
        assert!(get_threshold_params(2, 9).is_none()); // N > 8
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

#[cfg(test)]
mod coverage_tests {
    use super::*;

    #[test]
    fn unsupported_configurations_and_binomial_edges() {
        assert!(get_threshold_params(1, 2).is_none());
        assert!(get_threshold_params(3, 2).is_none());
        assert!(get_threshold_params(2, 9).is_none());
        assert!(get_threshold_params(5, 8).is_none()); // impractical K
        assert_eq!(binomial(2, 3), 0);
        assert_eq!(binomial(5, 0), 1);
        assert_eq!(binomial(5, 5), 1);
        assert_eq!(binomial(6, 4), 15);
        let mut p = get_threshold_params(2, 3).unwrap();
        assert!(p.validate().is_ok());
        p.nu = f64::NAN;
        assert!(p.validate().is_err());
        p = get_threshold_params(2, 3).unwrap();
        p.t = 7;
        assert!(p.validate().is_err());
    }
}

#[cfg(test)]
mod review_tests {
    use super::*;

    #[test]
    fn n7_n8_rows_exist_except_impractical_ones() {
        for n in 7..=8u8 {
            for t in 2..=n {
                let unsupported = n == 8 && (5..=7).contains(&t);
                let p = get_threshold_params(t, n);
                assert_eq!(p.is_none(), unsupported, "({t}, {n})");
                if let Some(p) = p {
                    assert!(p.validate().is_ok() && p.r < p.r1 && p.k_reps >= 4);
                }
            }
        }
    }

    #[test]
    fn num_subsets_handles_invalid_input() {
        assert_eq!(num_subsets(2, 5), 0);
        assert_eq!(num_subsets(3, 0), 0);
        assert_eq!(num_subsets(6, 3), 15); // C(6, 4)
        assert_eq!(num_subsets(6, 4), 20); // C(6, 3)
    }
}
