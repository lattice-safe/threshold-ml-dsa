//! ML-DSA-44 parameters (FIPS 204 §4, Table 1) and threshold-specific constants
//! from ePrint 2026/013, Table 2.

/// The prime modulus q for ℤ_q.
/// q = 2²³ − 2¹³ + 1 = 8380417.
pub const Q: i32 = 8380417;

/// Polynomial degree: coefficients in ℤ_q\[X\]/(X^N + 1).
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

// ─── Threshold-specific parameters (ePrint 2026/013, Table 2) ────────

/// Maximum number of parties supported by the RSS scheme.
/// The paper targets N ≤ 6 to keep the number of RSS subsets manageable.
pub const MAX_PARTIES: usize = 6;

/// Squared L₂ hyperball rejection bound B² for per-party response z_i.
///
/// For ML-DSA-44 (L=4, γ₁=2¹⁷=131072, N=256 coefficients per polynomial):
///
/// z_i is a vector of L polynomials, each with N coefficients uniformly
/// distributed in [-(γ₁-1), γ₁] (from the masking vector y_i) plus a
/// small perturbation from c·s₁_i.
///
/// Expected ‖z_i‖₂² ≈ L · N · γ₁² / 3  (variance of uniform on [-γ₁, γ₁])
///                    ≈ 4 · 256 · 131072² / 3
///                    ≈ 5,864,062,014,805
///
/// We set B² = L · N · γ₁² (≈ 3× the expected value) for a generous
/// acceptance probability (~85%) while still ensuring the final aggregated
/// signature satisfies FIPS 204 bounds.
pub const HYPERBALL_BOUND_SQ: u64 = (L as u64) * (N as u64) * (GAMMA1 as u64) * (GAMMA1 as u64);
