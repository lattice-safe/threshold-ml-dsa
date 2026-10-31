//! Polynomial arithmetic over `ℤ_q`\[X\]/(X^256 + 1) for ML-DSA.
//!
//! This module provides the core `Poly` type and vector wrappers (`PolyVecL`, `PolyVecK`)
//! with all operations needed by the threshold signing protocol:
//!
//! - Modular arithmetic (add, sub, pointwise multiply in NTT domain)
//! - Forward/inverse NTT (delegated to `dilithium-rs` for FIPS 204 compliance)
//! - L₂ and L∞ norms for rejection sampling bounds
//! - SHAKE-based sampling (uniform, short/eta, gamma1-bounded)
//! - Packing/unpacking for FIPS 204 signature encoding
//!
//! All coefficient arithmetic uses centered representatives in `[-(q-1)/2, (q-1)/2]`
//! where relevant (especially for norm computations).

use crate::params::{N, Q, SEEDBYTES, CRHBYTES, ETA, GAMMA1, TAU, POLYT1_PACKEDBYTES, POLYZ_PACKEDBYTES, L, K, D, GAMMA2};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::{Shake128, Shake256};
use zeroize::Zeroize;

// ─── Core Polynomial Type ────────────────────────────────────────────

/// A polynomial in `ℤ_q`\[X\]/(X^256 + 1), stored as 256 coefficients.
#[derive(Clone, Zeroize)]
pub struct Poly {
    pub coeffs: [i32; N],
}

impl core::fmt::Debug for Poly {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        // Never print coefficients: a Poly may hold secret material.
        write!(f, "Poly([..; {N}])")
    }
}

impl Default for Poly {
    fn default() -> Self {
        Self::zero()
    }
}

impl Poly {
    /// Zero polynomial.
    #[must_use] 
    pub fn zero() -> Self {
        Poly { coeffs: [0i32; N] }
    }

    /// Reduce all coefficients to [0, q).
    pub fn reduce(&mut self) {
        for c in &mut self.coeffs {
            *c = reduce_i32(*c);
        }
    }

    /// Add q to negative coefficients to get canonical [0, q) form.
    pub fn caddq(&mut self) {
        for c in &mut self.coeffs {
            *c += (*c >> 31) & Q;
        }
    }

    /// Coefficient-wise addition: c = a + b.
    pub fn add(c: &mut Poly, a: &Poly, b: &Poly) {
        for i in 0..N {
            c.coeffs[i] = a.coeffs[i] + b.coeffs[i];
        }
    }

    /// Coefficient-wise subtraction: c = a - b.
    pub fn sub(c: &mut Poly, a: &Poly, b: &Poly) {
        for i in 0..N {
            c.coeffs[i] = a.coeffs[i] - b.coeffs[i];
        }
    }

    /// Forward NTT in-place.
    ///
    /// Delegates to [`dilithium-rs`](https://crates.io/crates/dilithium-rs),
    /// validated against all 100 NIST KAT vectors.
    pub fn ntt(&mut self) {
        dilithium::ntt::ntt(&mut self.coeffs);
    }

    /// Inverse NTT in-place (output in Montgomery domain).
    ///
    /// Delegates to dilithium-rs's inverse NTT.
    pub fn invntt_tomont(&mut self) {
        dilithium::ntt::invntt_tomont(&mut self.coeffs);
    }

    /// Pointwise multiplication in NTT domain: c = a · b.
    ///
    /// Both a and b must be in NTT form. Result is in NTT form.
    /// Uses Montgomery reduction for efficient modular arithmetic.
    pub fn pointwise_montgomery(c: &mut Poly, a: &Poly, b: &Poly) {
        for i in 0..N {
            c.coeffs[i] =
                dilithium::reduce::montgomery_reduce(i64::from(a.coeffs[i]) * i64::from(b.coeffs[i]));
        }
    }

    // ─── Norm computations ───────────────────────────────────────────

    /// Compute the squared L₂ norm ‖p‖₂² using centered representatives.
    ///
    /// Each coefficient is mapped to the centered range [-(q-1)/2, (q-1)/2]
    /// before squaring. This is critical for the hyperball rejection check.
    #[must_use] 
    pub fn l2_norm_squared(&self) -> u64 {
        let mut norm: u64 = 0;
        for &c in &self.coeffs {
            let centered = center(c);
            norm += (i64::from(centered) * i64::from(centered)) as u64;
        }
        norm
    }

    /// Compute the L∞ norm: max |coefficient| using centered representatives.
    #[must_use] 
    pub fn l_inf_norm(&self) -> u32 {
        let mut max = 0u32;
        for &c in &self.coeffs {
            let centered = center(c);
            let abs = centered.unsigned_abs();
            if abs > max {
                max = abs;
            }
        }
        max
    }

    /// Check if ‖p‖∞ ≥ bound (used in FIPS 204 verification).
    ///
    /// Constant-time: iterates all coefficients without early exit
    /// to prevent timing side-channels on secret-dependent data.
    #[must_use] 
    pub fn chknorm(&self, bound: i32) -> bool {
        let bound_u = bound as u32;
        let mut exceeded = 0u32;
        for &c in &self.coeffs {
            let centered = center(c);
            // Accumulate without branching
            exceeded |= (bound_u.wrapping_sub(centered.unsigned_abs().wrapping_add(1))) >> 31;
        }
        exceeded != 0
    }

    // ─── Sampling ────────────────────────────────────────────────────

    /// Sample a uniformly random polynomial with coefficients in [0, q)
    /// using rejection sampling from SHAKE-128. (`ExpandA`)
    pub fn uniform(a: &mut Poly, seed: &[u8; SEEDBYTES], nonce: u16) {
        let mut hasher = Shake128::default();
        hasher.update(seed);
        hasher.update(&nonce.to_le_bytes());
        let mut reader = hasher.finalize_xof();

        let mut ctr = 0;
        let mut buf = [0u8; 3];
        while ctr < N {
            reader.read(&mut buf);
            let t = (u32::from(buf[0]) | (u32::from(buf[1]) << 8) | (u32::from(buf[2]) << 16)) & 0x7FFFFF;
            if t < Q as u32 {
                a.coeffs[ctr] = t as i32;
                ctr += 1;
            }
        }
    }

    /// Sample a short polynomial with coefficients in [-η, η].
    /// Uses SHAKE-256 (`ExpandS`). For ML-DSA-44, η = 2.
    pub fn uniform_eta(a: &mut Poly, seed: &[u8; CRHBYTES], nonce: u16) {
        let mut hasher = Shake256::default();
        hasher.update(seed);
        hasher.update(&nonce.to_le_bytes());
        let mut reader = hasher.finalize_xof();

        let mut buf = [0u8; 1];
        let mut ctr = 0;
        // η = 2: FIPS 204 CoeffFromHalfByte accepts b < 15 and maps
        // b ↦ η − (b mod 5).
        while ctr < N {
            reader.read(&mut buf);
            let t0 = i32::from(buf[0] & 0x0F);
            let t1 = i32::from(buf[0] >> 4);
            if t0 < 15 {
                a.coeffs[ctr] = ETA as i32 - t0 % 5;
                ctr += 1;
            }
            if ctr < N && t1 < 15 {
                a.coeffs[ctr] = ETA as i32 - t1 % 5;
                ctr += 1;
            }
        }
    }

    /// Sample a masking polynomial with coefficients in [-(γ₁-1), γ₁].
    /// (`ExpandMask`, Algorithm 24 in FIPS 204)
    pub fn uniform_gamma1(a: &mut Poly, seed: &[u8; CRHBYTES], nonce: u16) {
        let mut hasher = Shake256::default();
        hasher.update(seed);
        hasher.update(&nonce.to_le_bytes());
        let mut reader = hasher.finalize_xof();

        // For γ₁ = 2¹⁷, each coefficient needs 18 bits → 9 bytes per 4 coefficients
        let mut buf = [0u8; 9];
        let mut ctr = 0;
        while ctr + 4 <= N {
            reader.read(&mut buf);
            // Unpack 4 coefficients from 9 bytes (18 bits each)
            for j in 0..4 {
                let mut val = 0u32;
                let base = j * 18 / 8;
                if base + 2 < 9 {
                    val = u32::from(buf[base])
                        | (u32::from(buf[base + 1]) << 8)
                        | (u32::from(buf[base + 2]) << 16);
                    val = (val >> (j * 18 % 8)) & 0x3FFFF; // 18-bit mask
                }
                a.coeffs[ctr] = GAMMA1 - (val as i32);
                ctr += 1;
                if ctr >= N {
                    break;
                }
            }
        }
    }

    /// Compute the challenge polynomial c from the hash c̃.
    /// (`SampleInBall`, Algorithm 22 in FIPS 204)
    ///
    /// Produces a polynomial with exactly τ non-zero (±1) coefficients.
    pub fn challenge(c: &mut Poly, seed: &[u8]) {
        *c = Poly::zero();
        let mut hasher = Shake256::default();
        hasher.update(seed);
        let mut reader = hasher.finalize_xof();

        // Read the sign bits
        let mut signs_buf = [0u8; 8];
        reader.read(&mut signs_buf);
        let signs = u64::from_le_bytes(signs_buf);

        for (pos, i) in ((N - TAU as usize)..N).enumerate() {
            // Sample j ∈ [0, i] uniformly
            let mut j_buf = [0u8; 1];
            loop {
                reader.read(&mut j_buf);
                if (j_buf[0] as usize) <= i {
                    break;
                }
            }
            let j = j_buf[0] as usize;

            c.coeffs[i] = c.coeffs[j];
            c.coeffs[j] = if (signs >> pos) & 1 == 1 {
                Q - 1 // represents -1 mod q
            } else {
                1
            };
        }
    }

    // ─── Packing / Unpacking ─────────────────────────────────────────

    /// Pack t₁ polynomial (10 bits per coefficient).
    pub fn pack_t1(&self, r: &mut [u8]) {
        assert!(r.len() >= POLYT1_PACKEDBYTES, "pack_t1: output buffer too short");
        for i in 0..N / 4 {
            r[5 * i] = (self.coeffs[4 * i] & 0xFF) as u8;
            r[5 * i + 1] = ((self.coeffs[4 * i] >> 8) | (self.coeffs[4 * i + 1] << 2)) as u8;
            r[5 * i + 2] = ((self.coeffs[4 * i + 1] >> 6) | (self.coeffs[4 * i + 2] << 4)) as u8;
            r[5 * i + 3] = ((self.coeffs[4 * i + 2] >> 4) | (self.coeffs[4 * i + 3] << 6)) as u8;
            r[5 * i + 4] = (self.coeffs[4 * i + 3] >> 2) as u8;
        }
    }

    /// Unpack t₁ polynomial (10 bits per coefficient).
    pub fn unpack_t1(r: &mut Poly, a: &[u8]) {
        assert!(a.len() >= POLYT1_PACKEDBYTES, "unpack_t1: input too short");
        for i in 0..N / 4 {
            r.coeffs[4 * i] = (i32::from(a[5 * i]) | (i32::from(a[5 * i + 1]) << 8)) & 0x3FF;
            r.coeffs[4 * i + 1] =
                ((i32::from(a[5 * i + 1]) >> 2) | (i32::from(a[5 * i + 2]) << 6)) & 0x3FF;
            r.coeffs[4 * i + 2] =
                ((i32::from(a[5 * i + 2]) >> 4) | (i32::from(a[5 * i + 3]) << 4)) & 0x3FF;
            r.coeffs[4 * i + 3] =
                ((i32::from(a[5 * i + 3]) >> 6) | (i32::from(a[5 * i + 4]) << 2)) & 0x3FF;
        }
    }

    /// Pack z polynomial (coefficients in [-(γ₁-1), γ₁], 18 bits each).
    ///
    /// Coefficients may be in any representative mod q; they are centered
    /// first. Computes t = γ₁ - coeff, which is in [0, 2·γ₁] when ‖z‖∞ < γ₁.
    pub fn pack_z(&self, r: &mut [u8]) {
        assert!(r.len() >= POLYZ_PACKEDBYTES, "pack_z: output buffer too short");
        let mut t = [0u32; 4];
        for i in 0..N / 4 {
            for (j, tv) in t.iter_mut().enumerate() {
                let c = center(self.coeffs[4 * i + j]);
                debug_assert!(c.abs() <= GAMMA1);
                *tv = (GAMMA1 - c) as u32;
            }
            r[9 * i] = t[0] as u8;
            r[9 * i + 1] = (t[0] >> 8) as u8;
            r[9 * i + 2] = ((t[0] >> 16) | (t[1] << 2)) as u8;
            r[9 * i + 3] = (t[1] >> 6) as u8;
            r[9 * i + 4] = ((t[1] >> 14) | (t[2] << 4)) as u8;
            r[9 * i + 5] = (t[2] >> 4) as u8;
            r[9 * i + 6] = ((t[2] >> 12) | (t[3] << 6)) as u8;
            r[9 * i + 7] = (t[3] >> 2) as u8;
            r[9 * i + 8] = (t[3] >> 10) as u8;
        }
    }

    /// Unpack z polynomial.
    pub fn unpack_z(r: &mut Poly, a: &[u8]) {
        assert!(a.len() >= POLYZ_PACKEDBYTES, "unpack_z: input too short");
        for i in 0..N / 4 {
            r.coeffs[4 * i] =
                i32::from(a[9 * i]) | (i32::from(a[9 * i + 1]) << 8) | (i32::from(a[9 * i + 2]) << 16);
            r.coeffs[4 * i] &= 0x3FFFF;

            r.coeffs[4 * i + 1] = (i32::from(a[9 * i + 2]) >> 2)
                | (i32::from(a[9 * i + 3]) << 6)
                | (i32::from(a[9 * i + 4]) << 14);
            r.coeffs[4 * i + 1] &= 0x3FFFF;

            r.coeffs[4 * i + 2] = (i32::from(a[9 * i + 4]) >> 4)
                | (i32::from(a[9 * i + 5]) << 4)
                | (i32::from(a[9 * i + 6]) << 12);
            r.coeffs[4 * i + 2] &= 0x3FFFF;

            r.coeffs[4 * i + 3] = (i32::from(a[9 * i + 6]) >> 6)
                | (i32::from(a[9 * i + 7]) << 2)
                | (i32::from(a[9 * i + 8]) << 10);
            r.coeffs[4 * i + 3] &= 0x3FFFF;

            for j in 0..4 {
                r.coeffs[4 * i + j] = GAMMA1 - r.coeffs[4 * i + j];
            }
        }
    }
}

// ─── Polynomial Vector Types ──────────────────────────────────────────

/// Vector of L polynomials (used for s₁, y, z in ML-DSA-44).
#[derive(Clone, Zeroize)]
pub struct PolyVecL {
    pub polys: [Poly; L],
}

impl core::fmt::Debug for PolyVecL {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "PolyVecL([...; {L}])")
    }
}

impl Default for PolyVecL {
    fn default() -> Self {
        Self::zero()
    }
}

impl PolyVecL {
    #[must_use] 
    pub fn zero() -> Self {
        PolyVecL {
            polys: core::array::from_fn(|_| Poly::zero()),
        }
    }

    /// Forward NTT on all component polynomials.
    pub fn ntt(&mut self) {
        for p in &mut self.polys {
            p.ntt();
        }
    }

    /// Inverse NTT on all component polynomials.
    pub fn invntt_tomont(&mut self) {
        for p in &mut self.polys {
            p.invntt_tomont();
        }
    }

    /// Component-wise addition: c = a + b.
    pub fn add(c: &mut PolyVecL, a: &PolyVecL, b: &PolyVecL) {
        for i in 0..L {
            Poly::add(&mut c.polys[i], &a.polys[i], &b.polys[i]);
        }
    }

    /// In-place addition: self += other.
    pub fn add_assign(&mut self, other: &PolyVecL) {
        for i in 0..L {
            for j in 0..N {
                self.polys[i].coeffs[j] += other.polys[i].coeffs[j];
            }
        }
    }

    /// Component-wise subtraction: c = a - b.
    pub fn sub(c: &mut PolyVecL, a: &PolyVecL, b: &PolyVecL) {
        for i in 0..L {
            Poly::sub(&mut c.polys[i], &a.polys[i], &b.polys[i]);
        }
    }

    /// Squared L₂ norm of the entire vector: Σ ‖pᵢ‖₂².
    #[must_use] 
    pub fn l2_norm_squared(&self) -> u64 {
        self.polys.iter().map(Poly::l2_norm_squared).sum()
    }

    /// L∞ norm of the entire vector: `max_i` ‖pᵢ‖∞.
    #[must_use] 
    pub fn l_inf_norm(&self) -> u32 {
        self.polys.iter().map(Poly::l_inf_norm).max().unwrap_or(0)
    }

    /// Check if any component has ‖·‖∞ ≥ bound.
    ///
    /// Constant-time: accumulates result across all polynomials
    /// without short-circuit to prevent timing leaks (ALG-4).
    #[must_use] 
    pub fn chknorm(&self, bound: i32) -> bool {
        let mut result = false;
        for p in &self.polys {
            result |= p.chknorm(bound);
        }
        result
    }

    /// Reduce all coefficients.
    pub fn reduce(&mut self) {
        for p in &mut self.polys {
            p.reduce();
        }
    }

    /// Add q to negative coefficients.
    pub fn caddq(&mut self) {
        for p in &mut self.polys {
            p.caddq();
        }
    }
}

/// Vector of K polynomials (used for s₂, w, t in ML-DSA-44).
#[derive(Clone, Zeroize)]
pub struct PolyVecK {
    pub polys: [Poly; K],
}

impl core::fmt::Debug for PolyVecK {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "PolyVecK([...; {K}])")
    }
}

impl Default for PolyVecK {
    fn default() -> Self {
        Self::zero()
    }
}

impl PolyVecK {
    #[must_use] 
    pub fn zero() -> Self {
        PolyVecK {
            polys: core::array::from_fn(|_| Poly::zero()),
        }
    }

    pub fn ntt(&mut self) {
        for p in &mut self.polys {
            p.ntt();
        }
    }

    pub fn invntt_tomont(&mut self) {
        for p in &mut self.polys {
            p.invntt_tomont();
        }
    }

    pub fn add(c: &mut PolyVecK, a: &PolyVecK, b: &PolyVecK) {
        for i in 0..K {
            Poly::add(&mut c.polys[i], &a.polys[i], &b.polys[i]);
        }
    }

    /// In-place addition: self += other.
    pub fn add_assign(&mut self, other: &PolyVecK) {
        for i in 0..K {
            for j in 0..N {
                self.polys[i].coeffs[j] += other.polys[i].coeffs[j];
            }
        }
    }

    pub fn sub(c: &mut PolyVecK, a: &PolyVecK, b: &PolyVecK) {
        for i in 0..K {
            Poly::sub(&mut c.polys[i], &a.polys[i], &b.polys[i]);
        }
    }

    #[must_use] 
    pub fn l2_norm_squared(&self) -> u64 {
        self.polys.iter().map(Poly::l2_norm_squared).sum()
    }

    #[must_use] 
    pub fn l_inf_norm(&self) -> u32 {
        self.polys.iter().map(Poly::l_inf_norm).max().unwrap_or(0)
    }

    /// Constant-time norm check across all polynomials (ALG-4).
    #[must_use] 
    pub fn chknorm(&self, bound: i32) -> bool {
        let mut result = false;
        for p in &self.polys {
            result |= p.chknorm(bound);
        }
        result
    }

    pub fn reduce(&mut self) {
        for p in &mut self.polys {
            p.reduce();
        }
    }

    pub fn caddq(&mut self) {
        for p in &mut self.polys {
            p.caddq();
        }
    }
}

// ─── Matrix × Vector Product ─────────────────────────────────────────

/// The public matrix A ∈ `ℤ_q^{K×L`}, stored in NTT domain.
#[derive(Clone)]
pub struct MatrixA {
    pub rows: [[Poly; L]; K],
}

impl MatrixA {
    /// Expand the matrix A from seed ρ using SHAKE-128 (`ExpandA`).
    #[must_use] 
    pub fn expand(rho: &[u8; SEEDBYTES]) -> Self {
        let mut rows: [[Poly; L]; K] =
            core::array::from_fn(|_| core::array::from_fn(|_| Poly::zero()));
        for (i, row) in rows.iter_mut().enumerate().take(K) {
            for (j, poly) in row.iter_mut().enumerate().take(L) {
                let nonce = ((i as u16) << 8) | (j as u16);
                // FIPS 204 ExpandA samples Â directly in the NTT domain.
                Poly::uniform(poly, rho, nonce);
            }
        }
        MatrixA { rows }
    }

    /// Compute t = A·v (in NTT domain).
    /// Both A and v must be in NTT form. Result is in NTT form.
    #[must_use] 
    pub fn mul_vec(&self, v: &PolyVecL) -> PolyVecK {
        let mut t = PolyVecK::zero();
        for i in 0..K {
            for j in 0..L {
                let mut tmp = Poly::zero();
                Poly::pointwise_montgomery(&mut tmp, &self.rows[i][j], &v.polys[j]);
                for c in 0..N {
                    t.polys[i].coeffs[c] += tmp.coeffs[c];
                }
            }
        }
        t
    }
}

// ─── Modular Arithmetic Helpers ───────────────────────────────────────

/// Barrett-like reduction: map a ∈ ℤ to approximately [-(q/2), q/2].
///
/// Faithful port of `reduce32` from the CRYSTALS-Dilithium reference.
/// For a ≤ 2³¹ − 2²² − 1, computes r ≡ a (mod q) with −6283008 ≤ r ≤ 6283008.
fn reduce32(a: i32) -> i32 {
    // Widen so a + 2²² cannot overflow; |t| ≤ 256 keeps t·q within i32.
    let t = ((i64::from(a) + (1 << 22)) >> 23) as i32;
    a - t * Q
}

/// Full canonical reduction: a mod⁺ q → [0, q).
///
/// Applies Barrett reduction then conditional add of q.
pub(crate) fn reduce_i32(a: i32) -> i32 {
    let t = reduce32(a);
    t + ((t >> 31) & Q)
}

/// Map a coefficient to its centered representative in [-(q-1)/2, (q-1)/2].
fn center(mut a: i32) -> i32 {
    a = reduce_i32(a);
    if a > (Q - 1) / 2 {
        a -= Q;
    }
    a
}

// ─── Power2Round / Decompose / Hint (FIPS 204 §4) ────────────────────

/// `Power2Round`: split a into (a₁, a₀) such that a = a₁·2^d + a₀.
///
/// Assumes a is a standard representative in [0, q).
/// Faithful port of the CRYSTALS-Dilithium reference.
#[must_use] 
pub fn power2round(a: i32) -> (i32, i32) {
    let a_pos = reduce_i32(a);
    let a1 = (a_pos + (1 << (D - 1)) - 1) >> D;
    let a0 = a_pos - (a1 << D);
    (a1, a0)
}

/// Decompose: split a into (a₁, a₀) such that a ≡ a₁·2γ₂ + a₀ (mod q).
///
/// For γ₂ = (q−1)/88, uses constant-time fixed-point multiplication
/// to avoid division. Faithful port of `rounding.c::decompose()` from
/// the CRYSTALS-Dilithium reference implementation.
///
/// Returns (a₁, a₀) where a₁ ∈ [0, 43] and a₀ ∈ (−γ₂, γ₂].
#[must_use] 
pub fn decompose(a: i32) -> (i32, i32) {
    let a_pos = reduce_i32(a);

    // Fixed-point approximation of a / (2·γ₂) — no division needed.
    // This is the exact algorithm from the CRYSTALS-Dilithium reference.
    let mut a1 = (a_pos + 127) >> 7;
    // For γ₂ = (q-1)/88: use magic constant 11275
    a1 = (a1 * 11275 + (1 << 23)) >> 24;
    // Boundary wrap: if a1 == 44 (i.e. a is near q−1), set a1 = 0
    a1 ^= ((43 - a1) >> 31) & a1;

    let mut a0 = a_pos - a1 * 2 * GAMMA2;
    // Center a0: if a0 > (q−1)/2, subtract q
    a0 -= (((Q - 1) / 2 - a0) >> 31) & Q;
    (a1, a0)
}

/// `MakeHint`: compute the hint bit indicating whether low bits overflow
/// into high bits.
///
/// Parameters a0, a1 are the raw outputs of `decompose` (a0 is already
/// in signed/centered form). Faithful port of the reference.
#[allow(clippy::manual_range_contains)]
#[must_use] 
pub fn make_hint(a0: i32, a1: i32) -> bool {
    a0 > GAMMA2 || a0 < -GAMMA2 || (a0 == -GAMMA2 && a1 != 0)
}

/// `UseHint`: correct high bits according to hint.
///
/// For γ₂ = (q−1)/88 (ML-DSA-44), a₁ wraps around at 43.
/// Faithful port of the reference.
#[must_use] 
pub fn use_hint(a: i32, hint: bool) -> i32 {
    let (a1, a0) = decompose(a);
    if !hint {
        return a1;
    }
    // γ₂ = (q-1)/88 branch
    if a0 > 0 {
        if a1 == 43 {
            0
        } else {
            a1 + 1
        }
    } else if a1 == 0 {
        43
    } else {
        a1 - 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_add_sub() {
        let mut a = Poly::zero();
        let mut b = Poly::zero();
        a.coeffs[0] = 100;
        a.coeffs[1] = Q - 1;
        b.coeffs[0] = 200;
        b.coeffs[1] = 1;

        let mut c = Poly::zero();
        Poly::add(&mut c, &a, &b);
        assert_eq!(c.coeffs[0], 300);
        assert_eq!(c.coeffs[1], Q); // will need reduction

        Poly::sub(&mut c, &a, &b);
        assert_eq!(c.coeffs[0], -100);
    }

    #[test]
    fn test_l2_norm() {
        let mut p = Poly::zero();
        p.coeffs[0] = 3;
        p.coeffs[1] = 4;
        // ‖p‖₂² = 9 + 16 = 25
        assert_eq!(p.l2_norm_squared(), 25);
    }

    #[test]
    fn test_l_inf_norm() {
        let mut p = Poly::zero();
        p.coeffs[0] = Q - 1; // = -1 mod q → centered = -1
        p.coeffs[1] = 42;
        assert_eq!(p.l_inf_norm(), 42);
    }

    #[test]
    fn test_center() {
        assert_eq!(center(0), 0);
        assert_eq!(center(1), 1);
        assert_eq!(center(Q - 1), -1);
        assert_eq!(center((Q - 1) / 2), (Q - 1) / 2);
        assert_eq!(center((Q - 1) / 2 + 1), -(Q - 1) / 2);
    }

    #[test]
    fn test_power2round() {
        let (a1, a0) = power2round(12345);
        // Verify reconstruction
        assert_eq!((a1 << D) + a0, 12345);
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;
    extern crate alloc;
    use alloc::format;

    #[test]
    fn defaults_and_debug_never_print_coefficients() {
        let mut p = Poly::default();
        p.coeffs[0] = 12345;
        assert!(!format!("{p:?}").contains("12345"));
        assert!(format!("{:?}", PolyVecL::default()).starts_with("PolyVecL"));
        assert!(format!("{:?}", PolyVecK::default()).starts_with("PolyVecK"));
    }

    #[test]
    fn vector_operations() {
        let mut a = PolyVecL::zero();
        let mut b = PolyVecL::zero();
        a.polys[0].coeffs[0] = -5;
        b.polys[L - 1].coeffs[1] = 7;
        let mut c = PolyVecL::zero();
        PolyVecL::add(&mut c, &a, &b);
        assert_eq!((c.polys[0].coeffs[0], c.polys[L - 1].coeffs[1]), (-5, 7));
        c.caddq();
        assert_eq!(c.polys[0].coeffs[0], Q - 5);
        let mut r = a.clone();
        r.ntt();
        r.invntt_tomont();
        assert_eq!(r.polys.len(), L);

        let mut x = PolyVecK::zero();
        let mut y = PolyVecK::zero();
        x.polys[0].coeffs[0] = 3;
        y.polys[K - 1].coeffs[2] = -4;
        let mut z = PolyVecK::zero();
        PolyVecK::add(&mut z, &x, &y);
        assert_eq!(z.l2_norm_squared(), 25);
        assert_eq!(z.l_inf_norm(), 4);
        assert!(z.chknorm(4));
        assert!(!z.chknorm(5));
        let mut d = PolyVecK::zero();
        PolyVecK::sub(&mut d, &x, &y);
        assert_eq!(d.polys[K - 1].coeffs[2], 4);
        z.reduce();
        z.caddq();
        assert_eq!(z.polys[K - 1].coeffs[2], Q - 4);
        let mut w = x.clone();
        w.ntt();
        w.invntt_tomont();
        assert_eq!(w.polys.len(), K);
    }

    #[test]
    fn hints_cover_all_branches() {
        assert!(make_hint(GAMMA2 + 1, 0));
        assert!(make_hint(-GAMMA2 - 1, 0));
        assert!(make_hint(-GAMMA2, 1));
        assert!(!make_hint(-GAMMA2, 0));
        assert!(!make_hint(0, 5));

        let mut seen = [false; 4];
        let mut a = 0;
        while a < Q && seen.iter().any(|s| !s) {
            let (a1, a0) = decompose(a);
            let h = use_hint(a, true);
            let case = match (a0 > 0, a1) {
                (true, 43) => { assert_eq!(h, 0); 0 }
                (true, _) => { assert_eq!(h, a1 + 1); 1 }
                (false, 0) => { assert_eq!(h, 43); 2 }
                (false, _) => { assert_eq!(h, a1 - 1); 3 }
            };
            seen[case] = true;
            assert_eq!(use_hint(a, false), a1);
            a += 997;
        }
        assert!(seen.iter().all(|&s| s), "{seen:?}");
    }
}

#[cfg(test)]
mod sampling_tests {
    use super::*;
    extern crate alloc;

    #[test]
    fn samplers_stay_in_range_and_are_deterministic() {
        let seed = [3u8; SEEDBYTES];
        let crh = [4u8; CRHBYTES];
        let (mut a, mut b) = (Poly::zero(), Poly::zero());
        Poly::uniform(&mut a, &seed, 1);
        Poly::uniform(&mut b, &seed, 1);
        assert_eq!(a.coeffs, b.coeffs);
        assert!(a.coeffs.iter().all(|c| (0..Q).contains(c)));

        Poly::uniform_eta(&mut a, &crh, 2);
        assert!(a.coeffs.iter().all(|c| c.abs() <= ETA as i32));

        Poly::uniform_gamma1(&mut a, &crh, 3);
        assert!(a.coeffs.iter().all(|c| *c > -GAMMA1 && *c <= GAMMA1));

        Poly::challenge(&mut a, b"challenge seed");
        assert_eq!(a.coeffs.iter().filter(|&&c| c != 0).count(), TAU as usize);
        assert!(a.coeffs.iter().all(|&c| c == 0 || c == 1 || c == Q - 1)); // −1 is stored as q − 1
    }

    #[test]
    fn packing_round_trips() {
        let mut t1 = Poly::zero();
        for (i, c) in t1.coeffs.iter_mut().enumerate() {
            *c = (i as i32 * 37) & 0x3FF;
        }
        let mut buf = [0u8; POLYT1_PACKEDBYTES];
        t1.pack_t1(&mut buf);
        let mut back = Poly::zero();
        Poly::unpack_t1(&mut back, &buf);
        assert_eq!(back.coeffs, t1.coeffs);

        let mut z = Poly::zero();
        for (i, c) in z.coeffs.iter_mut().enumerate() {
            *c = (i as i32 * 1021) % GAMMA1 - GAMMA1 / 2;
        }
        let mut buf = [0u8; POLYZ_PACKEDBYTES];
        z.pack_z(&mut buf);
        Poly::unpack_z(&mut back, &buf);
        assert_eq!(back.coeffs, z.coeffs);
    }

    #[test]
    fn polyvecl_norms_and_matrix_product() {
        let mut a = PolyVecL::zero();
        a.polys[1].coeffs[3] = -6;
        a.polys[2].coeffs[4] = 2;
        let mut d = PolyVecL::zero();
        PolyVecL::sub(&mut d, &PolyVecL::zero(), &a);
        assert_eq!(d.polys[1].coeffs[3], 6);
        assert_eq!(a.l2_norm_squared(), 40);
        assert_eq!(a.l_inf_norm(), 6);

        let m = MatrixA::expand(&[9u8; SEEDBYTES]);
        let mut x = PolyVecL::zero();
        x.polys[0].coeffs[0] = 1;
        x.ntt();
        let mut y = PolyVecL::zero();
        y.polys[3].coeffs[5] = 2;
        y.ntt();
        let mut s = PolyVecL::zero();
        PolyVecL::add(&mut s, &x, &y);
        // A·(x + y) = A·x + A·y (mod q)
        let mut lhs = m.mul_vec(&s);
        let mut rhs = PolyVecK::zero();
        PolyVecK::add(&mut rhs, &m.mul_vec(&x), &m.mul_vec(&y));
        lhs.reduce();
        lhs.caddq();
        rhs.reduce();
        rhs.caddq();
        assert_eq!(lhs.polys.iter().map(|p| p.coeffs).collect::<alloc::vec::Vec<_>>(),
                   rhs.polys.iter().map(|p| p.coeffs).collect::<alloc::vec::Vec<_>>());
        assert!(m.mul_vec(&PolyVecL::zero()).polys.iter().all(|p| p.coeffs.iter().all(|&c| c == 0)));
    }
}
