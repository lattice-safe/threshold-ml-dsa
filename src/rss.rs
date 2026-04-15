//! Replicated Secret Sharing (RSS) key generation for threshold ML-DSA.
//!
//! Implements the paper's Figure 4 (Keygen): fresh independent secrets
//! are sampled per (N-T+1)-subset, and the public key is derived from
//! their sum. This matches the Go reference implementation exactly.
//!
//! ## Key Differences from v0.2 (Security-Critical)
//!
//! v0.2 decomposed an existing ML-DSA key via additive sharing with a
//! remainder piece. This violated the paper's assumption that **every**
//! subset secret `s_I` is independently drawn from `χ_s = U([-η, η])`.
//!
//! v0.3 follows the paper exactly:
//! 1. Each subset secret `s_I` is independently sampled from `χ_s`
//! 2. The composite secret `s = Σ s_I` is computed as the sum
//! 3. The public key `t = [A|I]·s` is derived from the sum
//! 4. No existing key is decomposed — keys are always fresh
//!
//! ## Subset Encoding
//!
//! Subsets are encoded as bitmasks: bit `i` is set if party `i` is a
//! member. This matches the Go reference. For example, with N=4:
//! - `{0,1,2}` → bitmask `0b0111` = 7
//! - `{1,3}`   → bitmask `0b1010` = 10

#[cfg(not(feature = "std"))]
use alloc::{collections::BTreeMap, vec::Vec};
#[cfg(feature = "std")]
use std::collections::BTreeMap;

use crate::params::*;
use crate::poly::{Poly, PolyVecK, PolyVecL};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use zeroize::Zeroize;

/// A single subset's secret key: `(s₁_I, s₂_I)` where I ⊂ [N], |I| = N-T+1.
/// Both vectors have coefficients drawn from `χ_s = U([-η, η])`.
#[derive(Clone)]
pub struct Share {
    pub s1: PolyVecL,
    pub s2: PolyVecK,
    /// NTT-domain cache of s₁
    pub s1h: PolyVecL,
    /// NTT-domain cache of s₂
    pub s2h: PolyVecK,
}

impl Share {
    fn zero() -> Self {
        Self {
            s1: PolyVecL::zero(),
            s2: PolyVecK::zero(),
            s1h: PolyVecL::zero(),
            s2h: PolyVecK::zero(),
        }
    }
}

// Implement Zeroize manually since PolyVecL/K already implement it
impl Zeroize for Share {
    fn zeroize(&mut self) {
        self.s1.zeroize();
        self.s2.zeroize();
        self.s1h.zeroize();
        self.s2h.zeroize();
    }
}

impl Drop for Share {
    fn drop(&mut self) {
        self.zeroize();
    }
}

/// A party's private key material in the threshold scheme.
///
/// Contains this party's ID, the shared public data (ρ, tr),
/// per-party randomness `key`, and a map from subset bitmask
/// to the corresponding `Share`.
#[derive(Clone)]
pub struct ThresholdPrivateKey {
    /// This party's ID ∈ [0, N)
    pub id: u8,
    /// ρ: the public seed for A = ExpandA(ρ)
    pub rho: [u8; 32],
    /// Per-party random key for hedged nonce derivation
    pub key: [u8; 32],
    /// tr = CRH(pk) — hash of the public key
    pub tr: [u8; TRBYTES],
    /// Map from subset bitmask → Share
    /// (party holds all shares for subsets it belongs to)
    pub shares: BTreeMap<u8, Share>,
}

impl Zeroize for ThresholdPrivateKey {
    fn zeroize(&mut self) {
        self.key.zeroize();
        self.tr.zeroize();
        for (_, share) in self.shares.iter_mut() {
            share.zeroize();
        }
    }
}

impl Drop for ThresholdPrivateKey {
    fn drop(&mut self) {
        self.zeroize();
    }
}

impl core::fmt::Debug for ThresholdPrivateKey {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "ThresholdPrivateKey {{ id: {}, shares: {} }}",
            self.id,
            self.shares.len()
        )
    }
}

/// Threshold key generation following ePrint 2026/013, Figure 4.
///
/// Generates a fresh threshold ML-DSA-44 keypair from a 32-byte seed.
/// This is a deterministic function — the same seed always produces
/// the same keys.
///
/// # Process (matching Go reference `NewThresholdKeysFromSeed`)
///
/// 1. Derive ρ from seed via SHAKE-256
/// 2. Compute A = ExpandA(ρ)
/// 3. For each (N-T+1)-subset I (enumerated via Gosper's hack):
///    a. Derive s_I = (s1_I, s2_I) ← χ_s independently
///    b. Distribute s_I to all parties in I
/// 4. Compute composite s = Σ s_I
/// 5. Compute t = As₁ + s₂, Power2Round → (t₀, t₁)
/// 6. Public key = (ρ, t₁), tr = CRH(pk)
///
/// # Returns
/// `(public_key_bytes, Vec<ThresholdPrivateKey>)` — the packed public key
/// and one private key per party.
pub fn keygen_from_seed(
    seed: &[u8; 32],
    params: &ThresholdParams,
) -> ([u8; PK_BYTES], Vec<ThresholdPrivateKey>) {
    let n = params.n;
    let t = params.t;
    assert!(
        t >= 2 && t <= n,
        "invalid threshold parameters: t={}, n={}",
        t,
        n
    );
    assert!(
        (n as usize) <= MAX_PARTIES,
        "unsupported party count n={}, max={}",
        n,
        MAX_PARTIES
    );

    // Expand seed via SHAKE-256
    let mut h = Shake256::default();
    h.update(seed);
    // Match Go: if NIST mode, also absorb (K, L)
    // For ML-DSA-44 (FIPS mode), we skip this
    let mut reader = h.finalize_xof();

    // 1. Derive ρ (32 bytes)
    let mut rho = [0u8; 32];
    reader.read(&mut rho);

    // 2. Initialize private keys
    let mut sks: Vec<ThresholdPrivateKey> = Vec::with_capacity(n as usize);
    for i in 0..n {
        let mut key = [0u8; 32];
        reader.read(&mut key);
        sks.push(ThresholdPrivateKey {
            id: i,
            rho,
            key,
            tr: [0u8; TRBYTES],
            shares: BTreeMap::new(),
        });
    }

    // 3. Accumulator for the composite secret s = Σ s_I
    let mut s1_total = PolyVecL::zero();
    let mut s2_total = PolyVecK::zero();
    let mut s1h_total = PolyVecL::zero();
    let mut s2h_total = PolyVecK::zero();

    // Enumerate all (N-T+1)-subsets using Gosper's hack on bitmasks.
    // Start with the lexicographically smallest bitmask with (N-T+1) bits set.
    let subset_size = n - t + 1;
    let mut mask: u16 = (1u16 << subset_size) - 1; // e.g., for size=2: 0b11
    let limit: u16 = 1u16 << n;

    while mask < limit {
        let subset_mask = mask as u8;

        // Sample a fresh independent secret s_I from χ_s
        let mut sseed = [0u8; 64];
        reader.read(&mut sseed);

        let mut share = Share::zero();

        // Sample s1_I ← U([-η, η]^{L×N})
        for j in 0..L {
            sample_leq_eta(&mut share.s1.polys[j], &sseed, j as u16);
        }
        // Sample s2_I ← U([-η, η]^{K×N})
        for j in 0..K {
            sample_leq_eta(&mut share.s2.polys[j], &sseed, (j + L) as u16);
        }

        // Pre-compute NTT forms
        share.s1h = share.s1.clone();
        share.s1h.ntt();
        share.s2h = share.s2.clone();
        share.s2h.ntt();

        // 4. Distribute to all parties in this subset
        for i in 0..n {
            if subset_mask & (1 << i) != 0 {
                sks[i as usize].shares.insert(subset_mask, share.clone());
            }
        }

        // 5. Accumulate into composite secret
        s1_total.add_assign(&share.s1);
        s1h_total.add_assign(&share.s1h);
        s2_total.add_assign(&share.s2);
        s2h_total.add_assign(&share.s2h);

        // Gosper's hack: advance to next subset bitmask with same popcount
        let c = mask & mask.wrapping_neg();
        if c == 0 {
            break;
        }
        let r = mask.wrapping_add(c);
        if r >= limit {
            break;
        }
        mask = (((r ^ mask) >> 2) / c) | r;
    }

    // 6. Normalize the accumulated sums
    s1_total.reduce();
    s1h_total.reduce();
    s2_total.reduce();
    s2h_total.reduce();

    // 7. Compute public key: t = A·s₁ + s₂, Power2Round → (t₀, t₁)
    let pk_bytes = compute_public_key(&rho, &s1h_total, &s2_total);

    // 8. Compute tr = CRH(pk)
    let mut tr = [0u8; TRBYTES];
    let mut h2 = Shake256::default();
    h2.update(&pk_bytes);
    let mut tr_reader = h2.finalize_xof();
    tr_reader.read(&mut tr);

    // Set tr in all private keys
    for sk in sks.iter_mut() {
        sk.tr = tr;
    }

    (pk_bytes, sks)
}

/// Sample a polynomial with coefficients uniformly in `[-η, η]`.
///
/// Uses SHAKE-256 rejection sampling (matching the Go reference
/// `PolyDeriveUniformLeqEta`).
fn sample_leq_eta(p: &mut Poly, seed: &[u8; 64], nonce: u16) {
    let mut h = Shake256::default();
    let mut iv = [0u8; 66];
    iv[..64].copy_from_slice(seed);
    iv[64] = nonce as u8;
    iv[65] = (nonce >> 8) as u8;
    h.update(&iv);

    let mut reader = h.finalize_xof();
    let mut buf = [0u8; 136]; // SHAKE-256 rate
    let mut i = 0;

    while i < N {
        reader.read(&mut buf);
        for byte in &buf {
            if i >= N {
                break;
            }
            let t1 = (*byte & 0x0F) as u32;
            let t2 = (*byte >> 4) as u32;

            // η=2: accept if < 15, reduce mod 5
            if ETA == 2 {
                if t1 <= 14 {
                    let reduced = t1 - ((205 * t1) >> 10) * 5;
                    p.coeffs[i] = Q + ETA as i32 - reduced as i32;
                    i += 1;
                }
                if t2 <= 14 && i < N {
                    let reduced = t2 - ((205 * t2) >> 10) * 5;
                    p.coeffs[i] = Q + ETA as i32 - reduced as i32;
                    i += 1;
                }
            }
        }
    }
}

/// Compute the public key bytes: `pk = (ρ ‖ t₁_packed)`.
///
/// Internally: t = A·NTT⁻¹(s₁_hat) + s₂, Power2Round(t) → (t₀, t₁)
fn compute_public_key(rho: &[u8; 32], s1h_total: &PolyVecL, s2_total: &PolyVecK) -> [u8; PK_BYTES] {
    use dilithium::{
        packing::pack_pk,
        poly::Poly as DPoly,
        polyvec::{
            matrix_expand, matrix_pointwise_montgomery, polyveck_add, polyveck_caddq,
            polyveck_reduce, PolyVecK as DPolyVecK, PolyVecL as DPolyVecL,
        },
        ML_DSA_44,
    };

    let mode = ML_DSA_44;

    // Expand A from ρ
    let mut mat: [DPolyVecL; K] = core::array::from_fn(|_| DPolyVecL::default());
    matrix_expand(mode, &mut mat, rho);

    // Convert our s1h to dilithium-rs PolyVecL
    let mut s1h_d = DPolyVecL::default();
    for i in 0..L {
        for j in 0..N {
            s1h_d.vec[i].coeffs[j] = s1h_total.polys[i].coeffs[j];
        }
    }

    // t = A·s₁ (in NTT domain, then InvNTT + Montgomery)
    let mut t = DPolyVecK::default();
    matrix_pointwise_montgomery(mode, &mut t, &mat, &s1h_d);

    // InvNTT each component of t
    for i in 0..K {
        t.vec[i].invntt_tomont();
    }

    // t += s₂
    let mut s2_d = DPolyVecK::default();
    for i in 0..K {
        for j in 0..N {
            s2_d.vec[i].coeffs[j] = s2_total.polys[i].coeffs[j];
        }
    }
    let t_copy = t.clone();
    polyveck_add(mode, &mut t, &t_copy, &s2_d);
    // Threshold keygen accumulates subset secrets modulo q; after adding s2,
    // we must perform a full reduction before Power2Round.
    polyveck_reduce(mode, &mut t);
    polyveck_caddq(mode, &mut t);

    // Power2Round(t) → (t0, t1)
    let mut t0 = DPolyVecK::default();
    let mut t1 = DPolyVecK::default();
    for i in 0..K {
        DPoly::power2round(&mut t1.vec[i], &mut t0.vec[i], &t.vec[i]);
    }

    // Pack: pk = ρ ‖ t₁
    let mut pk = [0u8; PK_BYTES];
    pack_pk(mode, &mut pk, rho, &t1);

    pk
}

/// Legacy subset enumeration (preserved for backward compatibility).
///
/// Returns a list of subsets, where each subset is a sorted Vec of party indices.
pub fn enumerate_subsets(n: usize, t: usize) -> Vec<Vec<usize>> {
    if n == 0 || t == 0 || t > n || n < 2 {
        return Vec::new();
    }
    let subset_size = n - t + 1;
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(subset_size);
    enumerate_subsets_helper(n, subset_size, 0, &mut current, &mut result);
    result
}

fn enumerate_subsets_helper(
    n: usize,
    t: usize,
    start: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == t {
        result.push(current.clone());
        return;
    }
    let remaining = t - current.len();
    for i in start..=(n - remaining) {
        current.push(i);
        enumerate_subsets_helper(n, t, i + 1, current, result);
        current.pop();
    }
}

// ─── Legacy types for backward compatibility ───────────────────────────

/// Legacy party key share type (v0.2 API). Superseded by ThresholdPrivateKey.
#[derive(Clone, Zeroize)]
#[zeroize(drop)]
pub struct PartyKeyShare {
    #[zeroize(skip)]
    pub party_id: usize,
    #[zeroize(skip)]
    pub n: usize,
    #[zeroize(skip)]
    pub t: usize,
    pub s1_shares: Vec<PolyVecL>,
    pub s2_shares: Vec<PolyVecK>,
    #[zeroize(skip)]
    pub subset_indices: Vec<usize>,
}

impl core::fmt::Debug for PartyKeyShare {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "PartyKeyShare {{ party_id: {}, n: {}, t: {}, shares: {} }}",
            self.party_id,
            self.n,
            self.t,
            self.s1_shares.len()
        )
    }
}

impl PartyKeyShare {
    /// Compute this party's "effective" s₁ share by summing all share pieces.
    pub fn effective_s1(&self) -> PolyVecL {
        let mut result = PolyVecL::zero();
        for share in &self.s1_shares {
            result.add_assign(share);
        }
        result
    }

    /// Compute this party's "effective" s₂ share by summing all share pieces.
    pub fn effective_s2(&self) -> PolyVecK {
        let mut result = PolyVecK::zero();
        for share in &self.s2_shares {
            result.add_assign(share);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_keygen_produces_valid_keys() {
        let seed = [42u8; 32];
        let params = get_threshold_params(2, 3).unwrap();
        let (pk, sks) = keygen_from_seed(&seed, &params);

        // Should produce N=3 private keys
        assert_eq!(sks.len(), 3);

        // Each party holds shares for subsets it belongs to
        for sk in &sks {
            assert!(!sk.shares.is_empty());
            // Every share bitmask should include this party
            for (&mask, _) in &sk.shares {
                assert!(
                    mask & (1 << sk.id) != 0,
                    "Party {} has share for subset 0b{:06b} which doesn't include it",
                    sk.id,
                    mask
                );
            }
        }

        // tr should be non-zero and identical across all parties
        assert_ne!(sks[0].tr, [0u8; TRBYTES]);
        assert_eq!(sks[0].tr, sks[1].tr);
        assert_eq!(sks[1].tr, sks[2].tr);

        // Public key should start with ρ
        assert_eq!(&pk[..32], &sks[0].rho);
    }

    #[test]
    fn test_keygen_deterministic() {
        let seed = [7u8; 32];
        let params = get_threshold_params(2, 2).unwrap();
        let (pk1, _) = keygen_from_seed(&seed, &params);
        let (pk2, _) = keygen_from_seed(&seed, &params);
        assert_eq!(pk1, pk2);
    }

    #[test]
    fn test_keygen_all_configs() {
        let configs: &[(u8, u8)] = &[(2, 2), (2, 3), (3, 3), (2, 4), (3, 4), (4, 4)];
        for &(t, n) in configs {
            let seed = [t + n; 32];
            let params = get_threshold_params(t, n).unwrap();
            let (_, sks) = keygen_from_seed(&seed, &params);
            assert_eq!(sks.len(), n as usize, "Wrong key count for ({}, {})", t, n);
        }
    }

    #[test]
    fn test_subset_count_matches_binomial() {
        // For (T=2, N=3): C(3,2) = 3 subsets
        let seed = [1u8; 32];
        let params = get_threshold_params(2, 3).unwrap();
        let (_, sks) = keygen_from_seed(&seed, &params);

        // Count total unique subset bitmasks
        let mut all_masks: Vec<u8> = Vec::new();
        for sk in &sks {
            for &mask in sk.shares.keys() {
                if !all_masks.contains(&mask) {
                    all_masks.push(mask);
                }
            }
        }
        assert_eq!(all_masks.len(), num_subsets(3, 2));
    }
}
