//! Replicated Secret Sharing (RSS) key generation for threshold ML-DSA.
//!
//! Implements the paper's keygen (Fig. 5, USENIX version): fresh independent secrets
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
//! 4. No existing key is decomposed — keys are always fresh (see
//!    [`crate::aposteriori`] for sharing an existing key)
//!
//! ## Subset Encoding
//!
//! Subsets are encoded as bitmasks: bit `i` is set if party `i` is a
//! member. This matches the Go reference. For example, with N=4:
//! - `{0,1,2}` → bitmask `0b0111` = 7
//! - `{1,3}`   → bitmask `0b1010` = 10

#[cfg(not(feature = "std"))]
use alloc::{collections::BTreeMap, vec::Vec};
use zeroize::Zeroizing;
#[cfg(feature = "std")]
use std::collections::BTreeMap;

use crate::error::Error;
use crate::params::{get_threshold_params, TRBYTES, ThresholdParams, PK_BYTES, MAX_PARTIES, L, K, N, ETA};
use crate::poly::{Poly, PolyVecK, PolyVecL};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use zeroize::Zeroize;
use crate::aposteriori::SHARE_BOUND;

/// A single subset's secret key: `(s₁_I, s₂_I)` where I ⊂ [N], |I| = N-T+1.
/// Fresh and DKG keys draw coefficients from `χ_s = U([-η, η])`; a
/// posteriori sharing ([`crate::aposteriori`]) gives Gaussian coefficients
/// in `[-10, 10]`.
///
/// Fields are `pub(crate)` to prevent accidental exfiltration.
#[derive(Clone)]
pub struct Share {
    pub(crate) s1: PolyVecL,
    pub(crate) s2: PolyVecK,
    /// NTT-domain cache of s₁
    pub(crate) s1h: PolyVecL,
    /// NTT-domain cache of s₂
    pub(crate) s2h: PolyVecK,
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
///
/// Fields are `pub(crate)` to prevent accidental exfiltration
/// via logging, serialization, or debug plumbing. Not `Clone`, so secret
/// shares are not silently duplicated.
pub struct ThresholdPrivateKey {
    /// Threshold T this key was generated for
    pub(crate) t: u8,
    /// Number of parties N this key was generated for
    pub(crate) n: u8,
    /// This party's ID ∈ [0, N)
    pub(crate) id: u8,
    /// ρ: the public seed for A = ExpandA(ρ)
    pub(crate) rho: [u8; 32],
    /// Per-party random key for hedged nonce derivation
    pub(crate) key: [u8; 32],
    /// tr = CRH(pk) — hash of the public key
    pub(crate) tr: [u8; TRBYTES],
    /// Map from subset bitmask → Share
    /// (party holds all shares for subsets it belongs to)
    pub(crate) shares: BTreeMap<u8, Share>,
}

impl ThresholdPrivateKey {
    /// Returns this party's ID.
    #[must_use]
    pub fn id(&self) -> u8 {
        self.id
    }

    /// Returns the public key hash tr = CRH(pk).
    #[must_use]
    pub fn tr(&self) -> &[u8; TRBYTES] {
        &self.tr
    }

    /// Returns an iterator over subset bitmasks this party holds.
    pub fn share_masks(&self) -> impl Iterator<Item = &u8> {
        self.shares.keys()
    }

    /// Identifier of this party key for a [`crate::nonce::NonceGuard`]:
    /// `SHAKE256("th-ml-dsa-nonce-guard-key-v1" ‖ tr ‖ id)[..32]`.
    #[must_use]
    pub fn guard_key_id(&self) -> [u8; 32] {
        let mut h = Shake256::default();
        h.update(b"th-ml-dsa-nonce-guard-key-v1");
        h.update(&self.tr);
        h.update(&[self.id]);
        let mut out = [0u8; 32];
        h.finalize_xof().read(&mut out);
        out
    }

    /// Returns the (T, N) configuration this key belongs to.
    #[must_use]
    pub fn threshold(&self) -> (u8, u8) {
        (self.t, self.n)
    }

    /// Serialize this party's key for storage or transfer to the party.
    ///
    /// Format (version 2): `2 ‖ t ‖ n ‖ id ‖ ρ ‖ key ‖ tr ‖ count ‖
    /// (mask ‖ s₁ ‖ s₂)*`, each share coefficient as one byte `c + 10`
    /// (fresh keys have `|c| ≤ η`, a posteriori shares `|c| ≤ 10`).
    /// The output contains secret material and is wiped on drop; protect it
    /// at rest and in transit.
    pub fn to_bytes(&self) -> Result<Zeroizing<Vec<u8>>, Error> {
        let mut out = Zeroizing::new(Vec::with_capacity(
            KEY_HEADER_BYTES + self.shares.len() * KEY_SHARE_BYTES,
        ));
        out.extend_from_slice(&[KEY_FORMAT_VERSION, self.t, self.n, self.id]);
        out.extend_from_slice(&self.rho);
        out.extend_from_slice(&self.key);
        out.extend_from_slice(&self.tr);
        out.push(u8::try_from(self.shares.len()).map_err(|_| Error::InvalidShare)?);
        for (&mask, share) in &self.shares {
            out.push(mask);
            for p in share.s1.polys.iter().chain(share.s2.polys.iter()) {
                for &c in &p.coeffs {
                    if !(-SHARE_BOUND..=SHARE_BOUND).contains(&c) {
                        return Err(Error::InvalidShare);
                    }
                    out.push((c + SHARE_BOUND) as u8);
                }
            }
        }
        Ok(out)
    }

    /// Parse a key produced by [`ThresholdPrivateKey::to_bytes`].
    ///
    /// Rejects unknown versions, unsupported (T, N), wrong lengths, share
    /// sets that are not exactly the (N−T+1)-subsets containing this party,
    /// and coefficients outside `[-10, 10]`.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() < KEY_HEADER_BYTES || bytes[0] != KEY_FORMAT_VERSION {
            return Err(Error::InvalidParameters);
        }
        let (t, n, id) = (bytes[1], bytes[2], bytes[3]);
        get_threshold_params(t, n).ok_or(Error::InvalidParameters)?;
        if id >= n {
            return Err(Error::InvalidParameters);
        }
        let mut off = 4;
        let mut take = |len: usize| {
            let s = &bytes[off..off + len];
            off += len;
            s
        };
        let mut sk = ThresholdPrivateKey {
            t,
            n,
            id,
            rho: take(32).try_into().map_err(|_| Error::InvalidParameters)?,
            key: take(32).try_into().map_err(|_| Error::InvalidParameters)?,
            tr: take(TRBYTES).try_into().map_err(|_| Error::InvalidParameters)?,
            shares: BTreeMap::new(),
        };
        let count = take(1)[0] as usize;
        let subset_size = (n - t + 1) as u32;
        let expected: Vec<u8> = (0u16..(1 << n))
            .map(|m| m as u8)
            .filter(|m| m.count_ones() == subset_size && m & (1 << id) != 0)
            .collect();
        if count != expected.len() || bytes.len() != KEY_HEADER_BYTES + count * KEY_SHARE_BYTES {
            return Err(Error::InvalidParameters);
        }
        for &want in &expected {
            let mask = take(1)[0];
            if mask != want {
                return Err(Error::InvalidShare);
            }
            let mut share = Share::zero();
            for p in share.s1.polys.iter_mut().chain(share.s2.polys.iter_mut()) {
                for (c, &b) in p.coeffs.iter_mut().zip(take(N)) {
                    if i32::from(b) > 2 * SHARE_BOUND {
                        return Err(Error::InvalidShare);
                    }
                    *c = i32::from(b) - SHARE_BOUND;
                }
            }
            share.s1h = share.s1.clone();
            share.s1h.reduce();
            share.s1h.ntt();
            share.s2h = share.s2.clone();
            share.s2h.reduce();
            share.s2h.ntt();
            sk.shares.insert(mask, share);
        }
        Ok(sk)
    }
}

const KEY_FORMAT_VERSION: u8 = 2;
const KEY_HEADER_BYTES: usize = 4 + 32 + 32 + TRBYTES + 1;
const KEY_SHARE_BYTES: usize = 1 + (L + K) * N;

impl Zeroize for ThresholdPrivateKey {
    fn zeroize(&mut self) {
        self.key.zeroize();
        self.tr.zeroize();
        for share in self.shares.values_mut() {
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

/// Threshold key generation following ePrint 2026/013 (Fig. 5 in the USENIX version).
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
///    a. Derive `s_I` = (`s1_I`, `s2_I`) ← `χ_s` independently
///    b. Distribute `s_I` to all parties in I
/// 4. Compute composite s = Σ `s_I`
/// 5. Compute t = As₁ + s₂, `Power2Round` → (t₀, t₁)
/// 6. Public key = (ρ, t₁), tr = CRH(pk)
///
/// # Returns
/// `Ok((public_key_bytes, Vec<ThresholdPrivateKey>))` — the packed public key
/// and one private key per party.
///
/// Returns `Err(Error::InvalidParameters)` when `params` are out of range.
pub fn keygen_from_seed(
    seed: &[u8; 32],
    params: &ThresholdParams,
) -> Result<([u8; PK_BYTES], Vec<ThresholdPrivateKey>), Error> {
    let n = params.n;
    let t = params.t;
    if !(t >= 2 && t <= n) || (n as usize) > MAX_PARTIES {
        return Err(Error::InvalidParameters);
    }

    // Expand seed via SHAKE-256
    let mut h = Shake256::default();
    h.update(seed);
    // Match Go (thmldsa44 has NIST = true): absorb (K, L) after the seed, as
    // FIPS 204 ML-DSA.KeyGen_internal does.
    h.update(&[K as u8, L as u8]);
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
            t,
            n,
            id: i,
            rho,
            key,
            tr: [0u8; TRBYTES],
            shares: BTreeMap::new(),
        });
        key.zeroize();
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
        share.s1h.reduce();
        share.s1h.ntt();
        share.s2h = share.s2.clone();
        share.s2h.reduce();
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
        sseed.zeroize();

        // Gosper's hack: advance to next subset bitmask with same popcount
        // (mask ≠ 0, so c ≠ 0)
        let c = mask & mask.wrapping_neg();
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
    // The totals are the full secret (s₁, s₂); wipe them.
    s1_total.zeroize();
    s1h_total.zeroize();
    s2_total.zeroize();
    s2h_total.zeroize();

    // 8. Compute tr = CRH(pk)
    let mut tr = [0u8; TRBYTES];
    let mut h2 = Shake256::default();
    h2.update(&pk_bytes);
    let mut tr_reader = h2.finalize_xof();
    tr_reader.read(&mut tr);

    // Set tr in all private keys
    for sk in &mut sks {
        sk.tr = tr;
    }

    Ok((pk_bytes, sks))
}

/// Sample a polynomial with coefficients uniformly in `[-η, η]`.
///
/// Uses SHAKE-256 rejection sampling (matching the Go reference
/// `PolyDeriveUniformLeqEta`).
pub(crate) fn sample_leq_eta(p: &mut Poly, seed: &[u8; 64], nonce: u16) {
    let mut h = Shake256::default();
    let mut iv = [0u8; 66];
    iv[..64].copy_from_slice(seed);
    iv[64] = nonce as u8;
    iv[65] = (nonce >> 8) as u8;
    h.update(&iv);

    iv.zeroize();

    let mut reader = h.finalize_xof();
    let mut buf = [0u8; 136]; // SHAKE-256 rate
    let mut i = 0;

    while i < N {
        reader.read(&mut buf);
        for byte in &buf {
            if i >= N {
                break;
            }
            let t1 = u32::from(*byte & 0x0F);
            let t2 = u32::from(*byte >> 4);

            // η=2: accept if < 15, reduce mod 5
            if t1 <= 14 {
                let reduced = t1 - ((205 * t1) >> 10) * 5;
                p.coeffs[i] = ETA as i32 - reduced as i32;
                i += 1;
            }
            if t2 <= 14 && i < N {
                let reduced = t2 - ((205 * t2) >> 10) * 5;
                p.coeffs[i] = ETA as i32 - reduced as i32;
                i += 1;
            }
        }
    }
    buf.zeroize();
}

const _: () = assert!(ETA == 2, "sample_leq_eta assumes η = 2");

/// Compute the public key bytes: `pk = (ρ ‖ t₁_packed)`.
///
/// Internally: t = `A·NTT⁻¹(s₁_hat)` + s₂, Power2Round(t) → (t₀, t₁)
pub(crate) fn compute_public_key(rho: &[u8; 32], s1h_total: &PolyVecL, s2_total: &PolyVecK) -> [u8; PK_BYTES] {
    pk_from_t(rho, &compute_t(rho, s1h_total, s2_total))
}

/// Compute t = A·s₁ + s₂ with coefficients in `[0, q)`.
///
/// `s1h` is s₁ in the NTT domain; `s2` is in the normal domain.
pub(crate) fn compute_t(rho: &[u8; 32], s1h: &PolyVecL, s2: &PolyVecK) -> PolyVecK {
    use dilithium::{
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
        s1h_d.vec[i].coeffs = s1h.polys[i].coeffs;
    }

    // t = A·s₁ (in NTT domain, then InvNTT + Montgomery)
    let mut t = DPolyVecK::default();
    matrix_pointwise_montgomery(mode, &mut t, &mat, &s1h_d);
    for i in 0..K {
        t.vec[i].invntt_tomont();
    }

    // t += s₂
    let mut s2_d = DPolyVecK::default();
    for i in 0..K {
        s2_d.vec[i].coeffs = s2.polys[i].coeffs;
    }
    // t_copy = A·s₁, which together with the public t would give s₂.
    let mut t_copy = t.clone();
    polyveck_add(mode, &mut t, &t_copy, &s2_d);
    // Threshold keygen accumulates subset secrets modulo q; after adding s2,
    // we must perform a full reduction before Power2Round.
    polyveck_reduce(mode, &mut t);
    polyveck_caddq(mode, &mut t);

    let mut out = PolyVecK::zero();
    for i in 0..K {
        out.polys[i].coeffs = t.vec[i].coeffs;
    }
    for p in s1h_d.vec.iter_mut().chain(s2_d.vec.iter_mut()).chain(t_copy.vec.iter_mut()).chain(t.vec.iter_mut()) {
        p.coeffs.zeroize();
    }
    out
}

/// Power2Round(t) → t₁ and pack pk = ρ ‖ t₁. `t` must be in `[0, q)`.
pub(crate) fn pk_from_t(rho: &[u8; 32], t: &PolyVecK) -> [u8; PK_BYTES] {
    use dilithium::{packing::pack_pk, poly::Poly as DPoly, polyvec::PolyVecK as DPolyVecK, ML_DSA_44};

    let mut t0 = DPolyVecK::default();
    let mut t1 = DPolyVecK::default();
    for i in 0..K {
        let mut ti = DPoly::default();
        ti.coeffs = t.polys[i].coeffs;
        DPoly::power2round(&mut t1.vec[i], &mut t0.vec[i], &ti);
        ti.coeffs.zeroize();
    }
    // t₀ is part of a standard ML-DSA secret key; do not leave it behind.
    for p in &mut t0.vec {
        p.coeffs.zeroize();
    }

    // Pack: pk = ρ ‖ t₁
    let mut pk = [0u8; PK_BYTES];
    pack_pk(ML_DSA_44, &mut pk, rho, &t1);
    pk
}

/// tr = SHAKE256(pk), 64 bytes.
pub(crate) fn tr_from_pk(pk: &[u8; PK_BYTES]) -> [u8; TRBYTES] {
    let mut tr = [0u8; TRBYTES];
    let mut h = Shake256::default();
    h.update(pk);
    h.finalize_xof().read(&mut tr);
    tr
}

/// Derive a subset share (s₁_I, s₂_I) ← χ_s from its 64-byte seed, with
/// NTT-domain caches. Matches the per-subset sampling in `keygen_from_seed`.
pub(crate) fn share_from_seed(sseed: &[u8; 64]) -> Share {
    let mut share = Share::zero();
    for j in 0..L {
        sample_leq_eta(&mut share.s1.polys[j], sseed, j as u16);
    }
    for j in 0..K {
        sample_leq_eta(&mut share.s2.polys[j], sseed, (j + L) as u16);
    }
    share.s1h = share.s1.clone();
    share.s1h.reduce();
    share.s1h.ntt();
    share.s2h = share.s2.clone();
    share.s2h.reduce();
    share.s2h.ntt();
    share
}

/// All (N−T+1)-subset bitmasks of [N] in ascending order.
pub(crate) fn subset_masks(n: u8, t: u8) -> Vec<u8> {
    let size = u32::from(n - t + 1);
    (0u16..(1u16 << n))
        .filter(|m| m.count_ones() == size)
        .map(|m| m as u8)
        .collect()
}

/// Legacy subset enumeration (preserved for backward compatibility).
///
/// Returns a list of subsets, where each subset is a sorted Vec of party indices.
#[must_use] 
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{get_threshold_params, num_subsets};

    #[test]
    fn test_keygen_produces_valid_keys() {
        let seed = [42u8; 32];
        let params = get_threshold_params(2, 3).unwrap();
        let (pk, sks) = keygen_from_seed(&seed, &params).unwrap();

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
        let (pk1, _) = keygen_from_seed(&seed, &params).unwrap();
        let (pk2, _) = keygen_from_seed(&seed, &params).unwrap();
        assert_eq!(pk1, pk2);
    }

    #[test]
    fn test_keygen_all_configs() {
        let configs: &[(u8, u8)] = &[(2, 2), (2, 3), (3, 3), (2, 4), (3, 4), (4, 4)];
        for &(t, n) in configs {
            let seed = [t + n; 32];
            let params = get_threshold_params(t, n).unwrap();
            let (_, sks) = keygen_from_seed(&seed, &params).unwrap();
            assert_eq!(sks.len(), n as usize, "Wrong key count for ({}, {})", t, n);
        }
    }

    #[test]
    fn test_subset_count_matches_binomial() {
        // For (T=2, N=3): C(3,2) = 3 subsets
        let seed = [1u8; 32];
        let params = get_threshold_params(2, 3).unwrap();
        let (_, sks) = keygen_from_seed(&seed, &params).unwrap();

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

    #[test]
    fn test_keygen_invalid_params_rejected() {
        let seed = [0u8; 32];
        let bad = ThresholdParams {
            t: 3,
            n: 2,
            r: 0.0,
            r1: 0.0,
            k_reps: 1,
            nu: 3.0,
        };
        assert!(matches!(
            keygen_from_seed(&seed, &bad),
            Err(Error::InvalidParameters)
        ));
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;
    extern crate alloc;
    use alloc::format;

    #[test]
    fn serialization_errors_and_debug() {
        let params = get_threshold_params(2, 3).unwrap();
        let (_, mut sks) = keygen_from_seed(&[8u8; 32], &params).unwrap();
        let s = format!("{:?}", sks[0]);
        assert!(s.contains("id: 0") && !s.contains("rho"));
        assert_eq!(sks[1].tr().len(), TRBYTES);

        // Wrong subset mask in the encoding.
        let mut bytes = sks[0].to_bytes().unwrap().to_vec();
        let first_mask = KEY_HEADER_BYTES;
        bytes[first_mask] ^= 0b100;
        assert_eq!(ThresholdPrivateKey::from_bytes(&bytes).err(), Some(Error::InvalidShare));

        // Coefficient outside [-10, 10] cannot be serialized.
        let mask = *sks[0].share_masks().next().unwrap();
        sks[0].shares.get_mut(&mask).unwrap().s1.polys[0].coeffs[0] = SHARE_BOUND + 1;
        assert_eq!(sks[0].to_bytes().err(), Some(Error::InvalidShare));
    }

    #[test]
    fn invalid_keygen_and_subset_inputs() {
        let mut p = get_threshold_params(2, 3).unwrap();
        p.t = 4;
        assert_eq!(keygen_from_seed(&[0u8; 32], &p).err(), Some(Error::InvalidParameters));
        assert!(enumerate_subsets(0, 0).is_empty());
        assert!(enumerate_subsets(3, 4).is_empty());
        assert_eq!(enumerate_subsets(3, 2).len(), 3);
    }
}

#[cfg(test)]
mod parse_tests {
    use super::*;

    #[test]
    fn keys_round_trip_and_bad_headers_are_rejected() {
        let params = get_threshold_params(3, 5).unwrap();
        let (_, sks) = keygen_from_seed(&[15u8; 32], &params).unwrap();
        for sk in &sks {
            let bytes = sk.to_bytes().unwrap();
            let back = ThresholdPrivateKey::from_bytes(&bytes).unwrap();
            assert_eq!(back.to_bytes().unwrap()[..], bytes[..]);
            assert_eq!((back.t, back.n, back.id, back.rho, back.key, back.tr), (sk.t, sk.n, sk.id, sk.rho, sk.key, sk.tr));
            for (mask, share) in &sk.shares {
                let other = &back.shares[mask];
                assert!(share.s1h.polys.iter().zip(&other.s1h.polys).all(|(a, b)| a.coeffs == b.coeffs));
                assert!(share.s2h.polys.iter().zip(&other.s2h.polys).all(|(a, b)| a.coeffs == b.coeffs));
            }
        }
        let good = sks[0].to_bytes().unwrap();
        let mut bad = good.to_vec();
        bad[0] = 1; // old version
        assert_eq!(ThresholdPrivateKey::from_bytes(&bad).err(), Some(Error::InvalidParameters));
        let mut bad = good.to_vec();
        bad[2] = 7; // unsupported N
        assert_eq!(ThresholdPrivateKey::from_bytes(&bad).err(), Some(Error::InvalidParameters));
        let mut bad = good.to_vec();
        bad[3] = 5; // id ≥ N
        assert_eq!(ThresholdPrivateKey::from_bytes(&bad).err(), Some(Error::InvalidParameters));
        let mut bad = good.to_vec();
        bad[KEY_HEADER_BYTES - 1] += 1; // share count
        assert_eq!(ThresholdPrivateKey::from_bytes(&bad).err(), Some(Error::InvalidParameters));
        assert_eq!(ThresholdPrivateKey::from_bytes(&good[..good.len() - 1]).err(), Some(Error::InvalidParameters));
        let mut bad = good.to_vec();
        let last = bad.len() - 1;
        bad[last] = 21; // coefficient 11
        assert_eq!(ThresholdPrivateKey::from_bytes(&bad).err(), Some(Error::InvalidShare));
    }
}
