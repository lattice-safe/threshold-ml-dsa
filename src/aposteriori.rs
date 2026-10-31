//! A posteriori key sharing: split an existing ML-DSA-44 secret key into
//! threshold Replicated Secret Sharing (RSS) keys with the same public key.
//!
//! Follows §3.4 of the paper (USENIX Security '26 version of ePrint
//! 2026/013): "we sample a short secret sharing (as in the key generation of
//! 3.1) with the difference that the shares will be conditioned on their sum
//! being s."
//!
//! ## Algorithm
//!
//! For `S = C(N, N−T+1)` subsets, every coefficient position of `(s₁, s₂)`
//! is shared independently: `(x_1, …, x_S)` is drawn from the discrete
//! Gaussian of parameter σ = √2 (density ∝ `exp(−x²/4)`, the same variance
//! as `U([−2, 2])` in fresh keygen) restricted to the coset
//! `{x ∈ Z^S : Σ x_i = s}` of the root lattice `{Σ x_i = 0}`, with every
//! `|x_i| ≤ 10` (tail mass below 2⁻⁴⁴ per coordinate).
//!
//! The paper samples this coset Gaussian with a GPV sampler analysed under
//! Rényi divergence. This module samples it **exactly** by rejection:
//! propose `x_1..x_{S−1}` i.i.d. from the truncated `D_{Z,σ}`, set
//! `x_S = s − Σ x_i`, and accept with probability `exp(−x_S²/4)·[|x_S| ≤ 10]`.
//! The accepted vector has density ∝ `Π exp(−x_i²/4)` on the coset, which is
//! symmetric in the S shares: no share is a special "remainder".
//!
//! Sampling is constant-time in the secret: the Gaussian sampler is a
//! full-table CDT scan, the acceptance probability is looked up with a full
//! table scan, and each coefficient runs a fixed batch of
//! [`TRIALS`] proposals and selects the first accepted one without
//! branching. A batch fails for some coefficient with probability below
//! 2⁻³³ per key; the batch is then repeated with fresh randomness.
//!
//! ## Security
//!
//! A short sharing of a fixed secret cannot be perfectly hiding: T−1
//! parties see every share but one, which is a hint on `s`. The paper
//! captures this as hint-MLWE. Its formal result (Theorem 3.3) needs a
//! larger share width (at most 7 bits of security lost, up to 10× the
//! communication) and assumes Gaussian secrets, and the parameters for it
//! are only in the full version. With σ = √2, as here, signing parameters
//! and communication are unchanged, but the bound is heuristic: the paper
//! estimates the loss at **at most 12 bits** with the lattice estimator.
//! Prefer fresh keys ([`crate::rss::keygen_from_seed`] or [`crate::dkg`])
//! when the public key does not have to be kept.
//!
//! The function runs on the key owner's machine, which sees the whole
//! secret. Deliver each party's key over an authenticated, encrypted
//! channel and delete the original secret key afterwards.

#[cfg(not(feature = "std"))]
use alloc::{collections::BTreeMap, vec::Vec};
#[cfg(feature = "std")]
use std::collections::BTreeMap;

use crate::error::Error;
use crate::params::{ThresholdParams, ETA, K, L, MAX_PARTIES, N, PK_BYTES, SK_BYTES, TRBYTES};
use crate::poly::{PolyVecK, PolyVecL};
use crate::rss::{compute_public_key, subset_masks, Share, ThresholdPrivateKey};
use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use subtle::ConstantTimeEq;
use zeroize::{Zeroize, Zeroizing};

/// Largest share coefficient magnitude produced by a posteriori sharing.
pub const SHARE_BOUND: i32 = 10;

/// Proposals per coefficient in one constant-time batch.
pub const TRIALS: usize = 128;

/// Maximum number of subsets (C(6, 3) = 20).
const MAX_SUBSETS: usize = 20;

/// `CDT[i] = ⌊2⁶⁴ · P(X ≤ −10 + i)⌋` for `X ~ D_{Z,√2}` truncated to `[−10, 10]`.
const CDT: [u64; 20] = [
    0x0000_0000_044e_bd3c,
    0x0000_0001_f632_0f77,
    0x0000_008a_4ed6_54ec,
    0x0000_172f_e5f7_ef81,
    0x0002_5f41_b96c_55ec,
    0x0026_0fa4_9b64_4581,
    0x0178_ab43_285d_f6be,
    0x0915_390c_d667_794c,
    0x23a6_5816_4624_d361,
    0x5be4_515f_7ac7_82c8,
    0xa41b_aea0_8538_7d37,
    0xdc59_a7e9_b9db_2c9e,
    0xf6ea_c6f3_2998_86b3,
    0xfe87_54bc_d7a2_0941,
    0xffd9_f05b_649b_ba7e,
    0xfffd_a0be_4693_aa13,
    0xffff_e8d0_1a08_107e,
    0xffff_ff75_b129_ab13,
    0xffff_fffe_09cd_f088,
    0xffff_ffff_fbb1_42c3,
];

/// `ACCEPT[k] = ⌊2⁶⁴ · exp(−k²/4)⌋` for `k ∈ [0, 10]`.
const ACCEPT: [u128; 11] = [
    0x1_0000_0000_0000_0000,
    0xc75f_7cf5_6410_5743,
    0x5e2d_58d8_b3bc_df1a,
    0x1afb_718e_8457_f74e,
    0x04b0_556e_084f_3d1d,
    0x007e_83a5_61ab_af1a,
    0x0008_1679_1293_2a2c,
    0x0000_5047_ec36_81dc,
    0x0000_01e3_55bb_aee8,
    0x0000_0006_e4f7_8789,
    0x0000_0000_0f45_1bd2,
];

/// 1 if `a < b`, else 0, without branching (`a, b < 2¹²⁷`).
#[inline]
fn ct_lt(a: u128, b: u128) -> u128 {
    a.wrapping_sub(b) >> 127
}

/// All-ones if `a == b`, else 0, without branching.
#[inline]
fn ct_eq_mask(a: i32, b: i32) -> i32 {
    let d = (a ^ b) as u32;
    // (d | -d) has its top bit set iff d != 0.
    (((d | d.wrapping_neg()) >> 31) as i32) - 1
}

/// Sample from `D_{Z,√2}` truncated to `[−10, 10]` in constant time.
fn sample_gaussian(xof: &mut impl XofReader) -> i32 {
    let mut b = [0u8; 8];
    xof.read(&mut b);
    let u = u128::from(u64::from_le_bytes(b));
    let mut x = -SHARE_BOUND;
    for &c in &CDT {
        x += (1 - ct_lt(u, u128::from(c))) as i32;
    }
    b.zeroize();
    x
}

/// `ACCEPT[|x|]`, or 0 if `|x| > 10`, in constant time.
fn accept_threshold(x: i32) -> u128 {
    let sign = x >> 31;
    let a = (x ^ sign) - sign;
    let mut out = 0u128;
    for (k, &v) in ACCEPT.iter().enumerate() {
        out |= v & (ct_eq_mask(a, k as i32) as i64 as i128 as u128);
    }
    out
}

/// Share one coefficient `s` into `out[..num]` with `Σ out = s`, from the
/// σ = √2 Gaussian restricted to that coset. Returns `false` if no proposal
/// in the batch was accepted (`out` is then all zero).
fn share_coefficient(s: i32, out: &mut [i32], xof: &mut impl XofReader) -> bool {
    let num = out.len();
    let mut cand = [0i32; MAX_SUBSETS];
    out.iter_mut().for_each(|o| *o = 0);
    // All-ones once a proposal has been accepted.
    let mut done = 0i32;
    for _ in 0..TRIALS {
        let mut sum = 0i32;
        for c in cand.iter_mut().take(num - 1) {
            *c = sample_gaussian(xof);
            sum += *c;
        }
        cand[num - 1] = s - sum;
        let mut b = [0u8; 8];
        xof.read(&mut b);
        let u = u128::from(u64::from_le_bytes(b));
        let accept = -(ct_lt(u, accept_threshold(cand[num - 1])) as i32);
        let take = accept & !done;
        for (o, &c) in out.iter_mut().zip(cand.iter()) {
            *o |= c & take;
        }
        done |= accept;
    }
    cand.zeroize();
    done != 0
}

/// Share an existing ML-DSA-44 secret key among N parties with threshold T.
///
/// Returns the original packed public key and one private key per party.
/// See the module documentation for the algorithm and its security loss.
///
/// # Errors
/// - `Error::InvalidParameters` if `sk_bytes` has the wrong length, `params`
///   are not the canonical parameters for (T, N), the key's `s₁`/`s₂` are
///   not η-short, or its stored `tr` does not match its public key.
pub fn share_existing_key<R: RngCore + CryptoRng>(
    sk_bytes: &[u8],
    params: &ThresholdParams,
    rng: &mut R,
) -> Result<([u8; PK_BYTES], Vec<ThresholdPrivateKey>), Error> {
    if sk_bytes.len() != SK_BYTES {
        return Err(Error::InvalidParameters);
    }
    params.validate()?;
    let (t, n) = (params.t, params.n);
    if !(t >= 2 && t <= n) || (n as usize) > MAX_PARTIES {
        return Err(Error::InvalidParameters);
    }

    use dilithium::{
        packing::unpack_sk,
        polyvec::{PolyVecK as DPolyVecK, PolyVecL as DPolyVecL},
        ML_DSA_44,
    };

    let mut rho = [0u8; 32];
    let mut tr_orig = [0u8; TRBYTES];
    let mut key_orig = Zeroizing::new([0u8; 32]);
    let mut t0_d = DPolyVecK::default();
    let mut s1_d = DPolyVecL::default();
    let mut s2_d = DPolyVecK::default();
    unpack_sk(ML_DSA_44, &mut rho, &mut tr_orig, &mut key_orig, &mut t0_d, &mut s1_d, &mut s2_d, sk_bytes);

    let mut s1_orig = PolyVecL::zero();
    for i in 0..L {
        s1_orig.polys[i].coeffs.copy_from_slice(&s1_d.vec[i].coeffs[..N]);
    }
    let mut s2_orig = PolyVecK::zero();
    for i in 0..K {
        s2_orig.polys[i].coeffs.copy_from_slice(&s2_d.vec[i].coeffs[..N]);
    }
    for p in s1_d.vec.iter_mut().chain(s2_d.vec.iter_mut()).chain(t0_d.vec.iter_mut()) {
        p.coeffs.zeroize();
    }

    let result = share_vectors(&s1_orig, &s2_orig, &rho, &tr_orig, sk_bytes, params, rng);
    s1_orig.zeroize();
    s2_orig.zeroize();
    result
}

fn share_vectors<R: RngCore + CryptoRng>(
    s1_orig: &PolyVecL,
    s2_orig: &PolyVecK,
    rho: &[u8; 32],
    tr_orig: &[u8; TRBYTES],
    sk_bytes: &[u8],
    params: &ThresholdParams,
    rng: &mut R,
) -> Result<([u8; PK_BYTES], Vec<ThresholdPrivateKey>), Error> {
    let (t, n) = (params.t, params.n);

    // Reject malformed keys: s₁, s₂ must be η-short.
    let eta = ETA as i32;
    if s1_orig.polys.iter().chain(s2_orig.polys.iter()).any(|p| p.coeffs.iter().any(|&c| !(-eta..=eta).contains(&c))) {
        return Err(Error::InvalidParameters);
    }

    // The stored tr must match the public key the shares will produce.
    let mut s1h = s1_orig.clone();
    s1h.reduce();
    s1h.ntt();
    let pk = compute_public_key(rho, &s1h, s2_orig);
    s1h.zeroize();
    let mut tr = [0u8; TRBYTES];
    let mut h = Shake256::default();
    h.update(&pk);
    h.finalize_xof().read(&mut tr);
    if !bool::from(tr.ct_eq(tr_orig)) {
        return Err(Error::InvalidParameters);
    }

    // Hedged randomness: fresh RNG output plus the key itself.
    let mut fresh = Zeroizing::new([0u8; 64]);
    rng.fill_bytes(&mut fresh[..]);
    let mut h = Shake256::default();
    h.update(b"th-ml-dsa-aposteriori-gauss-v1");
    h.update(&[t, n]);
    h.update(&fresh[..]);
    h.update(sk_bytes);
    let mut xof = h.finalize_xof();

    let masks = subset_masks(n, t);
    let num = masks.len();
    if !(2..=MAX_SUBSETS).contains(&num) {
        return Err(Error::InvalidParameters);
    }

    let mut shares: Vec<Share> = (0..num)
        .map(|_| Share { s1: PolyVecL::zero(), s2: PolyVecK::zero(), s1h: PolyVecL::zero(), s2h: PolyVecK::zero() })
        .collect();
    let mut out = [0i32; MAX_SUBSETS];
    for poly in 0..L + K {
        for j in 0..N {
            let s = if poly < L { s1_orig.polys[poly].coeffs[j] } else { s2_orig.polys[poly - L].coeffs[j] };
            while !share_coefficient(s, &mut out[..num], &mut xof) {}
            for (share, &x) in shares.iter_mut().zip(out.iter()) {
                if poly < L {
                    share.s1.polys[poly].coeffs[j] = x;
                } else {
                    share.s2.polys[poly - L].coeffs[j] = x;
                }
            }
        }
    }
    out.zeroize();

    let mut sks: Vec<ThresholdPrivateKey> = (0..n)
        .map(|id| {
            let mut key = [0u8; 32];
            xof.read(&mut key);
            ThresholdPrivateKey { t, n, id, rho: *rho, key, tr, shares: BTreeMap::new() }
        })
        .collect();
    for (mask, mut share) in masks.into_iter().zip(shares) {
        share.s1h = share.s1.clone();
        share.s1h.reduce();
        share.s1h.ntt();
        share.s2h = share.s2.clone();
        share.s2h.reduce();
        share.s2h.ntt();
        for sk in sks.iter_mut().filter(|sk| mask & (1 << sk.id) != 0) {
            sk.shares.insert(mask, share.clone());
        }
    }

    Ok((pk, sks))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::get_threshold_params;
    use crate::verify;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use sha3::Shake128;

    fn xof() -> impl XofReader {
        let mut h = Shake128::default();
        h.update(b"aposteriori-test");
        h.finalize_xof()
    }

    #[test]
    fn gaussian_sampler_matches_sigma_sqrt2() {
        let mut x = xof();
        let trials = 200_000;
        let (mut sum, mut sq) = (0i64, 0i64);
        for _ in 0..trials {
            let v = i64::from(sample_gaussian(&mut x));
            sum += v;
            sq += v * v;
        }
        let mean = sum as f64 / trials as f64;
        let var = sq as f64 / trials as f64 - mean * mean;
        assert!(mean.abs() < 0.02, "mean {mean}");
        // Var(D_{Z,√2}) ≈ 2.0
        assert!((var - 2.0).abs() < 0.04, "var {var}");
    }

    #[test]
    fn accept_threshold_table() {
        assert_eq!(accept_threshold(0), 1u128 << 64);
        assert_eq!(accept_threshold(-3), ACCEPT[3]);
        assert_eq!(accept_threshold(10), ACCEPT[10]);
        assert_eq!(accept_threshold(11), 0);
        assert_eq!(accept_threshold(-40), 0);
    }

    #[test]
    fn coefficient_shares_sum_and_are_symmetric() {
        // S = 20 (T = 4, N = 6): first and last share must have the same
        // distribution, i.e. there is no "remainder" share.
        let mut x = xof();
        let num = 20;
        let trials = 20_000;
        let mut out = [0i32; MAX_SUBSETS];
        let (mut m_first, mut m_last, mut v_first, mut v_last) = (0f64, 0f64, 0f64, 0f64);
        for i in 0..trials {
            let s = i % 5 - 2;
            while !share_coefficient(s, &mut out[..num], &mut x) {}
            assert_eq!(out[..num].iter().sum::<i32>(), s);
            assert!(out[..num].iter().all(|c| c.abs() <= SHARE_BOUND));
            m_first += f64::from(out[0]);
            m_last += f64::from(out[num - 1]);
            v_first += f64::from(out[0] * out[0]);
            v_last += f64::from(out[num - 1] * out[num - 1]);
        }
        let n = f64::from(trials);
        let (vf, vl) = (v_first / n, v_last / n);
        assert!((m_first / n).abs() < 0.05 && (m_last / n).abs() < 0.05);
        // Conditioned variance ≈ σ²(1 − 1/S) + E[s²]/S² ≈ 1.9 for both.
        assert!((vf - 1.9).abs() < 0.1, "first share variance {vf}");
        assert!((vl - 1.9).abs() < 0.1, "last share variance {vl}");
    }

    #[test]
    fn sharing_preserves_pk_and_sums_to_secret() {
        let (pk_bytes, sk_bytes) = verify::keygen(&[42u8; 32]);
        let params = get_threshold_params(3, 5).unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        let (pk, sks) = share_existing_key(&sk_bytes, &params, &mut rng).unwrap();
        assert_eq!(pk[..], pk_bytes[..PK_BYTES]);
        assert_eq!(sks.len(), 5);

        // Σ_I s_I equals the original s₁ (take each subset once).
        let mut total = PolyVecL::zero();
        let mut seen = 0u32;
        for sk in &sks {
            for (&mask, share) in &sk.shares {
                if seen & (1 << mask) == 0 {
                    seen |= 1 << mask;
                    total.add_assign(&share.s1);
                    assert!(share.s1.polys.iter().all(|p| p.coeffs.iter().all(|c| c.abs() <= SHARE_BOUND)));
                }
            }
        }
        use dilithium::{packing::unpack_sk, polyvec::{PolyVecK as DK, PolyVecL as DL}, ML_DSA_44};
        let (mut rho, mut tr, mut key) = ([0u8; 32], [0u8; TRBYTES], [0u8; 32]);
        let (mut t0, mut s1, mut s2) = (DK::default(), DL::default(), DK::default());
        unpack_sk(ML_DSA_44, &mut rho, &mut tr, &mut key, &mut t0, &mut s1, &mut s2, &sk_bytes);
        for i in 0..L {
            assert_eq!(total.polys[i].coeffs[..], s1.vec[i].coeffs[..N]);
        }
    }

    #[test]
    fn sdk_threshold_sign_with_shared_key() {
        use crate::sdk::ThresholdMlDsa44Sdk;
        let (pk_bytes, sk_bytes) = verify::keygen(&[123u8; 32]);
        let mut rng = StdRng::seed_from_u64(999);
        let mut sdk = ThresholdMlDsa44Sdk::from_existing_key(&sk_bytes, 2, 3, 10, &mut rng).unwrap();
        assert_eq!(sdk.pk(), &pk_bytes[..PK_BYTES]);
        let msg = b"Threshold signing with a posteriori key sharing";
        let sig = sdk.threshold_sign(&[0, 1], msg, &mut rng).unwrap();
        assert!(verify::verify(&sig, msg, &pk_bytes));
    }
}
