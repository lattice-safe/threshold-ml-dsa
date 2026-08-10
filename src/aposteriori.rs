//! A Posteriori Key Sharing for Threshold ML-DSA.
//!
//! Implements ePrint 2026/013 §4.2: converting an existing, standard ML-DSA-44
//! secret key into threshold Replicated Secret Sharing (RSS) key material
//! while preserving the original public key.
//!
//! ## Protocol Overview
//!
//! Given an existing secret key `(s₁, s₂)`:
//! 1. Enumerate all C(N, N-T+1) subsets `I ⊂ [N]` of size `N-T+1`.
//! 2. For all subsets except the last, sample fresh short shares `s_I ← U([-η, η])`.
//! 3. The last subset receives the remainder share:
//!    `s_last = s_orig - Σ_{I ≠ last} s_I`.
//! 4. Distribute shares to parties as in standard RSS keygen.
//!
//! The resulting threshold secret key material sums to the original secret key `(s₁, s₂)`,
//! so the corresponding threshold public key is identical to the original public key.

#[cfg(not(feature = "std"))]
use alloc::{collections::BTreeMap, vec::Vec};
#[cfg(feature = "std")]
use std::collections::BTreeMap;

use crate::error::Error;
use crate::params::{
    ThresholdParams, K, L, MAX_PARTIES, N, PK_BYTES, SK_BYTES, TRBYTES,
};
use crate::poly::{PolyVecK, PolyVecL};
use crate::rss::{compute_public_key, sample_leq_eta, Share, ThresholdPrivateKey};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// Decompose an existing ML-DSA-44 secret key into threshold RSS private keys.
///
/// # Arguments
/// * `sk_bytes` — Packed ML-DSA-44 secret key (`SK_BYTES` = 2592 bytes).
/// * `seed` — 32-byte seed for generating non-remainder subset shares.
/// * `params` — Threshold parameters (T, N).
///
/// # Returns
/// `Ok((public_key_bytes, Vec<ThresholdPrivateKey>))` — The original packed public key
/// and per-party threshold private keys.
///
/// # Errors
/// Returns `Error::InvalidParameters` if `sk_bytes` length is invalid or `params` are out of range.
pub fn share_existing_key(
    sk_bytes: &[u8],
    seed: &[u8; 32],
    params: &ThresholdParams,
) -> Result<([u8; PK_BYTES], Vec<ThresholdPrivateKey>), Error> {
    if sk_bytes.len() != SK_BYTES {
        return Err(Error::InvalidParameters);
    }
    let n = params.n;
    let t = params.t;
    if !(t >= 2 && t <= n) || (n as usize) > MAX_PARTIES {
        return Err(Error::InvalidParameters);
    }

    use dilithium::{
        packing::unpack_sk,
        polyvec::{PolyVecK as DPolyVecK, PolyVecL as DPolyVecL},
        ML_DSA_44,
    };
    let mode = ML_DSA_44;

    // 1. Unpack original ML-DSA-44 secret key
    let mut rho = [0u8; 32];
    let mut tr_orig = [0u8; TRBYTES];
    let mut key_orig = [0u8; 32];
    let mut t0_d = DPolyVecK::default();
    let mut s1_d = DPolyVecL::default();
    let mut s2_d = DPolyVecK::default();

    unpack_sk(
        mode,
        &mut rho,
        &mut tr_orig,
        &mut key_orig,
        &mut t0_d,
        &mut s1_d,
        &mut s2_d,
        sk_bytes,
    );

    // Convert dilithium-rs vectors to our internal PolyVecL and PolyVecK
    let mut s1_orig = PolyVecL::zero();
    for i in 0..L {
        for j in 0..N {
            s1_orig.polys[i].coeffs[j] = s1_d.vec[i].coeffs[j];
        }
    }

    let mut s2_orig = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            s2_orig.polys[i].coeffs[j] = s2_d.vec[i].coeffs[j];
        }
    }

    // Expand randomness seed for generating non-remainder subset shares
    let mut h = Shake256::default();
    h.update(b"th-ml-dsa-aposteriori-v1");
    h.update(seed);
    h.update(sk_bytes);
    let mut reader = h.finalize_xof();

    // 2. Initialize per-party private keys
    let mut sks: Vec<ThresholdPrivateKey> = Vec::with_capacity(n as usize);
    for i in 0..n {
        let mut pkey = [0u8; 32];
        reader.read(&mut pkey);
        sks.push(ThresholdPrivateKey {
            id: i,
            rho,
            key: pkey,
            tr: [0u8; TRBYTES],
            shares: BTreeMap::new(),
        });
    }

    // Collect all subset bitmasks via Gosper's hack
    let subset_size = n - t + 1;
    let mut mask: u16 = (1u16 << subset_size) - 1;
    let limit: u16 = 1u16 << n;
    let mut subset_masks = Vec::new();

    while mask < limit {
        subset_masks.push(mask as u8);
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

    if subset_masks.is_empty() {
        return Err(Error::InvalidParameters);
    }

    // Accumulator for non-remainder sampled shares
    let mut s1_acc = PolyVecL::zero();
    let mut s2_acc = PolyVecK::zero();

    let num_subsets = subset_masks.len();

    // 3. Process all subsets except the last one
    for &subset_mask in &subset_masks[..num_subsets - 1] {
        let mut sseed = [0u8; 64];
        reader.read(&mut sseed);

        let mut share = Share {
            s1: PolyVecL::zero(),
            s2: PolyVecK::zero(),
            s1h: PolyVecL::zero(),
            s2h: PolyVecK::zero(),
        };

        for j in 0..L {
            sample_leq_eta(&mut share.s1.polys[j], &sseed, j as u16);
        }
        for j in 0..K {
            sample_leq_eta(&mut share.s2.polys[j], &sseed, (j + L) as u16);
        }

        // Accumulate
        s1_acc.add_assign(&share.s1);
        s2_acc.add_assign(&share.s2);

        // Pre-compute NTT domain representation
        share.s1h = share.s1.clone();
        share.s1h.reduce();
        share.s1h.ntt();
        share.s2h = share.s2.clone();
        share.s2h.reduce();
        share.s2h.ntt();

        // Distribute to parties in this subset
        for i in 0..n {
            if subset_mask & (1 << i) != 0 {
                sks[i as usize].shares.insert(subset_mask, share.clone());
            }
        }
    }

    // 4. Last subset gets the remainder: s_last = s_orig - s_acc
    let last_mask = subset_masks[num_subsets - 1];
    let mut share_last = Share {
        s1: PolyVecL::zero(),
        s2: PolyVecK::zero(),
        s1h: PolyVecL::zero(),
        s2h: PolyVecK::zero(),
    };

    share_last.s1 = PolyVecL::zero();
    PolyVecL::sub(&mut share_last.s1, &s1_orig, &s1_acc);
    share_last.s2 = PolyVecK::zero();
    PolyVecK::sub(&mut share_last.s2, &s2_orig, &s2_acc);

    share_last.s1h = share_last.s1.clone();
    share_last.s1h.reduce();
    share_last.s1h.ntt();
    share_last.s2h = share_last.s2.clone();
    share_last.s2h.reduce();
    share_last.s2h.ntt();

    for i in 0..n {
        if last_mask & (1 << i) != 0 {
            sks[i as usize].shares.insert(last_mask, share_last.clone());
        }
    }

    // 5. Derive the public key from s1_orig and s2_orig
    let mut s1_orig_ntt = s1_orig.clone();
    s1_orig_ntt.reduce();
    s1_orig_ntt.ntt();

    let pk_bytes = compute_public_key(&rho, &s1_orig_ntt, &s2_orig);

    // 6. Derive tr = CRH(pk)
    let mut tr = [0u8; TRBYTES];
    let mut h_tr = Shake256::default();
    h_tr.update(&pk_bytes);
    let mut tr_reader = h_tr.finalize_xof();
    tr_reader.read(&mut tr);

    for sk in &mut sks {
        sk.tr = tr;
    }

    Ok((pk_bytes, sks))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::get_threshold_params;
    use crate::verify;

    #[test]
    fn test_aposteriori_sharing_preserves_pk_and_verifies() {
        let seed = [42u8; 32];
        let (pk_bytes, sk_bytes) = verify::keygen(&seed);
        let mut sk_fixed = [0u8; SK_BYTES];
        sk_fixed.copy_from_slice(&sk_bytes);

        let params = get_threshold_params(2, 3).unwrap();
        let (th_pk, sks) = share_existing_key(&sk_fixed, &[99u8; 32], &params).unwrap();

        // The derived public key from a posteriori sharing MUST match the original pk
        assert_eq!(th_pk, pk_bytes[..PK_BYTES]);
        assert_eq!(sks.len(), 3);

        // Sign a message using standard ML-DSA to ensure keypair is valid
        let msg = b"A posteriori key sharing verification test";
        let std_sig = verify::sign_standard(msg, &sk_fixed, &[1u8; 32]);
        assert!(verify::verify(&std_sig, msg, &th_pk));
    }

    #[test]
    fn test_sdk_aposteriori_threshold_sign() {
        use crate::sdk::ThresholdMlDsa44Sdk;
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let seed = [123u8; 32];
        let (pk_bytes, sk_bytes) = verify::keygen(&seed);

        let sdk = ThresholdMlDsa44Sdk::from_existing_key(&sk_bytes, &[42u8; 32], 2, 3, 10).unwrap();
        assert_eq!(sdk.pk(), &pk_bytes[..PK_BYTES]);

        let mut rng = StdRng::seed_from_u64(999);
        let msg = b"Threshold signing with aposteriori key sharing";
        let active = [0, 1]; // party 0 and 1 active out of 3

        let sig = sdk.threshold_sign(&active, msg, &mut rng).unwrap();
        assert!(sdk.verify(msg, &sig));
        assert!(verify::verify(&sig, msg, &pk_bytes));
    }
}
