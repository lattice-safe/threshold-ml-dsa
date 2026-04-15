//! v0.3 integration tests for the paper-faithful ePrint 2026/013 implementation.

use threshold_ml_dsa::fvec::{sample_hyperball, FVec};
use threshold_ml_dsa::params::*;
use threshold_ml_dsa::partition;
use threshold_ml_dsa::rss;
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;

#[test]
fn test_keygen_from_seed_deterministic() {
    let seed = [42u8; 32];
    let params = get_threshold_params(2, 3).unwrap();
    let (pk1, _) = rss::keygen_from_seed(&seed, &params).unwrap();
    let (pk2, _) = rss::keygen_from_seed(&seed, &params).unwrap();
    assert_eq!(pk1, pk2, "Keygen should be deterministic");
}

#[test]
fn test_keygen_all_configs_produce_valid_keys() {
    let configs: &[(u8, u8)] = &[
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
        // (2, 6) through (6, 6) are slow — test the first few
    ];

    for &(t, n) in configs {
        let seed = [t.wrapping_add(n); 32];
        let params = get_threshold_params(t, n).unwrap();
        let (pk, sks) = rss::keygen_from_seed(&seed, &params).unwrap();

        // N private keys
        assert_eq!(sks.len(), n as usize, "Wrong key count for ({}, {})", t, n);

        // Each party holds correct subsets
        for sk in &sks {
            for (&mask, _) in &sk.shares {
                assert!(
                    mask & (1 << sk.id) != 0,
                    "Party {} has share for subset it doesn't belong to",
                    sk.id
                );
            }
        }

        // Public key has correct size
        assert_eq!(pk.len(), PK_BYTES);

        // All parties share the same tr
        let tr0 = sks[0].tr;
        for sk in &sks {
            assert_eq!(sk.tr, tr0, "All parties must share the same tr");
        }
    }
}

#[test]
fn test_partition_all_configs() {
    let configs: &[(u8, u8)] = &[
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

    for &(t, n) in configs {
        let active: Vec<u8> = (0..t).collect();
        let p = partition::rss_recover(&active, n, t).unwrap();
        assert_eq!(p.len(), t as usize, "partition size for ({}, {})", t, n);

        // All subsets should be covered exactly once
        let mut all: Vec<u8> = p.iter().flat_map(|v| v.iter().copied()).collect();
        all.sort();
        let before = all.len();
        all.dedup();
        assert_eq!(
            before,
            all.len(),
            "Duplicate in partition for ({}, {})",
            t,
            n
        );

        // Each party should only hold subsets it belongs to
        for (i, shares) in p.iter().enumerate() {
            let party = active[i];
            for &mask in shares {
                assert!(
                    mask & (1 << party) != 0,
                    "Party {} in partition has subset 0b{:b} it doesn't belong to",
                    party,
                    mask
                );
            }
        }
    }
}

#[test]
fn test_hyperball_sampling_norm_bound() {
    let rhop = [0u8; 64];
    for nonce in 0..10u16 {
        let mut fv = FVec::zero();
        sample_hyperball(&mut fv, 252778.0, 3.0, &rhop, nonce);

        // The ν-scaled norm should be ≤ radius
        assert!(
            !fv.excess(252778.0, 3.0),
            "Hyperball sample at nonce {} exceeds radius",
            nonce
        );
    }
}

#[test]
fn test_sdk_creation_and_verify_roundtrip() {
    let seed = [7u8; 32];
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 3, 10).unwrap();

    // Public key should verify with the standard ML-DSA-44 verifier
    // (We can't sign yet — protocol is not fully wired — but we can
    // verify that the SDK creates valid key material)
    assert_eq!(sdk.sks.len(), 3);
    assert_eq!(sdk.params.t, 2);
    assert_eq!(sdk.params.n, 3);
    assert_eq!(sdk.params.k_reps, 3); // From Figure 8: (2,3) → K=3
}

#[test]
fn test_threshold_params_match_paper() {
    // Verify specific values from Figure 8
    let p = get_threshold_params(2, 2).unwrap();
    assert_eq!(p.k_reps, 2);
    assert!((p.r - 252778.0).abs() < 0.01);
    assert!((p.r1 - 252833.0).abs() < 0.01);
    assert!((p.nu - 3.0).abs() < 0.01);

    let p = get_threshold_params(5, 6).unwrap();
    assert_eq!(p.k_reps, 100);
    assert!((p.r - 250590.0).abs() < 0.01);
}

#[test]
fn test_end_to_end_threshold_sign_2_2() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let seed = [42u8; 32];
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 2, 64).unwrap();

    let mut rng = StdRng::seed_from_u64(12345);
    let msg = b"Hello, threshold ML-DSA!";
    let active = [0u8, 1];

    let sig = sdk
        .threshold_sign(&active, msg, &mut rng)
        .expect("2-of-2 threshold_sign should produce a valid signature");
    assert!(
        sdk.verify(msg, &sig),
        "Threshold signature failed FIPS 204 verification"
    );
    assert!(!sdk.verify(b"wrong message", &sig));
}

#[test]
fn test_end_to_end_threshold_sign_2_3() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let seed = [7u8; 32];
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 3, 64).unwrap();

    let mut rng = StdRng::seed_from_u64(99999);
    let msg = b"threshold 2-of-3 test";
    let active = [0u8, 1];

    let sig = sdk
        .threshold_sign(&active, msg, &mut rng)
        .expect("2-of-3 threshold_sign should produce a valid signature");
    assert!(
        sdk.verify(msg, &sig),
        "2-of-3 threshold signature failed verification"
    );
}

#[test]
fn test_end_to_end_threshold_sign_4_6() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let seed = [46u8; 32];
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 4, 6, 32).unwrap();

    let mut rng = StdRng::seed_from_u64(464646);
    let msg = b"threshold 4-of-6 test";
    // Pick 4 arbitrary parties out of the 6
    let active = [1u8, 2, 4, 5];

    let sig = sdk
        .threshold_sign(&active, msg, &mut rng)
        .expect("4-of-6 threshold_sign should produce a valid signature");
    assert!(
        sdk.verify(msg, &sig),
        "4-of-6 threshold signature failed verification"
    );
}

#[test]
fn test_end_to_end_threshold_sign_5_6() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let seed = [56u8; 32];
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 5, 6, 32).unwrap();

    let mut rng = StdRng::seed_from_u64(565656);
    let msg = b"threshold 5-of-6 test with extremely large sub-shares";
    let active = [0u8, 1, 3, 4, 5];

    let sig = sdk
        .threshold_sign(&active, msg, &mut rng)
        .expect("5-of-6 threshold_sign should produce a valid signature");
    assert!(
        sdk.verify(msg, &sig),
        "5-of-6 threshold signature failed verification"
    );
}
