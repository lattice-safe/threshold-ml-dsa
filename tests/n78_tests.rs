//! N = 7, 8 support: parameters, partitions, keygen, DKG, a posteriori
//! sharing and signing.

use rand::rngs::StdRng;
use rand::SeedableRng;
use threshold_ml_dsa::params::get_threshold_params;
use threshold_ml_dsa::partition::rss_recover;
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use threshold_ml_dsa::{dkg, verify};

const SUPPORTED: [(u8, u8); 10] = [(2, 7), (3, 7), (4, 7), (5, 7), (6, 7), (7, 7), (2, 8), (3, 8), (4, 8), (8, 8)];

fn binom(n: u32, k: u32) -> u32 {
    (0..k).fold(1, |acc, i| acc * (n - i) / (i + 1))
}

#[test]
fn parameters_exist_and_validate() {
    for (t, n) in SUPPORTED {
        let p = get_threshold_params(t, n).unwrap_or_else(|| panic!("({t}, {n}) missing"));
        assert!(p.validate().is_ok());
        assert!(p.r1 > p.r && p.nu == 3.0, "({t}, {n})");
    }
}

#[test]
fn partitions_cover_every_subset_once_for_every_active_set() {
    for (t, n) in SUPPORTED {
        let size = u32::from(n - t + 1);
        let all: Vec<u8> = (0u16..1 << n).filter(|m| m.count_ones() == size).map(|m| m as u8).collect();
        let bound = binom(u32::from(n), u32::from(t) - 1).div_ceil(u32::from(t)) as usize;
        for act in (0u16..1 << n).filter(|a| a.count_ones() == u32::from(t)) {
            let active: Vec<u8> = (0..n).filter(|i| act & (1 << i) != 0).collect();
            let p = rss_recover(&active, n, t).unwrap();
            let mut seen: Vec<u8> = Vec::new();
            for (i, shares) in p.iter().enumerate() {
                assert!(shares.len() <= bound, "({t}, {n}) load {} > {bound}", shares.len());
                for &m in shares {
                    assert!(m & (1 << active[i]) != 0, "({t}, {n}) party not in subset");
                    seen.push(m);
                }
            }
            seen.sort_unstable();
            assert_eq!(seen, all, "({t}, {n}) active {active:?}");
        }
    }
}

fn sign_and_verify(t: u8, n: u8, seed: u64) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut sdk = ThresholdMlDsa44Sdk::from_seed(&[t * 16 + n; 32], t, n, 40).unwrap();
    let active: Vec<u8> = (n - t..n).collect(); // a non-canonical set
    let msg = format!("N = {n}, T = {t}");
    let sig = sdk.threshold_sign(&active, msg.as_bytes(), &mut rng).unwrap();
    assert!(verify::verify(&sig, msg.as_bytes(), sdk.pk()));
}

#[test]
fn small_k_configs_sign_and_verify() {
    for (i, (t, n)) in [(2u8, 7u8), (3, 7), (7, 7), (2, 8), (3, 8)].into_iter().enumerate() {
        sign_and_verify(t, n, 700 + i as u64);
    }
}

/// K = 227–1195 per attempt: slow in debug builds.
/// Run with `cargo test --release --test n78_tests -- --ignored`.
#[test]
#[ignore]
fn large_k_configs_sign_and_verify() {
    for (i, (t, n)) in [(8u8, 8u8), (4, 7), (5, 7), (6, 7), (4, 8)].into_iter().enumerate() {
        sign_and_verify(t, n, 800 + i as u64);
    }
}

#[test]
fn dkg_with_eight_parties() {
    let mut rng = StdRng::seed_from_u64(78);
    let (pk, sks, _) = dkg::simulate(3, 8, [0u8; 32], &mut rng).unwrap();
    assert_eq!(sks.len(), 8);
    let mut sdk = ThresholdMlDsa44Sdk::from_party_keys(pk, sks, 40).unwrap();
    let sig = sdk.threshold_sign(&[1, 4, 7], b"dkg n=8", &mut rng).unwrap();
    assert!(verify::verify(&sig, b"dkg n=8", &pk));
}

#[test]
fn aposteriori_sharing_with_more_than_twenty_subsets() {
    // (3, 8): C(8, 6) = 28 subsets, so the larger proposal batch is used.
    let (pk, sk) = verify::keygen(&[87u8; 32]);
    let mut rng = StdRng::seed_from_u64(87);
    let mut sdk = ThresholdMlDsa44Sdk::from_existing_key(&sk, 3, 8, 40, &mut rng).unwrap();
    assert_eq!(sdk.pk(), &pk[..]);
    let sig = sdk.threshold_sign(&[0, 3, 6], b"a posteriori n=8", &mut rng).unwrap();
    assert!(verify::verify(&sig, b"a posteriori n=8", &pk));
}
