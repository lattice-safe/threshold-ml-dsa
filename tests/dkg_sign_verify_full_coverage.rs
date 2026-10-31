//! Full Coverage Integration Tests: DKG + Threshold Sign + FIPS 204 Verify.
//!
//! Tests complete end-to-end flows for:
//! 1. DKG -> Threshold Sign -> FIPS 204 Verify across all (T, N) setups
//! 2. A Posteriori Key Sharing -> Threshold Sign -> FIPS 204 Verify
//! 3. DKG error handling, adversarial tampering and coalition privacy
//! 4. Active signer set validation & edge case handling
//! 5. Signature verification robustness (tampered signatures, wrong messages, wrong keys)
//! 6. FIPS 204 compliance across ML-DSA-44, ML-DSA-65, and ML-DSA-87

use rand::rngs::StdRng;
use rand::SeedableRng;
use threshold_ml_dsa::params::{PK_BYTES, SK_BYTES};
use threshold_ml_dsa::{
    dkg::{self, DkgConfig, DkgParty, DkgSeedMsg},
    error::Error,
    sdk::ThresholdMlDsa44Sdk,
    verify,
};

/// Helper: run the DKG for all N parties in-process and wrap the keys in an SDK.
fn run_dkg_protocol(t: u8, n: u8, rng: &mut StdRng) -> Result<ThresholdMlDsa44Sdk, Error> {
    let (pk, sks, _) = dkg::simulate(t, n, [0u8; 32], rng)?;
    ThresholdMlDsa44Sdk::from_party_keys(pk, sks, 100)
}

/// Run DKG rounds 0-1 for all parties.
fn dkg_until_round1(
    t: u8,
    n: u8,
    rng: &mut StdRng,
) -> (Vec<dkg::DkgAwaitingCoins>, Vec<dkg::DkgCoinMsg>, Vec<DkgSeedMsg>) {
    let mut parties = Vec::new();
    let mut r0 = Vec::new();
    for me in 0..n {
        let (p, m) = DkgParty::start(DkgConfig { t, n, me, context: [1u8; 32] }, rng).unwrap();
        parties.push(p);
        r0.push(m);
    }
    let mut st = Vec::new();
    let mut coins = Vec::new();
    let mut seeds = Vec::new();
    for p in parties {
        let (s, c, sm) = p.round1(&r0).unwrap();
        st.push(s);
        coins.push(c);
        seeds.extend(sm);
    }
    (st, coins, seeds)
}

/// Run DKG rounds 0-2 for all parties.
fn dkg_until_round2(t: u8, n: u8, rng: &mut StdRng) -> (Vec<dkg::DkgAwaitingCommits>, Vec<dkg::DkgCommitMsg>) {
    let (st, coins, seeds) = dkg_until_round1(t, n, rng);
    let mut st2 = Vec::new();
    let mut commits = Vec::new();
    for (me, s) in st.into_iter().enumerate() {
        let (s2, c) = s.round2(&coins, seeds_for(&seeds, me as u8), rng).unwrap();
        st2.push(s2);
        commits.push(c);
    }
    (st2, commits)
}

fn seeds_for(seeds: &[DkgSeedMsg], me: u8) -> Vec<DkgSeedMsg> {
    seeds
        .iter()
        .filter(|m| m.to() == me)
        .map(|m| DkgSeedMsg::from_bytes(&m.to_bytes()).unwrap())
        .collect()
}

// ─── 1. End-to-End DKG + Sign + Verify Tests ───────────────────────────────

#[test]
fn test_dkg_sign_verify_e2e_all_configs() {
    let mut rng = StdRng::seed_from_u64(42);

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
        let mut sdk = run_dkg_protocol(t, n, &mut rng)
            .unwrap_or_else(|e| panic!("DKG failed for ({}, {}): {:?}", t, n, e));

        assert_eq!(sdk.num_parties(), n as usize);
        assert_eq!(sdk.params().t, t);
        assert_eq!(sdk.params().n, n);

        // Sign message using canonical active set {0, ..., T-1}
        let active: Vec<u8> = (0..t).collect();
        let msg = format!("Message signed via DKG ({}, {})", t, n);

        let sig = sdk
            .threshold_sign(&active, msg.as_bytes(), &mut rng)
            .unwrap_or_else(|e| panic!("Signing failed for ({}, {}): {:?}", t, n, e));

        // Verify with SDK
        assert!(sdk.verify(msg.as_bytes(), &sig));

        // Verify with standalone FIPS 204 verifier
        assert!(verify::verify(&sig, msg.as_bytes(), sdk.pk()));
    }
}

#[test]
fn test_dkg_sign_verify_with_non_canonical_active_sets() {
    let mut rng = StdRng::seed_from_u64(999);

    // Test (T=3, N=5) with non-canonical active signer sets
    let mut sdk = run_dkg_protocol(3, 5, &mut rng).unwrap();

    let active_sets: &[&[u8]] = &[
        &[0, 1, 2],
        &[0, 2, 4],
        &[1, 3, 4],
        &[2, 3, 4],
        &[0, 3, 4],
    ];

    for &active in active_sets {
        let msg = format!("Non-canonical active set {:?}", active);
        let sig = sdk.threshold_sign(active, msg.as_bytes(), &mut rng).unwrap();

        assert!(sdk.verify(msg.as_bytes(), &sig));
        assert!(verify::verify(&sig, msg.as_bytes(), sdk.pk()));
    }
}

// ─── 2. A Posteriori Key Sharing + Sign + Verify Tests ───────────────────

#[test]
fn test_aposteriori_sign_verify_e2e_all_configs() {
    let mut rng = StdRng::seed_from_u64(777);

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
        let seed = [t * 10 + n; 32];
        let (pk_bytes, sk_bytes) = verify::keygen(&seed);
        let mut sk_fixed = [0u8; SK_BYTES];
        sk_fixed.copy_from_slice(&sk_bytes);

        // Share the existing key
        let mut sdk = ThresholdMlDsa44Sdk::from_existing_key(&sk_fixed, t, n, 100, &mut rng)
            .unwrap_or_else(|e| panic!("A posteriori sharing failed for ({}, {}): {:?}", t, n, e));

        // Public key must be identical to original FIPS 204 public key
        assert_eq!(sdk.pk(), &pk_bytes[..PK_BYTES]);

        // Sign with active set {0, ..., T-1}
        let active: Vec<u8> = (0..t).collect();
        let msg = format!("Message signed via aposteriori ({}, {})", t, n);

        let sig = sdk.threshold_sign(&active, msg.as_bytes(), &mut rng).unwrap();

        // Must verify against SDK, FIPS 204 verifier, and original public key
        assert!(sdk.verify(msg.as_bytes(), &sig));
        assert!(verify::verify(&sig, msg.as_bytes(), &pk_bytes));

        // Standard signature produced with original SK must also verify against same PK
        let std_sig = verify::sign_standard(msg.as_bytes(), &sk_fixed, &[88u8; 32]).unwrap();
        assert!(verify::verify(&std_sig, msg.as_bytes(), sdk.pk()));
    }
}

// ─── 3. DKG Error Handling and Adversarial Tests ─────────────

#[test]
fn test_dkg_keys_have_expected_structure() {
    let mut rng = StdRng::seed_from_u64(100);
    for n in 2..=6u8 {
        for t in 2..=n {
            let (pk, sks, _) = dkg::simulate(t, n, [0u8; 32], &mut rng).unwrap();
            for sk in &sks {
                let held: Vec<u8> = sk.share_masks().copied().collect();
                let expected: Vec<u8> = (0u16..1 << n)
                    .map(|m| m as u8)
                    .filter(|m| m.count_ones() == u32::from(n - t + 1) && m & (1 << sk.id()) != 0)
                    .collect();
                assert_eq!(held, expected, "({t}, {n}) party {}", sk.id());
                assert_eq!(sk.threshold(), (t, n));
            }
            // All parties agree on the key and their shares are consistent:
            // signing with any T of them works (checked in e2e tests).
            ThresholdMlDsa44Sdk::from_party_keys(pk, sks, 1).unwrap();
        }
    }
}

#[test]
fn test_dkg_coalition_never_receives_missing_subset() {
    let mut rng = StdRng::seed_from_u64(101);
    for (t, n) in [(2u8, 3u8), (3, 5), (4, 6), (3, 3)] {
        let (_, sks, transcript) = dkg::simulate(t, n, [0u8; 32], &mut rng).unwrap();
        // Every coalition C of size T-1 misses M = complement of C.
        for c in 0u16..1 << n {
            let c = c as u8;
            if c.count_ones() != u32::from(t - 1) {
                continue;
            }
            let full = ((1u16 << n) - 1) as u8;
            let missing = full & !c;
            if missing.count_ones() != u32::from(n - t + 1) {
                continue;
            }
            for sk in sks.iter().filter(|sk| c & (1 << sk.id()) != 0) {
                assert!(sk.share_masks().all(|&m| m != missing));
            }
            for m in &transcript.seeds {
                if c & (1 << m.to()) != 0 {
                    assert_ne!(m.mask(), missing, "seed of missing subset sent to coalition");
                }
            }
        }
    }
}

#[test]
fn test_dkg_tampered_reveal_is_detected() {
    let mut rng = StdRng::seed_from_u64(102);
    let (st2, commits) = dkg_until_round2(2, 3, &mut rng);
    let mut st3 = Vec::new();
    let mut reveals = Vec::new();
    for s in st2 {
        let (s3, r) = s.round3(&commits).unwrap();
        st3.push(s3);
        reveals.push(r);
    }
    reveals[1].t_parts[0].1[0] ^= 1;
    for s in st3 {
        assert_eq!(s.round4(&reveals).err(), Some(Error::DkgCommitmentMismatch { party: 1 }));
    }
}

#[test]
fn test_dkg_tampered_coin_is_detected() {
    let mut rng = StdRng::seed_from_u64(106);
    let (st, mut coins, seeds) = dkg_until_round1(2, 3, &mut rng);
    // A rushing party cannot change its coin share after seeing the others.
    coins[2].coin[0] ^= 1;
    for (me, s) in st.into_iter().enumerate() {
        let r = s.round2(&coins, seeds_for(&seeds, me as u8), &mut rng);
        let expected = if me == 2 { Error::DkgInconsistentTranscript } else { Error::DkgCommitmentMismatch { party: 2 } };
        assert_eq!(r.err(), Some(expected));
    }
}

#[test]
fn test_dkg_wrong_seed_is_detected_by_member() {
    let mut rng = StdRng::seed_from_u64(103);
    let (st, coins, seeds) = dkg_until_round1(2, 3, &mut rng);
    // Party 0 deals subset {0,1} (mask 0b011) to party 1; corrupt that key.
    let mut mine = seeds_for(&seeds, 1);
    let idx = mine.iter().position(|m| m.mask() == 0b011).unwrap();
    let mut bytes = mine[idx].to_bytes().to_vec();
    let last = bytes.len() - 1;
    bytes[last] ^= 1;
    mine[idx] = DkgSeedMsg::from_bytes(&bytes).unwrap();
    let s1 = st.into_iter().nth(1).unwrap();
    assert_eq!(
        s1.round2(&coins, mine, &mut rng).err(),
        Some(Error::DkgSeedMismatch { dealer: 0, mask: 0b011 })
    );
}

#[test]
fn test_dkg_dealer_commitment_checked_by_member() {
    let mut rng = StdRng::seed_from_u64(107);
    let (st2, mut commits) = dkg_until_round2(2, 3, &mut rng);
    // Party 0 commits to a t_I that does not match the key it sent:
    // party 1, a member of {0,1}, catches it before any reveal.
    commits[0].commits[0].1[0] ^= 1;
    let s1 = st2.into_iter().nth(1).unwrap();
    assert_eq!(s1.round3(&commits).err(), Some(Error::DkgCommitmentMismatch { party: 0 }));
}

#[test]
fn test_dkg_equivocation_is_detected_at_confirm() {
    let mut rng = StdRng::seed_from_u64(104);
    let (st2, commits) = dkg_until_round2(2, 3, &mut rng);
    // Party 0 equivocates: party 2 receives a different commitment for
    // subset {0,1} (mask 0b011), which party 2 is not a member of.
    let mut forked = commits.clone();
    forked[0].commits[0].1[0] ^= 1;
    let mut st3 = Vec::new();
    let mut reveals = Vec::new();
    for (me, s) in st2.into_iter().enumerate() {
        let view = if me == 2 { &forked } else { &commits };
        let (s3, r) = s.round3(view).unwrap();
        st3.push(s3);
        reveals.push(r);
    }
    // Party 2 detects the mismatch when the reveal arrives.
    let mut results = Vec::new();
    for s in st3 {
        results.push(s.round4(&reveals));
    }
    assert_eq!(results[2].as_ref().err(), Some(&Error::DkgCommitmentMismatch { party: 0 }));
    // Parties 0 and 1 proceed, but cannot finish without 2's confirmation.
    let mut confirms: Vec<_> = results.iter().filter_map(|r| r.as_ref().ok()).map(|(_, c)| c.clone()).collect();
    confirms.push(dkg::DkgConfirmMsg { from: 2, digest: [0u8; 32] });
    for r in results.into_iter().flatten() {
        assert_eq!(r.0.finalize(&confirms).err(), Some(Error::DkgInconsistentTranscript));
    }
}

#[test]
fn test_dkg_rejects_malformed_round_messages() {
    let mut rng = StdRng::seed_from_u64(105);
    // Invalid configuration
    assert_eq!(
        DkgParty::start(DkgConfig { t: 2, n: 3, me: 3, context: [0; 32] }, &mut rng).err(),
        Some(Error::InvalidParameters)
    );
    assert_eq!(
        DkgParty::start(DkgConfig { t: 2, n: 9, me: 0, context: [0; 32] }, &mut rng).err(),
        Some(Error::InvalidParameters)
    );

    // Duplicate sender / wrong key-commitment set / altered own message in round 0
    let mut parties = Vec::new();
    let mut r0 = Vec::new();
    for me in 0..3 {
        let (p, m) = DkgParty::start(DkgConfig { t: 2, n: 3, me, context: [0; 32] }, &mut rng).unwrap();
        parties.push(p);
        r0.push(m);
    }
    let mut dup = r0.clone();
    dup[1].from = 0;
    let mut extra = r0.clone();
    extra[2].key_commits.push((0b110, [0u8; 32]));
    let mut altered = r0.clone();
    altered[2].coin_commit[0] ^= 1;
    let mut parties = parties.into_iter();
    let p0 = parties.next().unwrap();
    let p1 = parties.next().unwrap();
    let p2 = parties.next().unwrap();
    assert_eq!(p0.round1(&dup).err(), Some(Error::DkgMalformedMessage));
    assert_eq!(p1.round1(&extra).err(), Some(Error::DkgMalformedMessage));
    assert_eq!(p2.round1(&altered).err(), Some(Error::DkgInconsistentTranscript));

    // Missing key / wrong session
    let (st, coins, seeds) = dkg_until_round1(2, 3, &mut rng);
    let mut st = st.into_iter();
    let _s0 = st.next().unwrap();
    let s1 = st.next().unwrap();
    let s2 = st.next().unwrap();
    // Party 2 deals nothing in (2, 3) and needs keys for {0,2} and {1,2}.
    let mut mine = seeds_for(&seeds, 2);
    assert_eq!(mine.len(), 2);
    mine.pop();
    assert_eq!(s2.round2(&coins, mine, &mut rng).err(), Some(Error::DkgMalformedMessage));
    let mut bad = coins.clone();
    bad[2].sid[0] ^= 1;
    assert_eq!(s1.round2(&bad, seeds_for(&seeds, 1), &mut rng).err(), Some(Error::DkgMalformedMessage));
    assert!(DkgSeedMsg::from_bytes(&[0u8; 10]).is_err());
}

// ─── 4. Threshold Signing Active Set Validation Tests ────────

#[test]
fn test_threshold_sign_active_set_validation() {
    let mut rng = StdRng::seed_from_u64(303);
    let mut sdk = run_dkg_protocol(3, 4, &mut rng).unwrap();
    let msg = b"Active set validation test";

    // Case 1: Wrong number of active signers (len != T)
    assert_eq!(
        sdk.threshold_sign(&[0, 1], msg, &mut rng),
        Err(Error::InvalidParameters)
    );
    assert_eq!(
        sdk.threshold_sign(&[0, 1, 2, 3], msg, &mut rng),
        Err(Error::InvalidParameters)
    );

    // Case 2: Duplicate signer ID
    assert_eq!(
        sdk.threshold_sign(&[0, 1, 1], msg, &mut rng),
        Err(Error::InvalidParameters)
    );

    // Case 3: Unsorted signer IDs
    assert_eq!(
        sdk.threshold_sign(&[0, 2, 1], msg, &mut rng),
        Err(Error::InvalidParameters)
    );

    // Case 4: Signer ID out of range (>= N)
    assert_eq!(
        sdk.threshold_sign(&[0, 1, 4], msg, &mut rng),
        Err(Error::InvalidParameters)
    );
}

// ─── 5. Verification Robustness & Cross-Verification Tests ───

#[test]
fn test_verification_fails_on_tampered_signature() {
    let mut rng = StdRng::seed_from_u64(404);
    let mut sdk = run_dkg_protocol(2, 3, &mut rng).unwrap();
    let msg = b"Tamper resistance test";

    let mut sig = sdk.threshold_sign(&[0, 1], msg, &mut rng).unwrap();

    // Verify valid signature
    assert!(sdk.verify(msg, &sig));

    // Flip bits in signature bytes
    sig[0] ^= 0x01;
    assert!(!sdk.verify(msg, &sig));
    assert!(!verify::verify(&sig, msg, sdk.pk()));

    // Flip bits in middle (z vector)
    sig[0] ^= 0x01; // restore
    sig[500] ^= 0x80;
    assert!(!sdk.verify(msg, &sig));
}

#[test]
fn test_verification_fails_on_wrong_message() {
    let mut rng = StdRng::seed_from_u64(505);
    let mut sdk = run_dkg_protocol(2, 3, &mut rng).unwrap();
    let msg1 = b"Original message";
    let msg2 = b"Tampered message";

    let sig = sdk.threshold_sign(&[0, 1], msg1, &mut rng).unwrap();

    assert!(sdk.verify(msg1, &sig));
    assert!(!sdk.verify(msg2, &sig));
}

#[test]
fn test_verification_fails_on_wrong_public_key() {
    let mut rng = StdRng::seed_from_u64(606);
    let mut sdk1 = run_dkg_protocol(2, 3, &mut rng).unwrap();
    let sdk2 = run_dkg_protocol(2, 3, &mut rng).unwrap();

    let msg = b"Cross-key verification test";
    let sig1 = sdk1.threshold_sign(&[0, 1], msg, &mut rng).unwrap();

    // Signature from sdk1 MUST NOT verify with sdk2's public key
    assert!(!sdk2.verify(msg, &sig1));
    assert!(!verify::verify(&sig1, msg, sdk2.pk()));
}

// ─── 6. Multi-Security-Level Verification Tests (ML-DSA-44/65/87) ───

#[test]
fn test_all_fips_204_security_levels_verify() {
    let seed = [77u8; 32];
    let msg = b"FIPS 204 multi-level verification check";
    let rng_seed = [88u8; 32];

    // Level 2: ML-DSA-44
    let (pk44, sk44) = verify::keygen(&seed);
    let sig44 = verify::sign_standard(msg, &sk44, &rng_seed).unwrap();
    assert!(verify::verify(&sig44, msg, &pk44));

    // Level 3: ML-DSA-65
    let (pk65, sk65) = verify::keygen_65(&seed);
    let sig65 = verify::sign_standard_65(msg, &sk65, &rng_seed).unwrap();
    assert!(verify::verify_65(&sig65, msg, &pk65));

    // Level 5: ML-DSA-87
    let (pk87, sk87) = verify::keygen_87(&seed);
    let sig87 = verify::sign_standard_87(msg, &sk87, &rng_seed).unwrap();
    assert!(verify::verify_87(&sig87, msg, &pk87));
}
