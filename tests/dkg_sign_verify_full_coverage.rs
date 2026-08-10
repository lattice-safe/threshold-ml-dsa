//! Full Coverage Integration Tests: DKG + Threshold Sign + FIPS 204 Verify.
//!
//! Tests complete end-to-end flows for:
//! 1. DKG -> Threshold Sign -> FIPS 204 Verify across various (T, N) threshold setups
//! 2. A Posteriori Key Sharing -> Threshold Sign -> FIPS 204 Verify
//! 3. Interactive DKG protocol error handling and adversarial tampering
//! 4. Active signer set validation & edge case handling
//! 5. Signature verification robustness (tampered signatures, wrong messages, wrong keys)
//! 6. FIPS 204 compliance across ML-DSA-44, ML-DSA-65, and ML-DSA-87

use rand::rngs::StdRng;
use rand::SeedableRng;
use threshold_ml_dsa::{
    dkg,
    error::Error,
    params::{PK_BYTES, SK_BYTES},
    sdk::ThresholdMlDsa44Sdk,
    verify,
};

/// Helper: Execute full DKG protocol for N parties and return the SDK instance.
fn run_dkg_protocol(
    t: u8,
    n: u8,
    rng: &mut StdRng,
) -> Result<ThresholdMlDsa44Sdk, Error> {
    // Round 1: Generate contributions and commitments for all N parties
    let mut contribs = Vec::new();
    let mut commits = Vec::new();

    for i in 0..n {
        let (c, h) = dkg::dkg_round1(i, rng)?;
        contribs.push((i, c));
        commits.push((i, h));
    }

    // Create SDK from DKG
    ThresholdMlDsa44Sdk::from_dkg(&contribs, &commits, t, n, 10)
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
        let sdk = run_dkg_protocol(t, n, &mut rng)
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
    let sdk = run_dkg_protocol(3, 5, &mut rng).unwrap();

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
        (3, 4),
        (4, 5),
        (5, 6),
    ];

    for &(t, n) in configs {
        let seed = [t * 10 + n; 32];
        let (pk_bytes, sk_bytes) = verify::keygen(&seed);
        let mut sk_fixed = [0u8; SK_BYTES];
        sk_fixed.copy_from_slice(&sk_bytes);

        // Share the existing key
        let sdk = ThresholdMlDsa44Sdk::from_existing_key(&sk_fixed, &[123u8; 32], t, n, 10)
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
        let std_sig = verify::sign_standard(msg.as_bytes(), &sk_fixed, &[88u8; 32]);
        assert!(verify::verify(&std_sig, msg.as_bytes(), sdk.pk()));
    }
}

// ─── 3. DKG Error Handling and Adversarial Tests ─────────────

#[test]
fn test_dkg_tampered_contribution_fails_commitment_check() {
    let mut rng = StdRng::seed_from_u64(101);
    let n = 4u8;

    let mut contribs = Vec::new();
    let mut commits = Vec::new();

    for i in 0..n {
        let (c, h) = dkg::dkg_round1(i, &mut rng).unwrap();
        contribs.push((i, c));
        commits.push((i, h));
    }

    // Party 2 swaps contribution seed after commitment broadcast
    let mut tampered_seed = *contribs[2].1.seed();
    tampered_seed[0] ^= 0x01;
    contribs[2] = (2, dkg::DkgContribution::from_raw(2, tampered_seed));

    assert_eq!(
        dkg::dkg_round2(n, &contribs, &commits),
        Err(Error::InvalidShare)
    );
}

#[test]
fn test_dkg_duplicate_or_invalid_party_ids() {
    let mut rng = StdRng::seed_from_u64(202);
    let n = 3u8;

    let mut contribs = Vec::new();
    let mut commits = Vec::new();

    for i in 0..n {
        let (c, h) = dkg::dkg_round1(i, &mut rng).unwrap();
        contribs.push((i, c));
        commits.push((i, h));
    }

    // Case A: Duplicate party ID in contributions
    let mut dup_contribs = contribs.clone();
    dup_contribs[1].0 = 0; // party 1 claims to be party 0

    assert_eq!(
        dkg::dkg_round2(n, &dup_contribs, &commits),
        Err(Error::InvalidParameters)
    );

    // Case B: Party ID out of range (>= MAX_PARTIES)
    assert_eq!(
        dkg::dkg_round1(8, &mut rng).err(),
        Some(Error::InvalidParameters)
    );
}

// ─── 4. Threshold Signing Active Set Validation Tests ────────

#[test]
fn test_threshold_sign_active_set_validation() {
    let mut rng = StdRng::seed_from_u64(303);
    let sdk = run_dkg_protocol(3, 4, &mut rng).unwrap();
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
    let sdk = run_dkg_protocol(2, 3, &mut rng).unwrap();
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
    let sdk = run_dkg_protocol(2, 3, &mut rng).unwrap();
    let msg1 = b"Original message";
    let msg2 = b"Tampered message";

    let sig = sdk.threshold_sign(&[0, 1], msg1, &mut rng).unwrap();

    assert!(sdk.verify(msg1, &sig));
    assert!(!sdk.verify(msg2, &sig));
}

#[test]
fn test_verification_fails_on_wrong_public_key() {
    let mut rng = StdRng::seed_from_u64(606);
    let sdk1 = run_dkg_protocol(2, 3, &mut rng).unwrap();
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
    let sig44 = verify::sign_standard(msg, &sk44, &rng_seed);
    assert!(verify::verify(&sig44, msg, &pk44));

    // Level 3: ML-DSA-65
    let (pk65, sk65) = verify::keygen_65(&seed);
    let sig65 = verify::sign_standard_65(msg, &sk65, &rng_seed);
    assert!(verify::verify_65(&sig65, msg, &pk65));

    // Level 5: ML-DSA-87
    let (pk87, sk87) = verify::keygen_87(&seed);
    let sig87 = verify::sign_standard_87(msg, &sk87, &rng_seed);
    assert!(verify::verify_87(&sig87, msg, &pk87));
}
