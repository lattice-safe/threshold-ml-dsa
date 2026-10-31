//! Regression tests for the v0.5.0 security review findings.

use rand::rngs::StdRng;
use rand::SeedableRng;
use threshold_ml_dsa::coordinator;
use threshold_ml_dsa::error::Error;
use threshold_ml_dsa::params::*;
use threshold_ml_dsa::poly::PolyVecK;
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use threshold_ml_dsa::nonce::InMemoryNonceGuard;
use threshold_ml_dsa::sign;

const MSG: &[u8] = b"security regression";
const SESSION: [u8; 32] = [7u8; 32];

/// Run Round 1 for both parties of a (2, 2) setup.
fn round1_both(sdk: &ThresholdMlDsa44Sdk, rng: &mut StdRng) -> (Vec<[u8; 32]>, Vec<sign::StRound1>) {
    let mut hashes = Vec::new();
    let mut states = Vec::new();
    for id in 0..2u8 {
        let sk = sdk.party_key(id as usize).unwrap();
        let (h, st) = sign::round1(sk, sdk.params(), 0b11, MSG, &SESSION, rng, &mut InMemoryNonceGuard::new()).unwrap();
        hashes.push(h);
        states.push(st);
    }
    (hashes, states)
}

#[test]
fn all_supported_configs_sign_and_verify() {
    let mut rng = StdRng::seed_from_u64(31);
    for n in 2..=6u8 {
        for t in 2..=n {
            let mut sdk = ThresholdMlDsa44Sdk::from_seed(&[t * 10 + n; 32], t, n, 200).unwrap();
            let active: Vec<u8> = (n - t..n).collect();
            let sig = sdk
                .threshold_sign(&active, MSG, &mut rng)
                .unwrap_or_else(|e| panic!("({t}, {n}) failed: {e:?}"));
            assert!(sdk.verify(MSG, &sig), "({t}, {n}) signature does not verify");
        }
    }
}

#[test]
fn unsupported_party_counts_are_rejected() {
    // (5, 8), (6, 8), (7, 8) would need K in the thousands to tens of thousands.
    for (t, n) in [(5u8, 8u8), (6, 8), (7, 8), (2, 9), (9, 9)] {
        assert!(get_threshold_params(t, n).is_none(), "({t}, {n}) must be unsupported");
    }
    assert_eq!(MAX_PARTIES, 8);
}

#[test]
fn round2_rejects_transcript_without_own_commitment() {
    let sdk = ThresholdMlDsa44Sdk::from_seed(&[1u8; 32], 2, 2, 1).unwrap();
    let mut rng = StdRng::seed_from_u64(1);
    let (hashes, mut states) = round1_both(&sdk, &mut rng);
    let sk0 = sdk.party_key(0).unwrap();

    // Coordinator replaces party 0's slot with a commitment of its choosing.
    let forged = vec![[0xAAu8; 32], hashes[1]];
    let st0 = states.remove(0);
    let res = sign::round2(sk0, 0b11, MSG, &SESSION, &forged, &hashes[0], st0, sdk.params());
    assert_eq!(res.err(), Some(Error::InvalidShare));
}

#[test]
fn round2_rejects_wrong_transcript_length() {
    let sdk = ThresholdMlDsa44Sdk::from_seed(&[2u8; 32], 2, 2, 1).unwrap();
    let mut rng = StdRng::seed_from_u64(2);
    let (hashes, mut states) = round1_both(&sdk, &mut rng);
    let sk0 = sdk.party_key(0).unwrap();

    let extra = vec![hashes[0], hashes[1], [0u8; 32]];
    let st0 = states.remove(0);
    let res = sign::round2(sk0, 0b11, MSG, &SESSION, &extra, &hashes[0], st0, sdk.params());
    assert_eq!(res.err(), Some(Error::InvalidParameters));
}

#[test]
fn round3_rejects_non_canonical_aggregate_with_same_packing() {
    let sdk = ThresholdMlDsa44Sdk::from_seed(&[3u8; 32], 2, 2, 1).unwrap();
    let k_reps = sdk.params().k_reps as usize;
    let mut rng = StdRng::seed_from_u64(3);
    let (hashes, states) = round1_both(&sdk, &mut rng);

    let mut reveals = Vec::new();
    let mut st2s = Vec::new();
    for (idx, st1) in states.into_iter().enumerate() {
        let sk = sdk.party_key(idx).unwrap();
        let (rev, st2) =
            sign::round2(sk, 0b11, MSG, &SESSION, &hashes, &hashes[idx], st1, sdk.params()).unwrap();
        reveals.push(rev);
        st2s.push(st2);
    }

    let size = sign::pack_w_single_size();
    let parsed: Vec<Vec<PolyVecK>> = reveals
        .iter()
        .map(|r| (0..k_reps).map(|k| sign::unpack_w_single(&r[k * size..(k + 1) * size])).collect())
        .collect();
    let mut wfinals = coordinator::aggregate_commitments(&parsed, k_reps).unwrap();

    // Build an aggregate that packs to the same bytes but decomposes
    // differently: carry 2^23 into coefficient 10 and borrow 1 from 11.
    let p = &mut wfinals[0].polys[0].coeffs;
    let j = (11..256).find(|&j| p[j] & 1 == 1).unwrap();
    p[j - 1] += 1 << 23;
    p[j] -= 1;
    assert_eq!(
        sign::pack_w_single(&wfinals[0]),
        sign::pack_w_single(&coordinator::aggregate_commitments(&parsed, k_reps).unwrap()[0]),
        "test setup: aliased aggregate must pack identically"
    );

    let sk0 = sdk.party_key(0).unwrap();
    let verified = sign::verify_all_round2_reveals(
        sk0.tr(),
        &[0, 1],
        0b11,
        MSG,
        &SESSION,
        &reveals,
        &hashes,
        k_reps,
    )
    .unwrap();
    let st2 = st2s.remove(0);
    let res = sign::round3(sk0, &wfinals, st2, sdk.params(), verified);
    assert_eq!(res.err(), Some(Error::InvalidParameters));
}

#[test]
fn poly_samplers_match_fips204_reference() {
    use dilithium::{poly::Poly as DPoly, polyvec::matrix_expand, polyvec::PolyVecL as DPolyVecL, ML_DSA_44};
    use threshold_ml_dsa::poly::{MatrixA, Poly};

    let seed = [9u8; CRHBYTES];
    for nonce in 0..8u16 {
        let mut ours = Poly::zero();
        Poly::uniform_eta(&mut ours, &seed, nonce);
        let mut theirs = DPoly::default();
        DPoly::uniform_eta(ML_DSA_44, &mut theirs, &seed, nonce);
        assert_eq!(ours.coeffs, theirs.coeffs, "uniform_eta differs at nonce {nonce}");
    }

    let rho = [5u8; SEEDBYTES];
    let ours = MatrixA::expand(&rho);
    let mut theirs: [DPolyVecL; K] = core::array::from_fn(|_| DPolyVecL::default());
    matrix_expand(ML_DSA_44, &mut theirs, &rho);
    for i in 0..K {
        for j in 0..L {
            assert_eq!(ours.rows[i][j].coeffs, theirs[i].vec[j].coeffs, "A[{i}][{j}] differs");
        }
    }
}

#[test]
fn pack_z_accepts_non_centered_representatives() {
    use threshold_ml_dsa::poly::Poly;
    let mut centered = Poly::zero();
    let mut shifted = Poly::zero();
    for i in 0..N {
        let c = (i as i32 * 977) % GAMMA1 - GAMMA1 / 2;
        centered.coeffs[i] = c;
        shifted.coeffs[i] = if c < 0 { c + Q } else { c };
    }
    let mut a = [0u8; POLYZ_PACKEDBYTES];
    let mut b = [0u8; POLYZ_PACKEDBYTES];
    centered.pack_z(&mut a);
    shifted.pack_z(&mut b);
    assert_eq!(a, b);
}

#[test]
fn poly_debug_does_not_print_coefficients() {
    use threshold_ml_dsa::poly::Poly;
    let mut p = Poly::zero();
    p.coeffs[0] = 123_456;
    assert!(!format!("{p:?}").contains("123456"));
}

#[test]
fn tampered_params_are_rejected() {
    let sdk = ThresholdMlDsa44Sdk::from_seed(&[4u8; 32], 2, 2, 1).unwrap();
    let sk0 = sdk.party_key(0).unwrap();
    let mut rng = StdRng::seed_from_u64(4);
    for tamper in [
        |p: &mut ThresholdParams| p.r = f64::NAN,
        |p: &mut ThresholdParams| p.r = f64::INFINITY,
        |p: &mut ThresholdParams| p.r1 *= 0.5,
        |p: &mut ThresholdParams| p.k_reps += 1,
        |p: &mut ThresholdParams| p.n = 9,
    ] {
        let mut params = *sdk.params();
        tamper(&mut params);
        assert_eq!(
            sign::round1(sk0, &params, 0b11, MSG, &SESSION, &mut rng, &mut InMemoryNonceGuard::new()).err(),
            Some(Error::InvalidParameters)
        );
    }
    // Key generated for (2, 2) must not be used with (2, 3) parameters.
    let other = get_threshold_params(2, 3).unwrap();
    assert_eq!(
        sign::round1(sk0, &other, 0b11, MSG, &SESSION, &mut rng, &mut InMemoryNonceGuard::new()).err(),
        Some(Error::InvalidParameters)
    );
}

#[test]
fn verify_and_standard_sign_reject_malformed_lengths() {
    use threshold_ml_dsa::verify;
    let (pk, sk) = verify::keygen(&[6u8; 32]);
    let sig = verify::sign_standard(MSG, &sk, &[0u8; 32]).unwrap();
    assert!(verify::verify(&sig, MSG, &pk));
    for bad_pk in [&pk[..0], &pk[..PK_BYTES - 1], &pk[..10]] {
        assert!(!verify::verify(&sig, MSG, bad_pk));
    }
    assert!(!verify::verify_65(&sig, MSG, &pk[..10]));
    assert!(!verify::verify_87(&sig, MSG, &pk[..10]));
    assert!(verify::sign_standard(MSG, &sk[..10], &[0u8; 32]).is_err());
    assert!(verify::sign_standard_65(MSG, &sk[..10], &[0u8; 32]).is_err());
    assert!(verify::sign_standard_87(MSG, &sk[..10], &[0u8; 32]).is_err());
}

#[test]
fn coordinator_rejects_out_of_range_party_input() {
    use threshold_ml_dsa::coordinator::Signature;
    use threshold_ml_dsa::poly::PolyVecL;

    let mut z = PolyVecL::zero();
    z.polys[0].coeffs[0] = i32::MAX;
    let honest = PolyVecL::zero();
    assert_eq!(
        coordinator::aggregate_responses(&[vec![z], vec![honest]], 1).err(),
        Some(Error::InvalidShare)
    );
    let mut w = PolyVecK::zero();
    w.polys[0].coeffs[0] = -1;
    assert_eq!(
        coordinator::aggregate_commitments(&[vec![w], vec![PolyVecK::zero()]], 1).err(),
        Some(Error::InvalidShare)
    );

    let sig = Signature { c_tilde: vec![0u8; 3], z: PolyVecL::zero(), h: vec![0u8; OMEGA + K] };
    assert!(sig.to_bytes().is_err());
    let sig = Signature { c_tilde: vec![0u8; CTILDEBYTES], z: PolyVecL::zero(), h: vec![0u8; 500] };
    assert!(sig.to_bytes().is_err());
}

#[test]
fn private_key_serialization_roundtrip_and_validation() {
    use threshold_ml_dsa::rss::ThresholdPrivateKey;

    let sdk = ThresholdMlDsa44Sdk::from_seed(&[8u8; 32], 3, 5, 1).unwrap();
    for id in 0..5 {
        let sk = sdk.party_key(id).unwrap();
        let bytes = sk.to_bytes().unwrap();
        let parsed = ThresholdPrivateKey::from_bytes(&bytes).unwrap();
        assert_eq!(parsed.id(), sk.id());
        assert_eq!(parsed.threshold(), (3, 5));
        assert_eq!(*parsed.to_bytes().unwrap(), *bytes);
    }

    // A parsed key must be usable for signing.
    let bytes = sdk.party_key(0).unwrap().to_bytes().unwrap();
    let parsed = ThresholdPrivateKey::from_bytes(&bytes).unwrap();
    let mut rng = StdRng::seed_from_u64(8);
    assert!(sign::round1(&parsed, sdk.params(), 0b111, MSG, &SESSION, &mut rng, &mut InMemoryNonceGuard::new()).is_ok());

    let good = sdk.party_key(1).unwrap().to_bytes().unwrap();
    let mut bad = good.to_vec();
    bad[0] = 1; // old version
    assert!(ThresholdPrivateKey::from_bytes(&bad).is_err());
    let mut bad = good.to_vec();
    bad[2] = 9; // N = 9 unsupported
    assert!(ThresholdPrivateKey::from_bytes(&bad).is_err());
    let mut bad = good.to_vec();
    bad[3] = 5; // id ≥ N
    assert!(ThresholdPrivateKey::from_bytes(&bad).is_err());
    let mut bad = good.to_vec();
    let last = bad.len() - 1;
    bad[last] = 21; // coefficient out of [-10, 10]
    assert!(ThresholdPrivateKey::from_bytes(&bad).is_err());
    assert!(ThresholdPrivateKey::from_bytes(&good[..good.len() - 1]).is_err());
    assert!(ThresholdPrivateKey::from_bytes(&good[..10]).is_err());
}

/// Vectors from the Go reference `NewThresholdKeysFromSeed`
/// (github.com/Threshold-ML-DSA/Threshold-ML-DSA, thmldsa44) with seed 00..1f.
#[test]
fn keygen_from_seed_matches_go_reference() {
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    fn hex(b: &[u8]) -> String {
        b.iter().map(|x| format!("{x:02x}")).collect()
    }
    let seed: [u8; 32] = core::array::from_fn(|i| i as u8);
    for (t, n, pk_hash) in [
        (2u8, 2u8, "0002e0024bde3b28b17462c09fceaa56"),
        (2, 3, "b30d2fb8901d860b2a1c9c8613cb92fa"),
        (3, 5, "ee87427885fc43881a024c49ce056bfc"),
    ] {
        let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, t, n, 1).unwrap();
        let pk = sdk.pk();
        assert_eq!(hex(&pk[..32]), "d7b2b47254aae0db45e7930d4a98d2c97d8f1397d1789dafa17024b316e9bec9");
        let mut h = sha3::Shake256::default();
        h.update(pk);
        let mut out = [0u8; 16];
        h.finalize_xof().read(&mut out);
        assert_eq!(hex(&out), pk_hash, "({t}, {n}) pk differs from Go");
        assert_eq!(hex(&sdk.party_key(0).unwrap().tr()[..16]), pk_hash);
    }
}

/// An RNG that always returns the same bytes (simulates a repeated RNG).
struct StuckRng;
impl rand::RngCore for StuckRng {
    fn next_u32(&mut self) -> u32 {
        7
    }
    fn next_u64(&mut self) -> u64 {
        7
    }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        dest.fill(7);
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        dest.fill(7);
        Ok(())
    }
}
impl rand::CryptoRng for StuckRng {}

#[test]
fn nonce_guard_refuses_replayed_session() {
    let sdk = ThresholdMlDsa44Sdk::from_seed(&[10u8; 32], 2, 2, 1).unwrap();
    let sk0 = sdk.party_key(0).unwrap();
    let mut guard = InMemoryNonceGuard::new();
    let mut rng = StdRng::seed_from_u64(10);
    sign::round1(sk0, sdk.params(), 0b11, MSG, &SESSION, &mut rng, &mut guard).unwrap();
    assert_eq!(
        sign::round1(sk0, sdk.params(), 0b11, MSG, &SESSION, &mut rng, &mut guard).err(),
        Some(Error::SessionReplay)
    );
    // A different key may use the same session ID.
    let sk1 = sdk.party_key(1).unwrap();
    assert!(sign::round1(sk1, sdk.params(), 0b11, MSG, &SESSION, &mut rng, &mut guard).is_ok());
}

#[test]
fn nonce_guard_prevents_nonce_reuse_with_stuck_rng() {
    use threshold_ml_dsa::nonce::NonceGuard;
    let sdk = ThresholdMlDsa44Sdk::from_seed(&[11u8; 32], 2, 2, 1).unwrap();
    let sk0 = sdk.party_key(0).unwrap();
    // Capacity 1: session A is evicted after B, so A is accepted again, but
    // the counter must still give a different nonce (and commitment).
    let mut guard = InMemoryNonceGuard::with_capacity(1);
    let a = [1u8; 32];
    let b = [2u8; 32];
    let (h1, _) = sign::round1(sk0, sdk.params(), 0b11, MSG, &a, &mut StuckRng, &mut guard).unwrap();
    sign::round1(sk0, sdk.params(), 0b11, MSG, &b, &mut StuckRng, &mut guard).unwrap();
    let (h2, _) = sign::round1(sk0, sdk.params(), 0b11, MSG, &a, &mut StuckRng, &mut guard).unwrap();
    assert_ne!(h1, h2, "same session + stuck RNG must not reproduce the nonce");

    let id = sk0.guard_key_id();
    let mut g = InMemoryNonceGuard::new().with_start_counter(id, 41);
    assert_eq!(g.reserve(&id, &[3u8; 32]).unwrap(), 41);
    assert_eq!(g.reserve(&id, &[4u8; 32]).unwrap(), 42);
    assert_eq!(g.counter(&id), 43);
}

#[test]
fn aposteriori_rejects_malformed_secret_keys() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use threshold_ml_dsa::aposteriori::share_existing_key;
    use threshold_ml_dsa::verify;
    let mut rng = StdRng::seed_from_u64(12);
    let params = get_threshold_params(2, 3).unwrap();
    let mut share = |sk: &[u8]| share_existing_key(sk, &params, &mut rng);
    let (_, sk) = verify::keygen(&[12u8; 32]);
    assert!(share(&sk).is_ok());

    // tr no longer matches the key (bytes 64..128 hold tr).
    let mut bad = sk.to_vec();
    bad[64] ^= 1;
    assert!(share(&bad).is_err());

    // s₁ coefficient outside [-η, η] (also breaks tr; the η check alone is
    // tested in aposteriori::tests).
    let mut bad = sk.to_vec();
    bad[128] = 0xFF; // first s₁ byte (sk = ρ ‖ K ‖ tr ‖ s₁ ‖ s₂ ‖ t₀)
    assert!(share(&bad).is_err());
}

#[test]
fn aposteriori_keys_serialize_and_sign_like_fresh_keys() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use threshold_ml_dsa::rss::ThresholdPrivateKey;
    use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
    use threshold_ml_dsa::verify;

    let mut rng = StdRng::seed_from_u64(2027);
    let (pk, sk) = verify::keygen(&[21u8; 32]);
    let params = get_threshold_params(3, 4).unwrap();
    let (tpk, sks) = threshold_ml_dsa::aposteriori::share_existing_key(&sk, &params, &mut rng).unwrap();
    assert_eq!(tpk[..], pk[..]);

    // Gaussian shares exceed η; the key format must round-trip them.
    let parsed: Vec<ThresholdPrivateKey> =
        sks.iter().map(|k| ThresholdPrivateKey::from_bytes(&k.to_bytes().unwrap()).unwrap()).collect();

    // Signing must succeed at about the same per-attempt rate as with fresh
    // keys: the paper keeps all signing parameters unchanged (σ = √2 has the
    // same variance as U([-2, 2])).
    let mut shared = ThresholdMlDsa44Sdk::from_party_keys(tpk, parsed, 1).unwrap();
    let mut fresh = ThresholdMlDsa44Sdk::from_seed(&[22u8; 32], 3, 4, 1).unwrap();
    let (mut ok_shared, mut ok_fresh) = (0i32, 0i32);
    let attempts = 200u32;
    for i in 0..attempts {
        let msg = i.to_le_bytes();
        if let Ok(sig) = shared.threshold_sign(&[0, 2, 3], &msg, &mut rng) {
            assert!(verify::verify(&sig, &msg, &pk));
            ok_shared += 1;
        }
        if fresh.threshold_sign(&[0, 2, 3], &msg, &mut rng).is_ok() {
            ok_fresh += 1;
        }
    }
    assert!(ok_shared > 0 && ok_fresh > 0);
    let diff = (ok_shared - ok_fresh).abs();
    assert!(diff <= 40, "a posteriori {ok_shared}/{attempts} vs fresh {ok_fresh}/{attempts}");
}
