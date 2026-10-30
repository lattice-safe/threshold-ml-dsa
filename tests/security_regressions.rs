//! Regression tests for the v0.5.0 security review findings.

use rand::rngs::StdRng;
use rand::SeedableRng;
use threshold_ml_dsa::coordinator;
use threshold_ml_dsa::error::Error;
use threshold_ml_dsa::params::*;
use threshold_ml_dsa::poly::PolyVecK;
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use threshold_ml_dsa::sign;

const MSG: &[u8] = b"security regression";
const SESSION: [u8; 32] = [7u8; 32];

/// Run Round 1 for both parties of a (2, 2) setup.
fn round1_both(sdk: &ThresholdMlDsa44Sdk, rng: &mut StdRng) -> (Vec<[u8; 32]>, Vec<sign::StRound1>) {
    let mut hashes = Vec::new();
    let mut states = Vec::new();
    for id in 0..2u8 {
        let sk = sdk.party_key(id as usize).unwrap();
        let (h, st) = sign::round1(sk, sdk.params(), 0b11, MSG, &SESSION, rng).unwrap();
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
            let sdk = ThresholdMlDsa44Sdk::from_seed(&[t * 10 + n; 32], t, n, 200).unwrap();
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
    for n in 7..=8u8 {
        for t in 2..=n {
            assert!(get_threshold_params(t, n).is_none(), "({t}, {n}) must be unsupported");
        }
    }
    assert_eq!(MAX_PARTIES, 6);
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
