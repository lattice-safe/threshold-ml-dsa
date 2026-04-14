//! Ported v0.2 coverage tests that exercise API-independent poly, verify, 
//! error, and params primitives.
//!
//! These tests were originally in full_coverage_tests.rs (v0.2) and
//! threshold_tests.rs (v0.2). They test invariants that are unchanged
//! across the v0.2 → v0.3 rewrite.

use threshold_ml_dsa::error::Error;
use threshold_ml_dsa::params::*;
use threshold_ml_dsa::poly::*;
use threshold_ml_dsa::rss;
use threshold_ml_dsa::verify;

// ═══════════════════════════════════════════════════════════════════════
//  POLY — NTT
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_ntt_identity_roundtrip() {
    let mut p = Poly::zero();
    p.coeffs[0] = 1;
    let mut ntt_p = p.clone();
    ntt_p.ntt();
    ntt_p.invntt_tomont();
    ntt_p.reduce();
    ntt_p.caddq();
    assert_ne!(ntt_p.coeffs[0], 0, "NTT roundtrip should produce non-zero");
}

#[test]
fn test_ntt_polynomial_product_linear() {
    let mut a = Poly::zero();
    let mut b = Poly::zero();
    a.coeffs[0] = 1;
    a.coeffs[1] = 2;
    b.coeffs[0] = 3;
    b.coeffs[1] = Q - 1;

    let mut a_ntt = a.clone();
    let mut b_ntt = b.clone();
    a_ntt.ntt();
    b_ntt.ntt();
    let mut c = Poly::zero();
    Poly::pointwise_montgomery(&mut c, &a_ntt, &b_ntt);
    c.invntt_tomont();
    c.reduce();
    c.caddq();

    assert_eq!(c.coeffs[0], 3);
    assert_eq!(c.coeffs[1], 5);
    assert_eq!(c.coeffs[2], Q - 2);
    for i in 3..N {
        assert_eq!(c.coeffs[i], 0, "coeff {i} should be 0");
    }
}

#[test]
fn test_ntt_product_constant() {
    let mut a = Poly::zero();
    let mut b = Poly::zero();
    a.coeffs[0] = 5;
    b.coeffs[0] = 7;

    let mut a_ntt = a.clone();
    let mut b_ntt = b.clone();
    a_ntt.ntt();
    b_ntt.ntt();
    let mut c = Poly::zero();
    Poly::pointwise_montgomery(&mut c, &a_ntt, &b_ntt);
    c.invntt_tomont();
    c.reduce();
    c.caddq();

    assert_eq!(c.coeffs[0], 35);
    for i in 1..N {
        assert_eq!(c.coeffs[i], 0, "coeff {i} should be 0");
    }
}

#[test]
fn test_ntt_product_with_zero() {
    let a = Poly::zero();
    let mut b = Poly::zero();
    b.coeffs[0] = 12345;
    b.coeffs[100] = 67890;

    let mut a_ntt = a.clone();
    let mut b_ntt = b.clone();
    a_ntt.ntt();
    b_ntt.ntt();
    let mut c = Poly::zero();
    Poly::pointwise_montgomery(&mut c, &a_ntt, &b_ntt);
    c.invntt_tomont();
    c.reduce();
    c.caddq();

    for i in 0..N {
        assert_eq!(c.coeffs[i], 0, "product with zero should be zero at {i}");
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  POLY — Arithmetic
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_poly_add_commutativity() {
    let mut a = Poly::zero();
    let mut b = Poly::zero();
    a.coeffs[0] = 100;
    a.coeffs[1] = 200;
    b.coeffs[0] = 300;
    b.coeffs[1] = 400;

    let mut c1 = Poly::zero();
    let mut c2 = Poly::zero();
    Poly::add(&mut c1, &a, &b);
    Poly::add(&mut c2, &b, &a);

    for i in 0..N {
        assert_eq!(c1.coeffs[i], c2.coeffs[i], "add not commutative at {i}");
    }
}

#[test]
fn test_poly_sub_inverse() {
    let mut a = Poly::zero();
    let mut b = Poly::zero();
    a.coeffs[0] = 500;
    a.coeffs[1] = 300;
    b.coeffs[0] = 200;
    b.coeffs[1] = 100;

    let mut sum = Poly::zero();
    Poly::add(&mut sum, &a, &b);
    let mut diff = Poly::zero();
    Poly::sub(&mut diff, &sum, &b);

    assert_eq!(diff.coeffs[0], a.coeffs[0]);
    assert_eq!(diff.coeffs[1], a.coeffs[1]);
}

#[test]
fn test_poly_reduce_to_range() {
    let mut p = Poly::zero();
    p.coeffs[0] = Q + 100;
    p.coeffs[1] = -50;
    p.coeffs[2] = 2 * Q + 7;
    p.reduce();
    assert_eq!(p.coeffs[0], 100);
    assert!(p.coeffs[1] >= 0 && p.coeffs[1] < Q);
    assert_eq!(p.coeffs[2], 7);
}

#[test]
fn test_poly_caddq_makes_positive() {
    let mut p = Poly::zero();
    p.coeffs[0] = -1;
    p.coeffs[1] = -Q;
    p.coeffs[2] = 0;
    p.coeffs[3] = 5;
    p.caddq();

    assert_eq!(p.coeffs[0], Q - 1);
    assert_eq!(p.coeffs[1], 0);
    assert_eq!(p.coeffs[2], 0);
    assert_eq!(p.coeffs[3], 5);
}

// ═══════════════════════════════════════════════════════════════════════
//  POLY — Norms
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_l2_norm_zero() {
    let p = Poly::zero();
    assert_eq!(p.l2_norm_squared(), 0);
}

#[test]
fn test_l2_norm_pythagorean() {
    let mut p = Poly::zero();
    p.coeffs[0] = 3;
    p.coeffs[1] = 4;
    assert_eq!(p.l2_norm_squared(), 25);
}

#[test]
fn test_l2_norm_negative_via_centering() {
    let mut p = Poly::zero();
    p.coeffs[0] = Q - 1;
    assert_eq!(p.l2_norm_squared(), 1);
}

#[test]
fn test_chknorm_within_bound() {
    let mut p = Poly::zero();
    p.coeffs[0] = 10;
    p.coeffs[1] = 20;
    assert!(!p.chknorm(30), "should not exceed bound 30");
}

#[test]
fn test_chknorm_exceeds_bound() {
    let mut p = Poly::zero();
    p.coeffs[0] = 50;
    assert!(p.chknorm(50), "50 >= 50 should exceed bound");
}

// ═══════════════════════════════════════════════════════════════════════
//  POLY — Sampling
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_uniform_sampling_deterministic() {
    let seed = [0u8; SEEDBYTES];
    let mut a1 = Poly::zero();
    let mut a2 = Poly::zero();
    Poly::uniform(&mut a1, &seed, 0);
    Poly::uniform(&mut a2, &seed, 0);
    assert_eq!(a1.coeffs, a2.coeffs, "same seed+nonce → same polynomial");
}

#[test]
fn test_uniform_sampling_different_nonces() {
    let seed = [0u8; SEEDBYTES];
    let mut a1 = Poly::zero();
    let mut a2 = Poly::zero();
    Poly::uniform(&mut a1, &seed, 0);
    Poly::uniform(&mut a2, &seed, 1);
    assert_ne!(a1.coeffs, a2.coeffs, "different nonces → different polynomials");
}

#[test]
fn test_uniform_sampling_in_range() {
    let seed = [42u8; SEEDBYTES];
    let mut a = Poly::zero();
    Poly::uniform(&mut a, &seed, 7);
    for (i, &c) in a.coeffs.iter().enumerate() {
        assert!((0..Q).contains(&c), "coefficient {i} = {c} out of [0,q)");
    }
}

#[test]
fn test_uniform_eta_bounded() {
    let seed = [0u8; CRHBYTES];
    let mut a = Poly::zero();
    Poly::uniform_eta(&mut a, &seed, 0);
    for (i, &c) in a.coeffs.iter().enumerate() {
        let centered = if c > Q / 2 { c - Q } else { c };
        assert!(centered.abs() <= ETA as i32, "eta sample at {i}: {centered} exceeds ±{ETA}");
    }
}

#[test]
fn test_uniform_gamma1_bounded() {
    let seed = [0u8; CRHBYTES];
    let mut a = Poly::zero();
    Poly::uniform_gamma1(&mut a, &seed, 0);
    for (i, &c) in a.coeffs.iter().enumerate() {
        let centered = if c > Q / 2 { c - Q } else { c };
        assert!(centered.abs() <= GAMMA1, "gamma1 sample at {i}: {centered} exceeds ±{GAMMA1}");
    }
}

#[test]
fn test_challenge_polynomial_weight() {
    let seed = [99u8; CTILDEBYTES];
    let mut c = Poly::zero();
    Poly::challenge(&mut c, &seed);

    let nonzero: usize = c.coeffs.iter().filter(|&&x| x != 0).count();
    assert_eq!(nonzero, TAU as usize, "challenge should have τ={TAU} nonzero coefficients");

    for &coeff in c.coeffs.iter() {
        assert!(coeff == 0 || coeff == 1 || coeff == Q - 1,
            "challenge coeff should be 0, 1, or -1 mod q; got {coeff}");
    }
}

#[test]
fn test_challenge_deterministic() {
    let seed = [7u8; CTILDEBYTES];
    let mut c1 = Poly::zero();
    let mut c2 = Poly::zero();
    Poly::challenge(&mut c1, &seed);
    Poly::challenge(&mut c2, &seed);
    assert_eq!(c1.coeffs, c2.coeffs, "same seed → same challenge");
}

// ═══════════════════════════════════════════════════════════════════════
//  POLY — Packing
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_pack_unpack_t1_roundtrip() {
    let mut p = Poly::zero();
    for i in 0..N {
        p.coeffs[i] = (i as i32 * 3 + 7) & 0x3FF;
    }
    let mut buf = vec![0u8; POLYT1_PACKEDBYTES];
    p.pack_t1(&mut buf);

    let mut q_poly = Poly::zero();
    Poly::unpack_t1(&mut q_poly, &buf);
    assert_eq!(p.coeffs, q_poly.coeffs, "t1 pack→unpack roundtrip failed");
}

#[test]
fn test_pack_unpack_z_roundtrip() {
    let mut p = Poly::zero();
    for i in 0..N {
        p.coeffs[i] = ((i as i32) * 137) % (2 * GAMMA1);
    }
    let mut buf = vec![0u8; POLYZ_PACKEDBYTES];
    p.pack_z(&mut buf);

    let mut q_poly = Poly::zero();
    Poly::unpack_z(&mut q_poly, &buf);

    for i in 0..N {
        let expected = p.coeffs[i] & 0x3FFFF;
        let got = q_poly.coeffs[i] & 0x3FFFF;
        assert_eq!(got, expected, "z pack→unpack mismatch at {i}");
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  POLY — Decompose / Power2Round / Hint
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_power2round_reconstruction() {
    for val in [0, 1, 100, 12345, Q / 2, Q - 1] {
        let (a1, a0) = power2round(val);
        assert_eq!((a1 << D) + a0, val, "Power2Round({val}) failed reconstruction");
    }
}

#[test]
fn test_decompose_values() {
    for val in [0, 1, 1000, Q / 3, Q / 2, Q - 1] {
        let (a1, a0) = decompose(val);
        let max_a1 = (Q - 1) / (2 * GAMMA2);
        assert!(a1 >= 0 && a1 <= max_a1, "a1={a1} out of range for val={val}");
        assert!(a0.abs() <= GAMMA2, "a0={a0} exceeds γ₂ for val={val}");
    }
}

#[test]
fn test_use_hint_no_hint() {
    let val = 100000;
    let (a1, _a0) = decompose(val);
    let a1_hint = use_hint(val, false);
    assert_eq!(a1, a1_hint, "use_hint with false should return HighBits");
}

#[test]
fn test_use_hint_with_hint() {
    let val = 100000;
    let a1_no = use_hint(val, false);
    let a1_yes = use_hint(val, true);
    assert_ne!(a1_no, a1_yes, "use_hint(true) should differ from use_hint(false)");
}

// ═══════════════════════════════════════════════════════════════════════
//  POLY — PolyVecL / PolyVecK
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_polyvecl_add_assign() {
    let mut a = PolyVecL::zero();
    let mut b = PolyVecL::zero();
    a.polys[0].coeffs[0] = 100;
    b.polys[0].coeffs[0] = 200;
    a.add_assign(&b);
    assert_eq!(a.polys[0].coeffs[0], 300);
}

#[test]
fn test_polyvecl_add_sub_inverse() {
    let mut a = PolyVecL::zero();
    let mut b = PolyVecL::zero();
    a.polys[0].coeffs[0] = 500;
    a.polys[1].coeffs[1] = 300;
    b.polys[0].coeffs[0] = 200;
    b.polys[1].coeffs[1] = 100;

    let mut sum = PolyVecL::zero();
    PolyVecL::add(&mut sum, &a, &b);
    let mut diff = PolyVecL::zero();
    PolyVecL::sub(&mut diff, &sum, &b);

    assert_eq!(diff.polys[0].coeffs[0], 500);
    assert_eq!(diff.polys[1].coeffs[1], 300);
}

#[test]
fn test_polyvecl_l2_norm() {
    let mut v = PolyVecL::zero();
    v.polys[0].coeffs[0] = 3;
    v.polys[0].coeffs[1] = 4;
    v.polys[1].coeffs[0] = 5;
    assert_eq!(v.l2_norm_squared(), 50);
}

#[test]
fn test_polyveck_add_assign() {
    let mut a = PolyVecK::zero();
    let mut b = PolyVecK::zero();
    a.polys[0].coeffs[0] = 100;
    b.polys[0].coeffs[0] = 200;
    a.add_assign(&b);
    assert_eq!(a.polys[0].coeffs[0], 300);
}

#[test]
fn test_polyveck_sub() {
    let mut a = PolyVecK::zero();
    let mut b = PolyVecK::zero();
    a.polys[0].coeffs[0] = 500;
    b.polys[0].coeffs[0] = 200;
    let mut c = PolyVecK::zero();
    PolyVecK::sub(&mut c, &a, &b);
    assert_eq!(c.polys[0].coeffs[0], 300);
}

// ═══════════════════════════════════════════════════════════════════════
//  POLY — MatrixA
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_matrix_expand_deterministic() {
    let rho = [42u8; SEEDBYTES];
    let a1 = MatrixA::expand(&rho);
    let a2 = MatrixA::expand(&rho);
    for i in 0..K {
        for j in 0..L {
            assert_eq!(a1.rows[i][j].coeffs, a2.rows[i][j].coeffs,
                "Matrix expansion should be deterministic");
        }
    }
}

#[test]
fn test_matrix_mul_vec_zero() {
    let rho = [0u8; SEEDBYTES];
    let a = MatrixA::expand(&rho);
    let v = PolyVecL::zero();
    let mut v_ntt = v.clone();
    v_ntt.ntt();
    let result = a.mul_vec(&v_ntt);
    for k in 0..K {
        for c in 0..N {
            assert_eq!(result.polys[k].coeffs[c], 0, "A*0 should be 0 at [{k}][{c}]");
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  RSS — enumerate_subsets (still public in v0.3)
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_enumerate_subsets_count() {
    let subs = rss::enumerate_subsets(3, 2);
    assert_eq!(subs.len(), 3); // C(3,2)=3

    let subs = rss::enumerate_subsets(4, 2);
    assert_eq!(subs.len(), 4); // C(4,3)=4

    let subs = rss::enumerate_subsets(5, 3);
    assert_eq!(subs.len(), 10); // C(5,3)=10

    let subs = rss::enumerate_subsets(6, 2);
    assert_eq!(subs.len(), 6); // C(6,5)=6
}

// ═══════════════════════════════════════════════════════════════════════
//  VERIFY MODULE — keygen, sign, verify
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_keygen_produces_valid_sizes() {
    let seed = [0u8; 32];
    let (pk, sk) = verify::keygen(&seed);
    assert_eq!(pk.len(), PK_BYTES, "pk size");
    assert!(!sk.is_empty(), "sk should be non-empty");
}

#[test]
fn test_keygen_deterministic() {
    let seed = [42u8; 32];
    let (pk1, sk1) = verify::keygen(&seed);
    let (pk2, sk2) = verify::keygen(&seed);
    assert_eq!(pk1, pk2, "same seed → same pk");
    assert_eq!(sk1, sk2, "same seed → same sk");
}

#[test]
fn test_keygen_different_seeds() {
    let (pk1, _) = verify::keygen(&[0u8; 32]);
    let (pk2, _) = verify::keygen(&[1u8; 32]);
    assert_ne!(pk1, pk2, "different seeds → different keys");
}

#[test]
fn test_sign_standard_produces_valid_size() {
    let seed = [42u8; 32];
    let (_, sk) = verify::keygen(&seed);
    let sig = verify::sign_standard(b"test", &sk, &[0u8; 32]);
    assert_eq!(sig.len(), SIG_BYTES, "signature size");
}

#[test]
fn test_standard_sign_verify_roundtrip() {
    let seed = [99u8; 32];
    let (pk, sk) = verify::keygen(&seed);
    let msg = b"roundtrip verification test";
    let sig = verify::sign_standard(msg, &sk, &[7u8; 32]);
    assert!(verify::verify(&sig, msg, &pk), "valid signature should verify");
}

#[test]
fn test_verify_wrong_message() {
    let seed = [50u8; 32];
    let (pk, sk) = verify::keygen(&seed);
    let sig = verify::sign_standard(b"correct", &sk, &[1u8; 32]);
    assert!(!verify::verify(&sig, b"wrong", &pk));
}

#[test]
fn test_verify_wrong_key() {
    let (pk1, sk1) = verify::keygen(&[10u8; 32]);
    let (pk2, _) = verify::keygen(&[20u8; 32]);
    let sig = verify::sign_standard(b"msg", &sk1, &[0u8; 32]);
    assert!(verify::verify(&sig, b"msg", &pk1), "correct key");
    assert!(!verify::verify(&sig, b"msg", &pk2), "wrong key");
}

#[test]
fn test_verify_corrupted_signature() {
    let seed = [60u8; 32];
    let (pk, sk) = verify::keygen(&seed);
    let mut sig = verify::sign_standard(b"msg", &sk, &[0u8; 32]);
    sig[10] ^= 0xFF;
    assert!(!verify::verify(&sig, b"msg", &pk), "corrupted sig should fail");
}

#[test]
fn test_verify_empty_message() {
    let seed = [70u8; 32];
    let (pk, sk) = verify::keygen(&seed);
    let msg: &[u8] = b"";
    let sig = verify::sign_standard(msg, &sk, &[0u8; 32]);
    assert!(verify::verify(&sig, msg, &pk), "empty message should work");
}

#[test]
fn test_sign_deterministic() {
    let seed = [80u8; 32];
    let (_, sk) = verify::keygen(&seed);
    let rng_seed = [5u8; 32];
    let sig1 = verify::sign_standard(b"test", &sk, &rng_seed);
    let sig2 = verify::sign_standard(b"test", &sk, &rng_seed);
    assert_eq!(sig1, sig2, "same inputs → same signature");
}

// ═══════════════════════════════════════════════════════════════════════
//  ERROR MODULE
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_error_display() {
    assert_eq!(format!("{}", Error::LocalRejectionAbort),
        "local hyperball rejection: ‖z_i‖₂² exceeded bound");
    assert_eq!(format!("{}", Error::InvalidShare), "invalid RSS share");
    assert_eq!(format!("{}", Error::InsufficientResponses),
        "insufficient valid responses from threshold parties");
    assert_eq!(format!("{}", Error::InvalidSignature),
        "signature verification failed");
    assert_eq!(format!("{}", Error::InvalidParameters),
        "invalid threshold parameters (N, T)");
    assert_eq!(format!("{}", Error::HintCheckFailed),
        "low-bits hint check failed");
}

#[test]
fn test_error_equality() {
    assert_eq!(Error::LocalRejectionAbort, Error::LocalRejectionAbort);
    assert_ne!(Error::LocalRejectionAbort, Error::InvalidShare);
    assert_ne!(Error::InvalidParameters, Error::InvalidSignature);
}

#[test]
fn test_error_clone() {
    let e = Error::InvalidShare;
    let cloned = e.clone();
    assert_eq!(e, cloned);
}

// ═══════════════════════════════════════════════════════════════════════
//  PARAMS — Constant Sanity Checks
// ═══════════════════════════════════════════════════════════════════════

#[test]
fn test_params_q_is_prime() {
    let q = Q as u64;
    for d in 2..2895 {
        assert!(d * (q / d) != q, "q={q} is divisible by {d}");
    }
}

#[test]
fn test_params_q_ntt_friendly() {
    assert_eq!(Q % 512, 1, "q should be 1 mod 512");
}

#[test]
fn test_params_consistency() {
    assert_eq!(BETA, (TAU * ETA) as i32);
    assert_eq!(GAMMA2, (Q - 1) / 88);
    assert_eq!(GAMMA1, 1 << 17);
    assert_eq!(PK_BYTES, SEEDBYTES + K * POLYT1_PACKEDBYTES);
    assert_eq!(SK_BYTES,
        3 * SEEDBYTES + TRBYTES
        + L * POLYETA_PACKEDBYTES + K * POLYETA_PACKEDBYTES
        + K * POLYT0_PACKEDBYTES);
    assert_eq!(SIG_BYTES, CTILDEBYTES + L * POLYZ_PACKEDBYTES + OMEGA + K);
}

#[test]
fn test_hyperball_bound() {
    let expected = (L as u64) * (N as u64) * (GAMMA1 as u64) * (GAMMA1 as u64);
    assert_eq!(HYPERBALL_BOUND_SQ, expected);
}
