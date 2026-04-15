use threshold_ml_dsa::poly::{decompose, Poly};

#[test]
fn test_invntt() {
    let mut p = Poly::zero();
    p.coeffs[0] = 1;
    p.ntt();
    p.invntt_tomont();
    p.reduce();
    p.caddq();
    let (a1, a0) = decompose(p.coeffs[0]);
    println!("a1={}, a0={}", a1, a0);
}
