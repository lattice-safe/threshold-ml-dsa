use threshold_ml_dsa::fvec::FVec;
use threshold_ml_dsa::poly::{PolyVecK, PolyVecL};

#[test]
fn test_fvec_roundtrip() {
    let mut z = PolyVecL::zero();
    let mut y = PolyVecK::zero();
    z.polys[0].coeffs[0] = 5;
    y.polys[0].coeffs[0] = -7;

    let fv = FVec::from_polyvecs(&z, &y);
    let mut z2 = PolyVecL::zero();
    let mut y2 = PolyVecK::zero();
    fv.round_to_polyvecs(&mut z2, &mut y2);

    println!("z2={}, y2={}", z2.polys[0].coeffs[0], y2.polys[0].coeffs[0]);
    assert_eq!(z2.polys[0].coeffs[0], 5);
    assert_eq!(y2.polys[0].coeffs[0], 8380410); // -7 mod Q
}
