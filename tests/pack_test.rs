use threshold_ml_dsa::poly::PolyVecK;

#[test]
fn test_pack_overflow() {
    let mut w = PolyVecK::zero();
    // 8380417 + 10000 = 8390417
    w.polys[0].coeffs[0] = 8390417;
    let buf = threshold_ml_dsa::sign::pack_w_single(&w);
    let w2 = threshold_ml_dsa::sign::unpack_w_single(&buf);
    println!("w2 = {}", w2.polys[0].coeffs[0]);
}
