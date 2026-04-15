use threshold_ml_dsa::fvec::{sample_hyperball, FVec};

#[test]
fn test_delta_norm() {
    let rhop = [0u8; 64];
    let mut fv = FVec::zero();
    sample_hyperball(&mut fv, 252833.0, 3.0, &rhop, 0);

    let mut max_e = 0.0;
    for i in 1024..2048 {
        if fv.coeffs[i].abs() > max_e {
            max_e = fv.coeffs[i].abs();
        }
    }
    println!("Max e component = {}", max_e);
}
