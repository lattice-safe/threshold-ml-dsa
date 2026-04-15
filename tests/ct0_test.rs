use threshold_ml_dsa::poly::Poly;
use threshold_ml_dsa::params::*;

#[test]
fn test_ct0_norm() {
    // Generate a random t0 in [-4096, 4096]
    use rand::Rng;
    let mut rng = rand::thread_rng();
    
    let mut max_all = 0;
    for _ in 0..1000 {
        let mut t0 = Poly::zero();
        for i in 0..N {
            t0.coeffs[i] = rng.gen_range(-4096..=4096);
        }
        
        let mut c_tilde = [0u8; CTILDEBYTES];
        rng.fill(&mut c_tilde);
        
        let mut c = Poly::zero();
        Poly::challenge(&mut c, &c_tilde);
        
        c.ntt();
        t0.ntt();
        let mut prod = Poly::zero();
        Poly::pointwise_montgomery(&mut prod, &c, &t0);
        prod.invntt_tomont();
        
        let mut m = 0;
        for i in 0..N {
            let mut v = prod.coeffs[i];
            if v > Q / 2 { v -= Q; }
            else if v < -Q/2 { v += Q; }
            if v.abs() > m { m = v.abs(); }
        }
        if m > max_all { max_all = m; }
    }
    println!("Max c*t0 norm over 1000 trials: {}", max_all);
}
