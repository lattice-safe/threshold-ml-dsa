use threshold_ml_dsa::verify;
use threshold_ml_dsa::sign::Party;
use threshold_ml_dsa::coordinator;
use threshold_ml_dsa::poly::{MatrixA, PolyVecK, PolyVecL};
use threshold_ml_dsa::rss;
use threshold_ml_dsa::params::*;
use rand::thread_rng;

fn main() {
    let mut rng = thread_rng();
    let seed = [42u8; 32];
    let (pk, sk) = verify::keygen(&seed);
    
    // Actually, in our scheme we need s1 and s2 to distribute.
    println!("OK");
}
