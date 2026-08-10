//! Benchmark suite for threshold ML-DSA (ePrint 2026/013).

use rand::rngs::StdRng;
use rand::SeedableRng;
use std::time::Instant;
use threshold_ml_dsa::{
    dkg,
    params::SK_BYTES,
    sdk::ThresholdMlDsa44Sdk,
    verify,
};

fn main() {
    println!("\n==========================================================");
    println!("        THRESHOLD ML-DSA (FIPS 204) BENCHMARK SUITE       ");
    println!("==========================================================");
    println!("Platform: Apple Silicon (macOS aarch64)");
    println!("Mode: Release (opt-level 3)\n");

    let mut rng = StdRng::seed_from_u64(42);
    let msg = b"Benchmark message for threshold ML-DSA signatures";

    // ─── 1. Key Generation Benchmarks ─────────────────────────────────────
    println!("--- 1. Fresh Key Generation (from_seed) ---");
    let configs: &[(u8, u8)] = &[
        (2, 2),
        (2, 3),
        (3, 3),
        (2, 4),
        (3, 4),
        (4, 4),
        (3, 5),
        (5, 5),
        (4, 6),
        (5, 6),
        (6, 6),
    ];

    for &(t, n) in configs {
        let iters = 100;
        let seed = [42u8; 32];
        let start = Instant::now();
        for _ in 0..iters {
            let _ = ThresholdMlDsa44Sdk::from_seed(&seed, t, n, 10).unwrap();
        }
        let elapsed = start.elapsed();
        let per_op = elapsed / iters as u32;
        println!("  Fresh Keygen ({}, {}): {:8.2?} / op", t, n, per_op);
    }

    // ─── 2. A Posteriori Key Sharing Benchmarks ───────────────────────────
    println!("\n--- 2. A Posteriori Key Sharing (from_existing_key) ---");
    let (pk_std, sk_std) = verify::keygen(&[99u8; 32]);
    let mut sk_fixed = [0u8; SK_BYTES];
    sk_fixed.copy_from_slice(&sk_std);

    for &(t, n) in &[(2, 2), (2, 3), (3, 4), (5, 6)] {
        let iters = 100;
        let start = Instant::now();
        for _ in 0..iters {
            let _ = ThresholdMlDsa44Sdk::from_existing_key(&sk_fixed, &[1u8; 32], t, n, 10).unwrap();
        }
        let elapsed = start.elapsed();
        let per_op = elapsed / iters as u32;
        println!("  A Posteriori Sharing ({}, {}): {:8.2?} / op", t, n, per_op);
    }

    // ─── 3. Distributed Key Generation (DKG) Benchmarks ──────────────────
    println!("\n--- 3. Distributed Key Generation (DKG) ---");
    for &(t, n) in &[(2, 3), (3, 4), (5, 6)] {
        let iters = 50;
        let start = Instant::now();
        for _ in 0..iters {
            let mut contribs = Vec::new();
            let mut commits = Vec::new();
            for i in 0..n {
                let (c, h) = dkg::dkg_round1(i, &mut rng).unwrap();
                contribs.push((i, c));
                commits.push((i, h));
            }
            let _ = ThresholdMlDsa44Sdk::from_dkg(&contribs, &commits, t, n, 10).unwrap();
        }
        let elapsed = start.elapsed();
        let per_op = elapsed / iters as u32;
        println!("  Full DKG Protocol ({}, {}): {:8.2?} / op", t, n, per_op);
    }

    // ─── 4. Threshold Signing Benchmarks ──────────────────────────────────
    println!("\n--- 4. Threshold Signing (threshold_sign) ---");
    for &(t, n) in configs {
        let seed = [t * 10 + n; 32];
        let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, t, n, 100).unwrap();
        let active: Vec<u8> = (0..t).collect();

        // Warmup
        let _ = sdk.threshold_sign(&active, msg, &mut rng).unwrap();

        let iters = 50;
        let start = Instant::now();
        for _ in 0..iters {
            let _ = sdk.threshold_sign(&active, msg, &mut rng).unwrap();
        }
        let elapsed = start.elapsed();
        let per_op = elapsed / iters as u32;
        println!("  Threshold Sign ({}, {}): {:8.2?} / op", t, n, per_op);
    }

    // ─── 5. FIPS 204 Signature Verification Benchmarks ─────────────────────
    println!("\n--- 5. Signature Verification (verify) ---");
    let sdk = ThresholdMlDsa44Sdk::from_seed(&[1u8; 32], 2, 3, 10).unwrap();
    let sig_th = sdk.threshold_sign(&[0, 1], msg, &mut rng).unwrap();

    let iters = 500;
    let start = Instant::now();
    for _ in 0..iters {
        let _ = sdk.verify(msg, &sig_th);
    }
    let elapsed = start.elapsed();
    let per_op = elapsed / iters as u32;
    println!("  Threshold Sig FIPS 204 Verify (Level 2 / 44): {:8.2?} / op", per_op);

    // Standard ML-DSA-65 & 87 verify
    let (pk65, sk65) = verify::keygen_65(&[2u8; 32]);
    let sig65 = verify::sign_standard_65(msg, &sk65, &[3u8; 32]);
    let start = Instant::now();
    for _ in 0..iters {
        let _ = verify::verify_65(&sig65, msg, &pk65);
    }
    println!("  Standard FIPS 204 Verify (Level 3 / 65):        {:8.2?} / op", start.elapsed() / iters as u32);

    let (pk87, sk87) = verify::keygen_87(&[4u8; 32]);
    let sig87 = verify::sign_standard_87(msg, &sk87, &[5u8; 32]);
    let start = Instant::now();
    for _ in 0..iters {
        let _ = verify::verify_87(&sig87, msg, &pk87);
    }
    println!("  Standard FIPS 204 Verify (Level 5 / 87):        {:8.2?} / op", start.elapsed() / iters as u32);

    println!("\n==========================================================");
}
