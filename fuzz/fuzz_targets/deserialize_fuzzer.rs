#![no_main]

use libfuzzer_sys::fuzz_target;
use threshold_ml_dsa::coordinator::aggregate_responses;
use threshold_ml_dsa::sign::PartialSignature;
use threshold_ml_dsa::poly::{PolyVecL, PolyVecK};
use threshold_ml_dsa::params::*;

fuzz_target!(|data: &[u8]| {
    // We need enough data to populate w, c_tilde, pk, and at least some partials
    let min_len = POLYK_SIZE + CTILDEBYTES + PK_BYTES + POLYL_SIZE;
    if data.len() < min_len {
        return;
    }

    let mut offset = 0;
    
    // We mock the structures with the raw fuzzer inputs
    let mut w = PolyVecK::zero();
    let mut z = PolyVecL::zero();
    
    // Fuzz dummy w
    for i in 0..K {
        for j in 0..N {
            w.polys[i].coeffs[j] = (data[offset % data.len()] as i32) * 1000;
            offset += 1;
        }
    }

    let mut c_tilde = vec![0u8; CTILDEBYTES];
    c_tilde.copy_from_slice(&data[offset..offset+CTILDEBYTES]);
    offset += CTILDEBYTES;

    let mut pk = vec![0u8; PK_BYTES];
    if offset + PK_BYTES <= data.len() {
        pk.copy_from_slice(&data[offset..offset+PK_BYTES]);
        offset += PK_BYTES;
    } else {
        return;
    }

    // Generate valid length PartialSignatures with raw z garbage
    for i in 0..L {
        for j in 0..N {
            z.polys[i].coeffs[j] = (data[offset % data.len()] as i32) * 100_000;
            offset += 1;
        }
    }
    
    let partials = vec![PartialSignature { party_id: 1, z }];

    // Run the aggregation core (which unpacks PK and evaluates MakeHint bounds)
    // The goal here is that it MUST NOT panic. It should just gracefully Err(InvalidSignature).
    let _ = aggregate_responses(&partials, &w, &c_tilde, 1, &pk);
});

// Size helpers since poly structure isn't direct
const POLYK_SIZE: usize = K * N;
const POLYL_SIZE: usize = L * N;
