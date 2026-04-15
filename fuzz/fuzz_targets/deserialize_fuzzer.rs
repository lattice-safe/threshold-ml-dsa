#![no_main]

use libfuzzer_sys::fuzz_target;
use threshold_ml_dsa::coordinator;
use threshold_ml_dsa::params::*;
use threshold_ml_dsa::poly::{PolyVecK, PolyVecL};
use threshold_ml_dsa::sign;

fn next_byte(data: &[u8], idx: &mut usize) -> u8 {
    if data.is_empty() {
        return 0;
    }
    let b = data[*idx % data.len()];
    *idx += 1;
    b
}

fn next_i32(data: &[u8], idx: &mut usize) -> i32 {
    let bytes = [
        next_byte(data, idx),
        next_byte(data, idx),
        next_byte(data, idx),
        next_byte(data, idx),
    ];
    i32::from_le_bytes(bytes)
}

fn fill_polyveck(data: &[u8], idx: &mut usize) -> PolyVecK {
    let mut out = PolyVecK::zero();
    for i in 0..K {
        for j in 0..N {
            out.polys[i].coeffs[j] = next_i32(data, idx);
        }
    }
    out
}

fn fill_polyvecl(data: &[u8], idx: &mut usize) -> PolyVecL {
    let mut out = PolyVecL::zero();
    for i in 0..L {
        for j in 0..N {
            out.polys[i].coeffs[j] = next_i32(data, idx);
        }
    }
    out
}

fuzz_target!(|data: &[u8]| {
    if data.is_empty() {
        return;
    }

    let mut idx = 0usize;

    // 1) Exercise pack/unpack robustness.
    let w = fill_polyveck(data, &mut idx);
    let packed = sign::pack_w_single(&w);
    let _ = sign::unpack_w_single(&packed);

    // Also try unpacking arbitrary-length attacker input directly.
    let _ = sign::unpack_w_single(data);

    // 2) Exercise combine() with adversarial commitments/responses.
    // Use a small config for speed while still hitting the full code path.
    let params = match get_threshold_params(2, 2) {
        Some(p) => p,
        None => return,
    };
    let k_reps = params.k_reps as usize;

    let mut pk = [0u8; PK_BYTES];
    for b in &mut pk {
        *b = next_byte(data, &mut idx);
    }

    let mut wfinals = Vec::with_capacity(k_reps);
    let mut zfinals = Vec::with_capacity(k_reps);
    for _ in 0..k_reps {
        wfinals.push(fill_polyveck(data, &mut idx));
    }
    for _ in 0..k_reps {
        zfinals.push(fill_polyvecl(data, &mut idx));
    }

    let msg = if idx < data.len() { &data[idx..] } else { b"" };
    let _ = coordinator::combine(&pk, msg, &wfinals, &zfinals, &params);
});
