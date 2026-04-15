//! E2E trace: step through one attempt to find where rejection occurs.

use threshold_ml_dsa::params::*;
use threshold_ml_dsa::poly::*;
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;

#[test]
fn trace_one_attempt() {
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    let seed = [42u8; 32];
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 2, 1).unwrap();

    let mut rng = StdRng::seed_from_u64(12345);
    let msg = b"Hello, threshold ML-DSA!";
    let active = [0u8, 1];
    let k_reps = sdk.params.k_reps as usize;

    let mut act: u8 = 0;
    for &id in &active {
        act |= 1 << id;
    }
    let session_id = [9u8; 32];

    // Round 1
    let mut rd1_hashes = Vec::new();
    let mut rd1_states = Vec::new();
    for &party_id in &active {
        let sk = &sdk.sks[party_id as usize];
        let (hash, st1) =
            threshold_ml_dsa::sign::round1(sk, &sdk.params, act, msg, &session_id, &mut rng)
                .unwrap();
        rd1_hashes.push(hash);
        rd1_states.push(st1);
    }

    // Round 2
    let mut rd2_reveals = Vec::new();
    let mut rd2_states = Vec::new();
    for (idx, &party_id) in active.iter().enumerate() {
        let sk = &sdk.sks[party_id as usize];
        let (reveal, st2) = threshold_ml_dsa::sign::round2(
            sk,
            act,
            msg,
            &session_id,
            &rd1_hashes,
            &rd1_hashes[idx],
            &rd1_states[idx],
            &sdk.params,
        )
        .unwrap();
        rd2_reveals.push(reveal);
        rd2_states.push(st2);
    }

    // Aggregate commitments
    let packed_size = threshold_ml_dsa::sign::pack_w_single_size();
    let mut all_reveals = Vec::new();
    for reveal in &rd2_reveals {
        let mut party_ws = Vec::new();
        for k in 0..k_reps {
            let start = k * packed_size;
            let end = start + packed_size;
            if end <= reveal.len() {
                party_ws.push(threshold_ml_dsa::sign::unpack_w_single(&reveal[start..end]));
            }
        }
        all_reveals.push(party_ws);
    }

    let wfinals = threshold_ml_dsa::coordinator::aggregate_commitments(&all_reveals, k_reps).unwrap();

    // Round 3
    let mut all_responses = Vec::new();
    for (st1, (&party_id, st2)) in rd1_states
        .into_iter()
        .zip(active.iter().zip(rd2_states.iter()))
    {
        let sk = &sdk.sks[party_id as usize];
        let zs = threshold_ml_dsa::sign::round3(
            sk,
            &wfinals,
            st1,
            st2,
            &sdk.params,
        )
        .unwrap();
        all_responses.push(zs);
    }

    let zfinals = threshold_ml_dsa::coordinator::aggregate_responses(&all_responses, k_reps).unwrap();

    // Check each slot
    for k in 0..k_reps {
        let z = &zfinals[k];
        let all_zero = z.polys.iter().all(|p| p.coeffs.iter().all(|&c| c == 0));
        let linf = z.l_inf_norm();

        if all_zero {
            println!("Slot {}: z is all-zero (hyperball rejected)", k);
            continue;
        }

        let norm_check = z.chknorm(GAMMA1 - BETA);
        println!(
            "Slot {}: L_inf={}, z_chknorm(gamma1-beta={})={}",
            k,
            linf,
            GAMMA1 - BETA,
            norm_check
        );

        if norm_check {
            println!("  -> REJECTED by z norm check");
            continue;
        }

        // If we got here, z passes norm check. Trace the hint computation.
        use dilithium::{
            packing::unpack_pk,
            poly::Poly as DPoly,
            polyvec::{
                matrix_expand, matrix_pointwise_montgomery, PolyVecK as DPolyVecK,
                PolyVecL as DPolyVecL,
            },
            polyvec::{polyveck_decompose, polyveck_pack_w1},
            rounding as d_rounding, ML_DSA_44,
        };

        let mode = ML_DSA_44;
        let mut rho = [0u8; 32];
        let mut t1 = DPolyVecK::default();
        unpack_pk(mode, &mut rho, &mut t1, &sdk.pk);

        let mut mat: [DPolyVecL; K] = core::array::from_fn(|_| DPolyVecL::default());
        matrix_expand(mode, &mut mat, &rho);

        // Convert z
        let mut zh_d = DPolyVecL::default();
        for i in 0..L {
            for j in 0..N {
                zh_d.vec[i].coeffs[j] = z.polys[i].coeffs[j];
            }
            zh_d.vec[i].ntt();
        }

        // Az
        let mut az = DPolyVecK::default();
        matrix_pointwise_montgomery(mode, &mut az, &mat, &zh_d);
        for i in 0..K {
            az.vec[i].reduce();
            az.vec[i].invntt_tomont();
        }

        // Challenge
        use sha3::{
            digest::{ExtendableOutput, Update, XofReader},
            Shake256,
        };
        use threshold_ml_dsa::poly::decompose;

        let mut h_tr = Shake256::default();
        h_tr.update(&sdk.pk);
        let mut tr = [0u8; TRBYTES];
        h_tr.finalize_xof().read(&mut tr);
        let mu = threshold_ml_dsa::sign::compute_mu(&tr, msg);

        let mut w1 = PolyVecK::zero();
        for i in 0..K {
            for j in 0..N {
                let a = wfinals[k].polys[i].coeffs[j];
                let a_pos = if a < 0 { a + Q } else { a };
                let (r1, _) = decompose(a_pos);
                w1.polys[i].coeffs[j] = r1;
            }
        }

        let mut h = Shake256::default();
        h.update(&mu);
        for poly in &w1.polys {
            let mut packed = [0u8; POLYW1_PACKEDBYTES];
            for kk in 0..N / 4 {
                packed[3 * kk] =
                    (poly.coeffs[4 * kk] as u8) | ((poly.coeffs[4 * kk + 1] as u8) << 6);
                packed[3 * kk + 1] =
                    ((poly.coeffs[4 * kk + 1] >> 2) as u8) | ((poly.coeffs[4 * kk + 2] as u8) << 4);
                packed[3 * kk + 2] =
                    ((poly.coeffs[4 * kk + 2] >> 4) as u8) | ((poly.coeffs[4 * kk + 3] as u8) << 2);
            }
            h.update(&packed);
        }
        let mut c_tilde = [0u8; CTILDEBYTES];
        h.finalize_xof().read(&mut c_tilde);

        // Compute challenge via dilithium polyveck_decompose+pack path as well.
        let mut w_ref = DPolyVecK::default();
        for i in 0..K {
            for j in 0..N {
                w_ref.vec[i].coeffs[j] = wfinals[k].polys[i].coeffs[j];
            }
        }
        let mut w1_ref = DPolyVecK::default();
        let mut w0_ref = DPolyVecK::default();
        polyveck_decompose(mode, &mut w1_ref, &mut w0_ref, &w_ref);
        let mut w1_packed_ref = vec![0u8; mode.k() * mode.polyw1_packedbytes()];
        polyveck_pack_w1(mode, &mut w1_packed_ref, &w1_ref);
        let mut h_ref = Shake256::default();
        h_ref.update(&mu);
        h_ref.update(&w1_packed_ref);
        let mut c_tilde_ref = [0u8; CTILDEBYTES];
        h_ref.finalize_xof().read(&mut c_tilde_ref);
        println!("  c_tilde_equal={}", c_tilde == c_tilde_ref);

        let mut ch = DPoly::default();
        DPoly::challenge(mode, &mut ch, &c_tilde);
        ch.ntt();

        // w_approx = az - ct1*2^d
        for i in 0..K {
            let mut ct1 = t1.vec[i].clone();
            ct1.shiftl();
            ct1.ntt();
            let ct1_copy = ct1.clone();
            DPoly::pointwise_montgomery(&mut ct1, &ct1_copy, &ch);
            ct1.reduce();
            ct1.invntt_tomont();

            let az_copy = az.vec[i].clone();
            DPoly::sub(&mut az.vec[i], &az_copy, &ct1);
            az.vec[i].caddq();
        }

        // delta = w - w_approx (centered)
        let mut delta_linf = 0i32;
        for i in 0..K {
            for j in 0..N {
                let mut d = wfinals[k].polys[i].coeffs[j] - az.vec[i].coeffs[j];
                if d > Q / 2 {
                    d -= Q;
                } else if d < -(Q / 2) {
                    d += Q;
                }
                let ad = d.abs();
                if ad > delta_linf {
                    delta_linf = ad;
                }
            }
        }
        println!("  delta_linf={} (gamma2={})", delta_linf, GAMMA2);

        // Count hints (same construction as coordinator::combine)
        let mut hint_pop = 0usize;
        let mut hint_fails = 0usize;
        let mut slot_valid = true;
        let mut w0 = PolyVecK::zero();
        for i in 0..K {
            for j in 0..N {
                let (_, r0) = d_rounding::decompose(mode, wfinals[k].polys[i].coeffs[j]);
                w0.polys[i].coeffs[j] = if r0 < 0 { r0 + Q } else { r0 };
                let target = w1.polys[i].coeffs[j];
                let approx = az.vec[i].coeffs[j];

                let mut f = (approx - wfinals[k].polys[i].coeffs[j]).rem_euclid(Q);
                if f > (Q - 1) / 2 {
                    f -= Q;
                }
                if f.unsigned_abs() > GAMMA2 as u32 {
                    slot_valid = false;
                    break;
                }

                let f_plus_q = if f < 0 { f + Q } else { f };
                let hint_input = (w0.polys[i].coeffs[j] + f_plus_q).rem_euclid(Q);
                let h = d_rounding::make_hint(mode, hint_input, target);
                if h {
                    hint_pop += 1;
                }
                if d_rounding::use_hint(mode, approx, h) != target {
                    hint_fails += 1;
                }
            }
            if !slot_valid || hint_pop > OMEGA {
                break;
            }
        }

        println!(
            "  hint_pop={}, hint_fails={}, slot_valid={}, OMEGA={}",
            hint_pop, hint_fails, slot_valid, OMEGA
        );
        if hint_pop <= OMEGA && hint_fails == 0 {
            println!("  -> WOULD PRODUCE VALID SIGNATURE!");
        }
    }
}
