//! High-level SDK API for threshold ML-DSA-44.
//!
//! v0.3.0: This module is being rewritten to use the paper-faithful keygen
//! from `rss::keygen_from_seed`. The old `from_secret_key` API is removed.

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
#[cfg(feature = "std")]
use std::vec::Vec;

use crate::error::Error;
use crate::params::*;
use crate::poly::{PolyVecK, PolyVecL};
use crate::rss;
use crate::verify;
use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use zeroize::Zeroize;

/// High-level threshold ML-DSA-44 SDK.
///
/// This provides a one-call interface for threshold key generation
/// and signing, following ePrint 2026/013 exactly.
pub struct ThresholdMlDsa44Sdk {
    /// Packed public key (ρ ‖ t₁)
    pub pk: [u8; PK_BYTES],
    /// Per-party private keys
    pub sks: Vec<rss::ThresholdPrivateKey>,
    /// Threshold parameters (T, N, K, r, r₁, ν)
    pub params: ThresholdParams,
    /// Maximum full-protocol retries
    pub max_retries: usize,
}

impl ThresholdMlDsa44Sdk {
    /// Create a threshold SDK from a 32-byte seed.
    ///
    /// Generates fresh threshold keys using the paper's keygen (Figure 4).
    /// Returns the SDK instance ready for signing.
    pub fn from_seed(seed: &[u8; 32], t: u8, n: u8, max_retries: usize) -> Result<Self, Error> {
        let params = get_threshold_params(t, n).ok_or(Error::InvalidParameters)?;

        if max_retries == 0 {
            return Err(Error::InvalidParameters);
        }

        let (pk, sks) = rss::keygen_from_seed(seed, &params);

        Ok(Self {
            pk,
            sks,
            params,
            max_retries,
        })
    }

    /// Sign a message using threshold signing.
    ///
    /// # Arguments  
    /// * `active` — sorted list of active signer IDs (length = T)
    /// * `msg` — message to sign
    /// * `rng` — CSPRNG for hedged nonce generation
    ///
    /// # Returns
    /// A valid FIPS 204 signature on success.
    pub fn threshold_sign<R: RngCore + CryptoRng>(
        &self,
        active: &[u8],
        msg: &[u8],
        rng: &mut R,
    ) -> Result<[u8; SIG_BYTES], Error> {
        use crate::coordinator;
        use crate::sign;

        if active.len() != self.params.t as usize {
            return Err(Error::InvalidParameters);
        }

        let t = self.params.t as usize;
        let k_reps = self.params.k_reps as usize;

        // Build active bitmask
        let mut act: u8 = 0;
        for &id in active {
            act |= 1 << id;
        }

        'attempt: for _attempt in 0..self.max_retries {
            // Session binding for this signing attempt.
            let mut session_entropy = [0u8; 32];
            rng.fill_bytes(&mut session_entropy);
            let mut h = Shake256::default();
            h.update(b"th-ml-dsa-session-v1");
            h.update(&session_entropy);
            h.update(&self.pk);
            h.update(&[act]);
            h.update(msg);
            let mut session_id = [0u8; 32];
            h.finalize_xof().read(&mut session_id);
            session_entropy.zeroize();

            // ── Round 1: Each party generates K commitments ──
            let mut rd1_hashes: Vec<[u8; 32]> = Vec::with_capacity(t);
            let mut rd1_states: Vec<sign::StRound1> = Vec::with_capacity(t);

            for &party_id in active {
                let sk = &self.sks[party_id as usize];
                // Clone rng bytes for each party (in production, each party has its own RNG)
                let (hash, st1) = sign::round1(sk, &self.params, act, msg, &session_id, rng)?;
                rd1_hashes.push(hash);
                rd1_states.push(st1);
            }

            // ── Round 2: Each party reveals commitments and computes μ ──
            let mut rd2_reveals: Vec<Vec<u8>> = Vec::with_capacity(t);
            let mut rd2_states: Vec<sign::StRound2> = Vec::with_capacity(t);

            for (idx, &party_id) in active.iter().enumerate() {
                let sk = &self.sks[party_id as usize];
                let (reveal, st2) = match sign::round2(
                    sk,
                    act,
                    msg,
                    &session_id,
                    &rd1_hashes,
                    &rd1_hashes[idx],
                    &rd1_states[idx],
                    &self.params,
                ) {
                    Ok(v) => v,
                    Err(_) => continue 'attempt,
                };
                rd2_reveals.push(reveal);
                rd2_states.push(st2);
            }

            // Verify each reveal against its Round 1 commitment hash.
            for (idx, &party_id) in active.iter().enumerate() {
                let tr = &self.sks[party_id as usize].tr;
                if !sign::verify_round2_reveal(
                    tr,
                    party_id,
                    act,
                    msg,
                    &session_id,
                    &rd2_reveals[idx],
                    &rd1_hashes[idx],
                ) {
                    continue 'attempt;
                }
            }

            // ── Aggregate commitments ──
            // Parse each party's reveal into K PolyVecK commitments
            let packed_size = sign::pack_w_single_size();
            let mut all_reveals: Vec<Vec<PolyVecK>> = Vec::with_capacity(t);
            for reveal in &rd2_reveals {
                let mut party_ws: Vec<PolyVecK> = Vec::with_capacity(k_reps);
                for k in 0..k_reps {
                    let start = k * packed_size;
                    let end = start + packed_size;
                    if end <= reveal.len() {
                        party_ws.push(sign::unpack_w_single(&reveal[start..end]));
                    }
                }
                all_reveals.push(party_ws);
            }

            let wfinals = coordinator::aggregate_commitments(&all_reveals, k_reps);

            // ── Round 3: Each party computes K responses ──
            let mut all_responses: Vec<Vec<PolyVecL>> = Vec::with_capacity(t);
            for (idx, &party_id) in active.iter().enumerate() {
                let sk = &self.sks[party_id as usize];
                let zs = sign::round3(
                    sk,
                    &wfinals,
                    &rd1_states[idx],
                    &rd2_states[idx],
                    &self.params,
                );
                all_responses.push(zs);
            }

            let zfinals = coordinator::aggregate_responses(&all_responses, k_reps);

            // ── Combine: try each of K slots ──
            match coordinator::combine(&self.pk, msg, &wfinals, &zfinals, &self.params) {
                Ok(sig) => {
                    // Fail-closed: verify before returning
                    if verify::verify(&sig, msg, &self.pk) {
                        return Ok(sig);
                    }
                    // Signature didn't verify — retry
                    continue;
                }
                Err(_) => continue,
            }
        }

        Err(Error::InsufficientResponses)
    }

    /// Verify a signature using the standard ML-DSA-44 verifier.
    pub fn verify(&self, msg: &[u8], sig: &[u8]) -> bool {
        verify::verify(sig, msg, &self.pk)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sdk_creation() {
        let seed = [42u8; 32];
        let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 3, 10).unwrap();
        assert_eq!(sdk.sks.len(), 3);
        assert_eq!(sdk.params.t, 2);
        assert_eq!(sdk.params.n, 3);
    }

    #[test]
    fn test_sdk_invalid_params() {
        let seed = [42u8; 32];
        // T < 2
        assert!(ThresholdMlDsa44Sdk::from_seed(&seed, 1, 2, 10).is_err());
        // T > N
        assert!(ThresholdMlDsa44Sdk::from_seed(&seed, 3, 2, 10).is_err());
        // N > 6
        assert!(ThresholdMlDsa44Sdk::from_seed(&seed, 2, 7, 10).is_err());
        // max_retries = 0
        assert!(ThresholdMlDsa44Sdk::from_seed(&seed, 2, 2, 0).is_err());
    }
}
