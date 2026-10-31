//! High-level SDK API for threshold ML-DSA-44.
//!
//! v0.3.0: This module is being rewritten to use the paper-faithful keygen
//! from `rss::keygen_from_seed`. The old `from_secret_key` API is removed.

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
#[cfg(feature = "std")]
use std::vec::Vec;

use crate::error::Error;
use crate::nonce::InMemoryNonceGuard;
use crate::params::{PK_BYTES, ThresholdParams, get_threshold_params, SIG_BYTES};
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
///
/// Fields are `pub(crate)` to prevent external mutation of secret
/// key material and parameter invariants. Use accessor methods
/// for read-only inspection.
pub struct ThresholdMlDsa44Sdk {
    /// Packed public key (ρ ‖ t₁)
    pub(crate) pk: [u8; PK_BYTES],
    /// Per-party private keys
    pub(crate) sks: Vec<rss::ThresholdPrivateKey>,
    /// Threshold parameters (T, N, K, r, r₁, ν)
    pub(crate) params: ThresholdParams,
    /// Maximum full-protocol retries
    pub(crate) max_retries: usize,
    /// Round 1 nonce guard shared by all in-process parties
    pub(crate) guard: InMemoryNonceGuard,
}

impl ThresholdMlDsa44Sdk {
    /// Returns a reference to the packed public key.
    #[must_use]
    pub fn pk(&self) -> &[u8; PK_BYTES] {
        &self.pk
    }

    /// Returns the number of parties (N).
    #[must_use]
    pub fn num_parties(&self) -> usize {
        self.sks.len()
    }

    /// Returns a reference to the threshold parameters.
    #[must_use]
    pub fn params(&self) -> &ThresholdParams {
        &self.params
    }

    /// Returns a reference to a party's private key by index.
    ///
    /// Intended for protocol-level integration tests that exercise
    /// the low-level 3-round API. Production callers should use
    /// `threshold_sign()` instead.
    #[must_use]
    pub fn party_key(&self, index: usize) -> Option<&rss::ThresholdPrivateKey> {
        self.sks.get(index)
    }


    /// Create a threshold SDK from a 32-byte seed.
    ///
    /// Generates fresh threshold keys using the paper's keygen (Fig. 5, USENIX version).
    /// Returns the SDK instance ready for signing.
    pub fn from_seed(seed: &[u8; 32], t: u8, n: u8, max_retries: usize) -> Result<Self, Error> {
        let params = get_threshold_params(t, n).ok_or(Error::InvalidParameters)?;

        if max_retries == 0 {
            return Err(Error::InvalidParameters);
        }

        let (pk, sks) = rss::keygen_from_seed(seed, &params)?;

        Ok(Self {
            pk,
            sks,
            params,
            max_retries,
            guard: InMemoryNonceGuard::new(),
        })
    }

    /// Create a threshold SDK by splitting an existing standard ML-DSA-44
    /// secret key with the paper's Gaussian a posteriori sharing (§3.4).
    /// The public key is unchanged.
    ///
    /// Short shares of a fixed key leak a hint about it to T−1 parties; the
    /// paper estimates the loss at up to 12 bits of security (heuristic).
    /// **Share a key at most once**: every further sharing adds another hint.
    /// See [`crate::aposteriori`].
    pub fn from_existing_key<R: RngCore + CryptoRng>(
        sk_bytes: &[u8],
        t: u8,
        n: u8,
        max_retries: usize,
        rng: &mut R,
    ) -> Result<Self, Error> {
        let params = get_threshold_params(t, n).ok_or(Error::InvalidParameters)?;

        if max_retries == 0 {
            return Err(Error::InvalidParameters);
        }

        let (pk, sks) = crate::aposteriori::share_existing_key(sk_bytes, &params, rng)?;

        Ok(Self {
            pk,
            sks,
            params,
            max_retries,
            guard: InMemoryNonceGuard::new(),
        })
    }

    /// Create an SDK from existing party keys (e.g. the outputs of
    /// [`crate::dkg`] or keys parsed with
    /// [`rss::ThresholdPrivateKey::from_bytes`]).
    ///
    /// `sks` must hold exactly one key for each party ID `0..N`, all for the
    /// same (T, N) and public key. Note that holding every party's key in one
    /// process means this process can sign alone.
    pub fn from_party_keys(
        pk: [u8; PK_BYTES],
        mut sks: Vec<rss::ThresholdPrivateKey>,
        max_retries: usize,
    ) -> Result<Self, Error> {
        if max_retries == 0 || sks.is_empty() {
            return Err(Error::InvalidParameters);
        }
        let (t, n) = sks[0].threshold();
        let params = get_threshold_params(t, n).ok_or(Error::InvalidParameters)?;
        if sks.len() != n as usize {
            return Err(Error::InvalidParameters);
        }
        sks.sort_by_key(rss::ThresholdPrivateKey::id);
        let tr = rss::tr_from_pk(&pk);
        for (i, sk) in sks.iter().enumerate() {
            if sk.id() as usize != i
                || sk.threshold() != (t, n)
                || sk.rho != pk[..32]
                || sk.tr != tr
            {
                return Err(Error::InvalidShare);
            }
        }
        // Every copy of a subset share must agree, and Σ_I s_I must give pk.
        let mut s1h_sum = PolyVecL::zero();
        let mut s2_sum = PolyVecK::zero();
        for mask in rss::subset_masks(n, t) {
            let mut holders = sks.iter().filter_map(|sk| sk.shares.get(&mask));
            let first = holders.next().ok_or(Error::InvalidShare)?;
            let mut diff = 0i32;
            for other in holders {
                for (a, b) in first.s1.polys.iter().chain(&first.s2.polys).zip(other.s1.polys.iter().chain(&other.s2.polys)) {
                    for (x, y) in a.coeffs.iter().zip(&b.coeffs) {
                        diff |= x ^ y;
                    }
                }
            }
            if diff != 0 {
                return Err(Error::InvalidShare);
            }
            s1h_sum.add_assign(&first.s1h);
            s2_sum.add_assign(&first.s2);
            s1h_sum.reduce();
            s2_sum.reduce();
        }
        let consistent = rss::pk_from_t(&pk_rho(&pk), &rss::compute_t(&pk_rho(&pk), &s1h_sum, &s2_sum)) == pk;
        s1h_sum.zeroize();
        s2_sum.zeroize();
        if !consistent {
            return Err(Error::InvalidShare);
        }
        Ok(Self {
            pk,
            sks,
            params,
            max_retries,
            guard: InMemoryNonceGuard::new(),
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
        &mut self,
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

        // Validate and build active bitmask.
        let mut act: u8 = 0;
        let mut prev: Option<u8> = None;
        for &id in active {
            if id >= self.params.n {
                return Err(Error::InvalidParameters);
            }
            if let Some(p) = prev {
                // Enforce strict ordering to reject duplicates and ambiguous sets.
                if id <= p {
                    return Err(Error::InvalidParameters);
                }
            }
            act |= 1 << id;
            prev = Some(id);
        }

        for _attempt in 0..self.max_retries {
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
                let (hash, st1) =
                    sign::round1(sk, &self.params, act, msg, &session_id, rng, &mut self.guard)?;
                rd1_hashes.push(hash);
                rd1_states.push(st1);
            }

            // ── Round 2: Each party reveals commitments and computes μ ──
            let mut rd2_reveals: Vec<Vec<u8>> = Vec::with_capacity(t);
            let mut rd2_states: Vec<sign::StRound2> = Vec::with_capacity(t);

            for (idx, (&party_id, st1)) in active.iter().zip(rd1_states).enumerate() {
                let sk = &self.sks[party_id as usize];
                let (reveal, st2) = sign::round2(
                    sk,
                    act,
                    msg,
                    &session_id,
                    &rd1_hashes,
                    &rd1_hashes[idx],
                    st1,
                    &self.params,
                )?;
                rd2_reveals.push(reveal);
                rd2_states.push(st2);
            }

            // ── Aggregate commitments ──
            // Parse each party's reveal into K PolyVecK commitments
            let active_ids: Vec<u8> = active.to_vec();
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

            let wfinals = coordinator::aggregate_commitments(&all_reveals, k_reps)?;

            // ── Round 3: Each party computes K responses ──
            // Each party independently verifies reveals and obtains its own
            // consumable VerifiedReveals witness. round3 takes the witness by
            // value to prevent cross-session replay, so each party needs a
            // fresh one (mirroring a real distributed deployment).
            // round3 consumes the Round 2 state (and with it the nonce state).
            let mut all_responses: Vec<Vec<PolyVecL>> = Vec::with_capacity(t);
            for (&party_id, st2) in active.iter().zip(rd2_states) {
                let sk = &self.sks[party_id as usize];

                // Each party verifies all reveals independently (as in a
                // distributed deployment) and gets its own consumable witness.
                let verified = sign::verify_all_round2_reveals(
                    &sk.tr,
                    &active_ids,
                    act,
                    msg,
                    &session_id,
                    &rd2_reveals,
                    &rd1_hashes,
                    k_reps,
                )?;

                // Rejected slots come back zeroed; rejection shows up in combine.
                let zs = sign::round3(sk, &wfinals, st2, &self.params, verified)?;
                all_responses.push(zs);
            }

            let zfinals = coordinator::aggregate_responses(&all_responses, k_reps)?;

            // ── Combine: try each of K slots ──
            // combine() only returns a signature that passes FIPS 204
            // verification under `self.pk`.
            if let Ok(sig) = coordinator::combine(&self.pk, msg, &wfinals, &zfinals, &self.params) {
                return Ok(sig);
            }
        }

        Err(Error::InsufficientResponses)
    }

    /// Verify a signature using the standard ML-DSA-44 verifier.
    #[must_use] 
    pub fn verify(&self, msg: &[u8], sig: &[u8]) -> bool {
        verify::verify(sig, msg, &self.pk)
    }
}

/// ρ, the first 32 bytes of a packed public key.
fn pk_rho(pk: &[u8; PK_BYTES]) -> [u8; 32] {
    let mut rho = [0u8; 32];
    rho.copy_from_slice(&pk[..32]);
    rho
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

    #[test]
    fn test_threshold_sign_rejects_duplicate_or_unsorted_active_set() {
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let seed = [99u8; 32];
        let mut sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 3, 2).unwrap();
        let mut rng = StdRng::seed_from_u64(1);
        let msg = b"active set validation";

        assert_eq!(
            sdk.threshold_sign(&[0, 0], msg, &mut rng),
            Err(Error::InvalidParameters)
        );
        assert_eq!(
            sdk.threshold_sign(&[2, 1], msg, &mut rng),
            Err(Error::InvalidParameters)
        );
        assert_eq!(
            sdk.threshold_sign(&[0, 3], msg, &mut rng),
            Err(Error::InvalidParameters)
        );
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn constructors_reject_bad_input() {
        let mut rng = StdRng::seed_from_u64(9);
        let (_, sk) = verify::keygen(&[9u8; 32]);
        assert_eq!(ThresholdMlDsa44Sdk::from_existing_key(&sk, 2, 3, 0, &mut rng).err(), Some(Error::InvalidParameters));
        assert_eq!(ThresholdMlDsa44Sdk::from_existing_key(&sk, 2, 7, 1, &mut rng).err(), Some(Error::InvalidParameters));
        assert_eq!(ThresholdMlDsa44Sdk::from_seed(&[0u8; 32], 2, 3, 0).err(), Some(Error::InvalidParameters));

        let params = get_threshold_params(2, 3).unwrap();
        let (pk, sks) = rss::keygen_from_seed(&[10u8; 32], &params).unwrap();
        assert_eq!(ThresholdMlDsa44Sdk::from_party_keys(pk, Vec::new(), 1).err(), Some(Error::InvalidParameters));
        let (_, mut two) = rss::keygen_from_seed(&[10u8; 32], &params).unwrap();
        two.pop();
        assert_eq!(ThresholdMlDsa44Sdk::from_party_keys(pk, two, 1).err(), Some(Error::InvalidParameters));
        let (other_pk, _) = rss::keygen_from_seed(&[11u8; 32], &params).unwrap();
        assert_eq!(ThresholdMlDsa44Sdk::from_party_keys(other_pk, sks, 1).err(), Some(Error::InvalidShare));

        let sdk = ThresholdMlDsa44Sdk::from_seed(&[12u8; 32], 2, 3, 1).unwrap();
        assert_eq!(sdk.num_parties(), 3);
        assert_eq!(sdk.params().t, 2);
        assert!(sdk.party_key(3).is_none());
    }
}

#[cfg(test)]
mod inconsistent_key_tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    // `from_party_keys` refuses such key sets; these tests build them
    // directly to check that signing fails loudly instead of retrying.

    #[test]
    fn swapped_party_keys_fail_round2() {
        let mut sdk = ThresholdMlDsa44Sdk::from_seed(&[13u8; 32], 2, 2, 3).unwrap();
        sdk.sks.swap(0, 1);
        let mut rng = StdRng::seed_from_u64(13);
        assert_eq!(sdk.threshold_sign(&[0, 1], b"m", &mut rng).err(), Some(Error::InvalidShare));
    }

    #[test]
    fn mismatched_tr_fails_reveal_verification() {
        let mut sdk = ThresholdMlDsa44Sdk::from_seed(&[14u8; 32], 2, 2, 3).unwrap();
        sdk.sks[1].tr[0] ^= 1;
        let mut rng = StdRng::seed_from_u64(14);
        assert_eq!(sdk.threshold_sign(&[0, 1], b"m", &mut rng).err(), Some(Error::InvalidShare));
    }
}

#[cfg(test)]
mod signing_tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn party_keys_sign_and_verify() {
        let params = get_threshold_params(2, 3).unwrap();
        let (pk, mut sks) = rss::keygen_from_seed(&[16u8; 32], &params).unwrap();
        sks.reverse(); // order does not matter
        let mut sdk = ThresholdMlDsa44Sdk::from_party_keys(pk, sks, 50).unwrap();
        assert_eq!(sdk.pk(), &pk);
        let mut rng = StdRng::seed_from_u64(16);
        let sig = sdk.threshold_sign(&[0, 2], b"msg", &mut rng).unwrap();
        assert!(sdk.verify(b"msg", &sig));
        assert!(!sdk.verify(b"other", &sig));
        assert_eq!(sdk.threshold_sign(&[0], b"msg", &mut rng).err(), Some(Error::InvalidParameters));
    }

    #[test]
    fn gives_up_after_max_retries() {
        // One attempt per call: about 40% of calls at (2, 2) hit rejection
        // in every slot and must report InsufficientResponses.
        let mut sdk = ThresholdMlDsa44Sdk::from_seed(&[17u8; 32], 2, 2, 1).unwrap();
        let mut rng = StdRng::seed_from_u64(17);
        let failures = (0..40u32)
            .filter(|i| sdk.threshold_sign(&[0, 1], &i.to_le_bytes(), &mut rng).err() == Some(Error::InsufficientResponses))
            .count();
        assert!(failures > 0);
    }
}

#[cfg(test)]
mod consistency_tests {
    use super::*;

    #[test]
    fn from_party_keys_rejects_inconsistent_shares() {
        let params = get_threshold_params(2, 3).unwrap();
        // One holder's copy of a share differs.
        let (pk, mut sks) = rss::keygen_from_seed(&[30u8; 32], &params).unwrap();
        let mask = *sks[0].share_masks().next().unwrap();
        sks[0].shares.get_mut(&mask).unwrap().s1.polys[0].coeffs[0] += 1;
        assert_eq!(ThresholdMlDsa44Sdk::from_party_keys(pk, sks, 1).err(), Some(Error::InvalidShare));

        // All copies changed consistently: shares no longer sum to pk. (A
        // change to s₁ moves t₁ everywhere; a small change to s₂ may only
        // move t₀, which pk does not contain.)
        let (pk, mut sks) = rss::keygen_from_seed(&[31u8; 32], &params).unwrap();
        for sk in &mut sks {
            if let Some(sh) = sk.shares.get_mut(&mask) {
                sh.s1.polys[0].coeffs[0] += 1;
                sh.s1h = sh.s1.clone();
                sh.s1h.reduce();
                sh.s1h.ntt();
            }
        }
        assert_eq!(ThresholdMlDsa44Sdk::from_party_keys(pk, sks, 1).err(), Some(Error::InvalidShare));
    }
}
