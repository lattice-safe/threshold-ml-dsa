//! Party (Signer) state machine for the threshold ML-DSA signing protocol.
//!
//! Implements the hardened 4-round flow:
//!
//! 0. **Pre-commit:** broadcast binding hash H(w_i ‖ party_id)
//! 1. **Reveal:** reveal full commitment w_i
//! 2. **Challenge:** coordinator derives c̃ from aggregate transcript
//! 3. **Sign:** return z_i = c·s_i + y_i with local hyperball rejection
//!
//! For transcript consistency, each round is bound to:
//! - the same signer set
//! - a coordinator-provided session context

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use crate::error::Error;
use crate::params::*;
use crate::poly::{MatrixA, Poly, PolyVecK, PolyVecL};
use crate::rss::{enumerate_subsets, PartyKeyShare};
use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use zeroize::Zeroize;

/// A binding hash for a commitment (Round 0 — pre-commitment).
///
/// This is broadcast before the full commitment w_i is revealed,
/// preventing a corrupt coordinator from grinding challenges by
/// selectively including/excluding parties (ADV-1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CommitmentHash {
    /// H(w_i_packed ‖ party_id) — 32 bytes.
    pub hash: [u8; 32],
    /// The party that produced this commitment.
    pub party_id: usize,
}

/// A partial commitment from Round 1 of the signing protocol.
///
/// Contains w_i = A·NTT(y_i) — the party's partial contribution
/// to the aggregate commitment.
#[derive(Clone, Zeroize)]
#[zeroize(drop)]
pub struct Commitment {
    /// The partial commitment vector w_i ∈ ℤ_q^K (in normal domain).
    pub w: PolyVecK,
    /// The binding hash for this commitment (ADV-1).
    #[zeroize(skip)]
    pub binding_hash: [u8; 32],
    /// The party that produced this commitment.
    #[zeroize(skip)]
    pub party_id: usize,
}

/// A partial signature from Round 3 of the signing protocol.
///
/// Contains the party's response z_i that passed the hyperball rejection check,
/// bound to the session via c_tilde_hash (ADV-6).
#[derive(Clone, Debug, Zeroize)]
#[zeroize(drop)]
pub struct PartialSignature {
    /// Party index.
    #[zeroize(skip)]
    pub party_id: usize,
    /// The partial response z_i ∈ ℤ_q^L.
    pub z: PolyVecL,
    /// Session binding: H(c̃ ‖ party_id) to prevent replay (ADV-6).
    #[zeroize(skip)]
    pub session_binding: [u8; 32],
}

/// A signing party in the threshold ML-DSA protocol.
///
/// Maintains per-round state across precommit/reveal/sign steps.
pub struct Party {
    /// Party index in [0, N).
    pub id: usize,

    /// This party's RSS key share.
    key_share: PartyKeyShare,

    /// Round-scoped additive share derived from RSS for the active signer set.
    signing_s1_share: Option<PolyVecL>,

    /// Local nonce counter mixed into hedged nonce derivation.
    nonce_counter: u64,

    /// The public matrix A (in NTT domain).
    a_hat: MatrixA,

    /// Ephemeral masking vector y_i (Round 1 state). Zeroized after use.
    y: Option<PolyVecL>,

    /// Partial commitment w_i (Round 1 output, carried to Round 3).
    w: Option<PolyVecK>,
}

impl Drop for Party {
    fn drop(&mut self) {
        self.key_share.zeroize();
        if let Some(ref mut s1) = self.signing_s1_share {
            s1.zeroize();
        }
        if let Some(ref mut y) = self.y {
            y.zeroize();
        }
    }
}

/// Compute a binding hash for a commitment: H(w_coefficients ‖ party_id).
fn compute_commitment_hash(w: &PolyVecK, party_id: usize) -> [u8; 32] {
    let mut hasher = Shake256::default();
    for i in 0..K {
        for j in 0..N {
            hasher.update(&w.polys[i].coeffs[j].to_le_bytes());
        }
    }
    hasher.update(&(party_id as u64).to_le_bytes());
    let mut hash = [0u8; 32];
    hasher.finalize_xof().read(&mut hash);
    hash
}

/// Compute a session binding tag: H(c̃ ‖ party_id).
fn compute_session_binding(c_tilde: &[u8], party_id: usize) -> [u8; 32] {
    let mut hasher = Shake256::default();
    hasher.update(c_tilde);
    hasher.update(&(party_id as u64).to_le_bytes());
    let mut binding = [0u8; 32];
    hasher.finalize_xof().read(&mut binding);
    binding
}

/// Verify a commitment against its pre-committed hash (ADV-1).
pub fn verify_commitment(commitment: &Commitment, expected_hash: &CommitmentHash) -> bool {
    if commitment.party_id != expected_hash.party_id {
        return false;
    }
    let actual_hash = compute_commitment_hash(&commitment.w, commitment.party_id);
    // Constant-time comparison
    use subtle::ConstantTimeEq;
    actual_hash.ct_eq(&expected_hash.hash).into()
}

/// Verify a partial signature's session binding (ADV-6).
pub fn verify_session_binding(partial: &PartialSignature, c_tilde: &[u8]) -> bool {
    let expected = compute_session_binding(c_tilde, partial.party_id);
    use subtle::ConstantTimeEq;
    partial.session_binding.ct_eq(&expected).into()
}

impl Party {
    /// Create a new signing party from its RSS key share and the public matrix.
    pub fn new(key_share: &PartyKeyShare, a_hat: MatrixA) -> Self {
        Party {
            id: key_share.party_id,
            key_share: key_share.clone(),
            signing_s1_share: None,
            nonce_counter: 0,
            a_hat,
            y: None,
            w: None,
        }
    }

    /// Return the `(n, t)` sharing parameters embedded in this party's key share.
    pub fn sharing_params(&self) -> (usize, usize) {
        (self.key_share.n, self.key_share.t)
    }

    /// Borrow this party's underlying RSS key share.
    pub(crate) fn key_share(&self) -> &PartyKeyShare {
        &self.key_share
    }

    /// Derive this party's additive signing share for a specific active signer set.
    ///
    /// Each RSS subset-piece is assigned to exactly one active signer (the smallest
    /// active party id in that subset). This avoids replicated-piece double counting.
    fn derive_signing_share(&self, active_party_ids: &[usize]) -> Result<PolyVecL, Error> {
        if !active_party_ids.contains(&self.id) {
            return Err(Error::InvalidParameters);
        }

        let mut active = [false; MAX_PARTIES];
        for &id in active_party_ids {
            if id >= MAX_PARTIES || active[id] {
                return Err(Error::InvalidParameters);
            }
            active[id] = true;
        }

        let subsets = enumerate_subsets(self.key_share.n, self.key_share.t);
        let mut share = PolyVecL::zero();

        for (local_idx, &global_subset_idx) in self.key_share.subset_indices.iter().enumerate() {
            let subset = subsets.get(global_subset_idx).ok_or(Error::InvalidShare)?;

            // Deterministically assign each subset-piece to one active signer.
            let owner = subset
                .iter()
                .filter(|&&member| member < MAX_PARTIES && active[member])
                .copied()
                .min()
                .ok_or(Error::InvalidParameters)?;

            if owner == self.id {
                share.add_assign(&self.key_share.s1_shares[local_idx]);
            }
        }

        Ok(share)
    }

    /// **Round 0 — Pre-commit**: Compute and return a binding hash H(w_i ‖ party_id).
    ///
    /// This hash must be broadcast before revealing the full commitment,
    /// preventing the coordinator from grinding challenges (ADV-1).
    /// The party internally runs the full commit logic and stores state.
    pub fn precommit<R: RngCore + CryptoRng>(
        &mut self,
        rng: &mut R,
        active_party_ids: &[usize],
        session_context: &[u8],
    ) -> Result<CommitmentHash, Error> {
        let commitment = self.commit_internal(rng, active_party_ids, session_context)?;
        let hash = commitment.binding_hash;
        let party_id = commitment.party_id;
        // Store w and keep y for sign()
        self.w = Some(commitment.w.clone());
        Ok(CommitmentHash { hash, party_id })
    }

    /// **Round 1 — Reveal**: Return the full commitment (after pre-commit).
    ///
    /// The coordinator verifies this against the pre-committed hash.
    pub fn reveal(&self) -> Option<Commitment> {
        self.w.as_ref().map(|w| {
            let binding_hash = compute_commitment_hash(w, self.id);
            Commitment {
                w: w.clone(),
                binding_hash,
                party_id: self.id,
            }
        })
    }

    /// **Round 1 — Commit**: Sample masking randomness and compute partial commitment.
    ///
    /// Uses **hedged nonce derivation** (ADV-4): the masking seed is derived from
    /// both RNG entropy AND the party's secret key material, so that even with a
    /// compromised RNG, different parties produce different nonces.
    ///
    /// Returns the partial commitment w_i with a binding hash.
    pub fn commit<R: RngCore + CryptoRng>(
        &mut self,
        rng: &mut R,
        active_party_ids: &[usize],
        session_context: &[u8],
    ) -> Result<Commitment, Error> {
        let commitment = self.commit_internal(rng, active_party_ids, session_context)?;
        self.w = Some(commitment.w.clone());
        Ok(commitment)
    }

    /// Internal commit logic shared by commit() and precommit().
    fn commit_internal<R: RngCore + CryptoRng>(
        &mut self,
        rng: &mut R,
        active_party_ids: &[usize],
        session_context: &[u8],
    ) -> Result<Commitment, Error> {
        // Build this round's additive signing share for the active signer set.
        let signing_s1_share = self.derive_signing_share(active_party_ids)?;

        // ── ADV-4: Hedged nonce derivation ──
        // seed = SHAKE-256(rng_entropy ‖ signing_share ‖ session_context ‖ counter)
        // This ensures that even if the RNG is compromised (e.g., VM snapshot),
        // the masking vector is still bound to:
        //   - this party's signing share
        //   - this signing session context
        //   - this party's per-round counter
        let mut rng_entropy = [0u8; CRHBYTES];
        rng.fill_bytes(&mut rng_entropy);

        let mut hasher = Shake256::default();
        hasher.update(&rng_entropy);
        // Mix in the party's secret key material for hedging
        for l in 0..L {
            for c in 0..N {
                hasher.update(&signing_s1_share.polys[l].coeffs[c].to_le_bytes());
            }
        }
        hasher.update(&(session_context.len() as u64).to_le_bytes());
        hasher.update(session_context);
        hasher.update(&self.id.to_le_bytes());
        hasher.update(&self.nonce_counter.to_le_bytes());

        let mut seed = [0u8; CRHBYTES];
        hasher.finalize_xof().read(&mut seed);
        self.nonce_counter = self.nonce_counter.wrapping_add(1);

        // SECURITY: zeroize raw entropy immediately
        rng_entropy.zeroize();

        // Sample masking vector y_i ∈ S_{γ₁}^L, then scale by active signer count.
        //
        // This keeps the aggregate y = Σ y_i within the FIPS-compatible γ₁ envelope.
        // Without this scaling, summing full-range party masks would overflow
        // the final z infinity bound almost always.
        let mut y = PolyVecL::zero();
        let scale = active_party_ids.len() as i32;
        for j in 0..L {
            Poly::uniform_gamma1(&mut y.polys[j], &seed, (self.id * L + j) as u16);
            for c in 0..N {
                y.polys[j].coeffs[c] /= scale;
            }
        }

        // SECURITY: seed must be zeroized to prevent masking vector recovery
        seed.zeroize();

        // Compute w_i = A · NTT(y_i)
        let mut y_hat = y.clone();
        y_hat.ntt();
        let mut w = self.a_hat.mul_vec(&y_hat);

        // SECURITY: y_hat contains NTT-domain secret masking vector
        y_hat.zeroize();

        w.reduce();
        w.invntt_tomont();
        w.reduce();
        w.caddq();

        // Store y for Round 3
        self.y = Some(y);
        self.signing_s1_share = Some(signing_s1_share);

        // Compute binding hash (ADV-1)
        let binding_hash = compute_commitment_hash(&w, self.id);

        Ok(Commitment {
            w,
            binding_hash,
            party_id: self.id,
        })
    }

    /// **Round 3 — Sign**: Compute partial response with hyperball rejection.
    ///
    /// Given the challenge polynomial c (derived from the aggregate commitment
    /// and message in Round 2), each party:
    ///
    /// 1. Computes z_i = c · s₁_i + y_i  (in NTT domain, then inverse)
    /// 2. **Hyperball check**: if ‖z_i‖₂² > B² → reject (return `LocalRejectionAbort`)
    /// 3. Returns z_i with session binding (ADV-6)
    ///
    /// This is the key innovation from ePrint 2026/013: rejection is performed
    /// **locally** using L₂-norm hyperballs instead of L∞-norm hypercubes,
    /// preventing exponential degradation of acceptance probability.
    pub fn sign(&mut self, challenge_hash: &[u8]) -> Result<PartialSignature, Error> {
        let y = self.y.take().ok_or(Error::InvalidShare)?;
        let mut signing_s1_share = self.signing_s1_share.take().ok_or(Error::InvalidShare)?;

        // Decode challenge polynomial c from the hash
        let mut c_poly = Poly::zero();
        Poly::challenge(&mut c_poly, challenge_hash);

        // Compute z_i = c · s₁_i + y_i
        // We need to do this in NTT domain: NTT(c) · NTT(s₁_i), then INTT
        let mut c_hat = c_poly.clone();
        c_hat.ntt();

        let mut z = PolyVecL::zero();
        let mut s1_hat = signing_s1_share.clone();
        s1_hat.ntt();

        for j in 0..L {
            // cs₁[j] = INTT(NTT(c) · NTT(s₁[j]))
            let mut cs1_j = Poly::zero();
            Poly::pointwise_montgomery(&mut cs1_j, &c_hat, &s1_hat.polys[j]);
            cs1_j.invntt_tomont();
            cs1_j.reduce();
            cs1_j.caddq();

            // z[j] = cs₁[j] + y[j]
            Poly::add(&mut z.polys[j], &cs1_j, &y.polys[j]);
            z.polys[j].reduce();
        }

        // ── Hyperball rejection check (L₂ norm) ──
        // Per ePrint 2026/013 §5.2: accept iff ‖z_i‖₂² ≤ B²
        // This is the ONLY per-party rejection check in the Mithril scheme.
        //
        // Note: We do NOT check ‖z_i‖∞ < γ₁ − β here because that bound
        // applies to the AGGREGATED z = Σ z_i, not individual partial responses.
        // Each party's z_i may exceed γ₁ − β because the RSS share s₁_i has
        // norm C(N-1,T-1)·η > η, making c·s₁_i larger than c·s₁ per-party.
        // The L∞ check is enforced at aggregation time instead.
        let l2_sq = z.l2_norm_squared();

        // SECURITY: zeroize NTT-domain secret material before returning
        c_hat.zeroize();
        s1_hat.zeroize();
        signing_s1_share.zeroize();

        if l2_sq > HYPERBALL_BOUND_SQ {
            return Err(Error::LocalRejectionAbort);
        }

        // Compute session binding (ADV-6): H(c̃ ‖ party_id)
        let session_binding = compute_session_binding(challenge_hash, self.id);

        Ok(PartialSignature {
            party_id: self.id,
            z,
            session_binding,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rss;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_party_commit() {
        let mut rng = StdRng::seed_from_u64(99);

        // Create a trivial secret
        let s1 = PolyVecL::zero();
        let s2 = PolyVecK::zero();

        let shares = rss::distribute_key(&s1, &s2, 3, 2, &mut rng).unwrap();

        let rho = [0u8; SEEDBYTES];
        let a_hat = MatrixA::expand(&rho);
        let mut party = Party::new(&shares[0], a_hat);

        let active = [0usize, 1usize];
        let commitment = party.commit(&mut rng, &active, b"test-session").unwrap();
        // Just verify it doesn't panic and produces non-trivial output
        assert_eq!(commitment.w.polys.len(), K);
        assert_eq!(commitment.party_id, 0);
        assert_ne!(commitment.binding_hash, [0u8; 32]);
    }

    #[test]
    fn test_commitment_binding_hash() {
        let mut rng = StdRng::seed_from_u64(101);
        let s1 = PolyVecL::zero();
        let s2 = PolyVecK::zero();
        let shares = rss::distribute_key(&s1, &s2, 3, 2, &mut rng).unwrap();
        let rho = [0u8; SEEDBYTES];
        let a_hat = MatrixA::expand(&rho);

        let mut party = Party::new(&shares[0], a_hat);
        let active = [0usize, 1usize];
        let commitment = party.commit(&mut rng, &active, b"test-session").unwrap();

        // Verify the binding hash matches
        let hash = CommitmentHash {
            hash: commitment.binding_hash,
            party_id: commitment.party_id,
        };
        assert!(verify_commitment(&commitment, &hash));

        // Verify wrong hash fails
        let wrong_hash = CommitmentHash {
            hash: [0u8; 32],
            party_id: commitment.party_id,
        };
        assert!(!verify_commitment(&commitment, &wrong_hash));
    }

    #[test]
    fn test_session_binding() {
        let mut rng = StdRng::seed_from_u64(102);
        let s1 = PolyVecL::zero();
        let s2 = PolyVecK::zero();
        let shares = rss::distribute_key(&s1, &s2, 3, 2, &mut rng).unwrap();
        let rho = [0u8; SEEDBYTES];
        let a_hat = MatrixA::expand(&rho);

        let mut party = Party::new(&shares[0], a_hat);
        let active = [0usize, 1usize];
        party.commit(&mut rng, &active, b"test-session").unwrap();
        let c_tilde = [42u8; CTILDEBYTES];

        if let Ok(partial) = party.sign(&c_tilde) {
            assert!(verify_session_binding(&partial, &c_tilde));
            // Wrong c_tilde should fail
            assert!(!verify_session_binding(&partial, &[0u8; CTILDEBYTES]));
        }
    }

    #[test]
    fn test_precommit_reveal_flow() {
        let mut rng = StdRng::seed_from_u64(103);
        let s1 = PolyVecL::zero();
        let s2 = PolyVecK::zero();
        let shares = rss::distribute_key(&s1, &s2, 3, 2, &mut rng).unwrap();
        let rho = [0u8; SEEDBYTES];
        let a_hat = MatrixA::expand(&rho);

        let mut party = Party::new(&shares[0], a_hat);

        // Pre-commit: get binding hash
        let active = [0usize, 1usize];
        let hash = party.precommit(&mut rng, &active, b"test-session").unwrap();

        // Reveal: get full commitment
        let commitment = party.reveal().expect("should have commitment");

        // Verify binding
        assert!(verify_commitment(&commitment, &hash));
    }
}
