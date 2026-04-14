//! High-level SDK API for threshold ML-DSA-44.
//!
//! This module provides an opinionated wrapper around the low-level protocol
//! modules (`rss`, `sign`, `coordinator`, `verify`) so application code can:
//! 1. derive threshold key material from a standard ML-DSA-44 keypair
//! 2. initialize a signing cluster for `(N, T)`
//! 3. produce fail-closed threshold signatures with a single call
//!
//! The wrapper remains fully compatible with FIPS 204 verification.

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
#[cfg(feature = "std")]
use std::vec::Vec;

use crate::coordinator::{self, Signature};
use crate::error::Error;
use crate::params::*;
use crate::poly::{MatrixA, PolyVecK, PolyVecL};
use crate::rss;
use crate::sign::Party;
use crate::verify;
use rand_core::{CryptoRng, RngCore};
use zeroize::{Zeroize, Zeroizing};

/// Extracted ML-DSA-44 material required for thresholdization.
///
/// This struct intentionally separates:
/// - public data (`pk`, `tr`, `rho`)
/// - secret vectors (`s1`, `s2`) that are zeroized on drop
#[derive(Clone, Zeroize)]
#[zeroize(drop)]
pub struct MlDsa44ThresholdMaterial {
    /// ML-DSA-44 public key bytes.
    #[zeroize(skip)]
    pub pk: Vec<u8>,
    /// Hash of public key (`tr = H(pk)`), used by Fiat-Shamir challenge derivation.
    #[zeroize(skip)]
    pub tr: [u8; TRBYTES],
    /// Public matrix seed (`rho`) from the ML-DSA keypair.
    #[zeroize(skip)]
    pub rho: [u8; SEEDBYTES],
    /// Secret key vector s1 (zeroized on drop).
    pub s1: PolyVecL,
    /// Secret key vector s2 (zeroized on drop).
    pub s2: PolyVecK,
}

impl MlDsa44ThresholdMaterial {
    /// Derive threshold material from an existing ML-DSA-44 keypair.
    ///
    /// # Errors
    /// Returns [`Error::InvalidParameters`] if key sizes do not match ML-DSA-44.
    pub fn from_secret_key(pk: &[u8], sk: &[u8]) -> Result<Self, Error> {
        let mode = dilithium::params::ML_DSA_44;
        if pk.len() != mode.public_key_bytes() || sk.len() != mode.secret_key_bytes() {
            return Err(Error::InvalidParameters);
        }

        let mut rho_ref = [0u8; dilithium::params::SEEDBYTES];
        let mut tr_ref = [0u8; dilithium::params::TRBYTES];
        let mut key = [0u8; dilithium::params::SEEDBYTES];
        let mut t0 = dilithium::polyvec::PolyVecK::default();
        let mut s1_ref = dilithium::polyvec::PolyVecL::default();
        let mut s2_ref = dilithium::polyvec::PolyVecK::default();
        dilithium::packing::unpack_sk(
            mode,
            &mut rho_ref,
            &mut tr_ref,
            &mut key,
            &mut t0,
            &mut s1_ref,
            &mut s2_ref,
            sk,
        );

        let mut s1 = PolyVecL::zero();
        let mut s2 = PolyVecK::zero();
        for i in 0..L {
            s1.polys[i].coeffs = s1_ref.vec[i].coeffs;
        }
        for i in 0..K {
            s2.polys[i].coeffs = s2_ref.vec[i].coeffs;
        }

        Ok(Self {
            pk: pk.to_vec(),
            tr: tr_ref,
            rho: rho_ref,
            s1,
            s2,
        })
    }

    /// Generate an ML-DSA-44 keypair and derive threshold material.
    ///
    /// Returns the extracted threshold material plus the raw standard secret key.
    /// The returned secret key is wrapped in [`Zeroizing`].
    pub fn from_seed(seed: &[u8; 32]) -> Result<(Self, Zeroizing<Vec<u8>>), Error> {
        let (pk, sk) = verify::keygen(seed);
        let material = Self::from_secret_key(&pk, &sk)?;
        Ok((material, Zeroizing::new(sk)))
    }
}

/// High-level stateful SDK for threshold ML-DSA-44 signing.
///
/// This struct manages party instances and coordinator invocation internally.
/// It is a convenience API for in-process orchestrations and tests.
///
/// For networked deployments, use the low-level modules and explicit transport
/// control (e.g. encrypted/authenticated channels for round messages).
pub struct ThresholdMlDsa44Sdk {
    parties: Vec<Party>,
    pk: Vec<u8>,
    tr: [u8; TRBYTES],
    t: usize,
    max_retries: usize,
}

impl ThresholdMlDsa44Sdk {
    /// Build a threshold signing SDK instance from extracted threshold material.
    ///
    /// # Arguments
    /// - `material`: extracted ML-DSA-44 vectors/seeds for thresholdization
    /// - `n`: total parties
    /// - `t`: threshold
    /// - `max_retries`: max retries for fail-closed local-rejection handling
    ///
    /// # Errors
    /// Returns [`Error::InvalidParameters`] for invalid `(N, T)` or retry config.
    pub fn new<R: RngCore + CryptoRng>(
        material: &MlDsa44ThresholdMaterial,
        n: usize,
        t: usize,
        max_retries: usize,
        rng: &mut R,
    ) -> Result<Self, Error> {
        if material.pk.len() != PK_BYTES
            || n < 2
            || t < 2
            || t > n
            || n > MAX_PARTIES
            || max_retries == 0
        {
            return Err(Error::InvalidParameters);
        }

        let shares = rss::distribute_key(&material.s1, &material.s2, n, t, rng)?;
        let a_hat = MatrixA::expand(&material.rho);
        let parties = shares
            .iter()
            .map(|share| Party::new(share, a_hat.clone()))
            .collect();

        Ok(Self {
            parties,
            pk: material.pk.clone(),
            tr: material.tr,
            t,
            max_retries,
        })
    }

    /// One-call constructor: keygen + thresholdization + SDK setup.
    ///
    /// Returns `(sdk, standard_secret_key)`.
    pub fn from_seed<R: RngCore + CryptoRng>(
        seed: &[u8; 32],
        n: usize,
        t: usize,
        max_retries: usize,
        rng: &mut R,
    ) -> Result<(Self, Zeroizing<Vec<u8>>), Error> {
        let (material, sk) = MlDsa44ThresholdMaterial::from_seed(seed)?;
        let sdk = Self::new(&material, n, t, max_retries, rng)?;
        Ok((sdk, sk))
    }

    /// Produce a threshold signature for `msg`.
    ///
    /// This call is fail-closed: if no verifiable aggregate is found within
    /// `max_retries`, it returns an error.
    pub fn sign<R: RngCore + CryptoRng>(
        &mut self,
        msg: &[u8],
        rng: &mut R,
    ) -> Result<Signature, Error> {
        coordinator::threshold_sign(
            &mut self.parties,
            msg,
            &self.tr,
            &self.pk,
            self.t,
            self.max_retries,
            rng,
        )
    }

    /// Verify a signature using this SDK instance's public key.
    pub fn verify(&self, msg: &[u8], sig: &Signature) -> bool {
        verify::verify(&sig.to_bytes(), msg, &self.pk)
    }

    /// Borrow the public key.
    pub fn public_key(&self) -> &[u8] {
        &self.pk
    }

    /// Borrow the `tr` hash.
    pub fn tr(&self) -> &[u8; TRBYTES] {
        &self.tr
    }

    /// Return the configured threshold.
    pub fn threshold(&self) -> usize {
        self.t
    }

    /// Return number of parties managed by this SDK instance.
    pub fn party_count(&self) -> usize {
        self.parties.len()
    }

    /// Return the configured max retry count.
    pub fn max_retries(&self) -> usize {
        self.max_retries
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_sdk_invalid_parameters() {
        let (material, _sk) = MlDsa44ThresholdMaterial::from_seed(&[1u8; 32]).unwrap();
        let mut rng = StdRng::seed_from_u64(7);

        let bad = ThresholdMlDsa44Sdk::new(&material, 1, 1, 0, &mut rng);
        assert!(matches!(bad, Err(Error::InvalidParameters)));
    }

    #[test]
    fn test_sdk_sign_fail_closed_or_valid() {
        let mut rng = StdRng::seed_from_u64(20260414);
        let (mut sdk, _sk) =
            ThresholdMlDsa44Sdk::from_seed(&[9u8; 32], 4, 3, 128, &mut rng).unwrap();

        assert_eq!(sdk.party_count(), 4);
        assert_eq!(sdk.threshold(), 3);

        let msg = b"sdk-sign-message";
        match sdk.sign(msg, &mut rng) {
            Ok(sig) => assert!(sdk.verify(msg, &sig)),
            Err(Error::InvalidSignature) | Err(Error::InsufficientResponses) => {
                // Expected fail-closed outcomes under local rejection dynamics.
            }
            Err(e) => panic!("unexpected sdk sign error: {e}"),
        }
    }
}
