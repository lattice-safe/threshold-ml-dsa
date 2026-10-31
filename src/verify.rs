//! ML-DSA (FIPS 204) verification — delegated to the `dilithium-rs` crate.
//!
//! The whole point of threshold ML-DSA is that the aggregated threshold signature
//! is bit-for-bit compatible with standard FIPS 204 verification. Rather than
//! re-implementing ML-DSA verification, we delegate to the
//! [`dilithium-rs`](https://crates.io/crates/dilithium-rs) crate, our pure-Rust
//! FIPS 204 implementation which has been validated against all 100 NIST KAT
//! vectors for ML-DSA-44/65/87.
//!
//! This provides the strongest possible proof of FIPS 204 interoperability:
//! if the threshold signature passes `dilithium::sign::verify()`, it is
//! indistinguishable from a standard ML-DSA signature.

extern crate alloc;
use alloc::{vec, vec::Vec};

use crate::error::Error;
use crate::params::{ml_dsa_65, ml_dsa_87, PK_BYTES, SIG_BYTES, SK_BYTES};
use zeroize::Zeroizing;

/// Verify an ML-DSA-44 signature against a public key and message.
///
/// This is a thin wrapper around `dilithium::sign::verify()` using
/// the ML-DSA-44 parameter set.
///
/// # Arguments
/// * `sig` — the encoded signature (`SIG_BYTES` bytes)
/// * `msg` — the signed message
/// * `pk` — the encoded public key (`PK_BYTES` bytes)
///
/// # Returns
/// `true` if the signature is valid, `false` otherwise.
#[must_use] 
pub fn verify(sig: &[u8], msg: &[u8], pk: &[u8]) -> bool {
    // dilithium-rs indexes the public key without a length check.
    if pk.len() != PK_BYTES || sig.len() != SIG_BYTES {
        return false;
    }
    dilithium::sign::verify(dilithium::params::ML_DSA_44, sig, msg, &[], pk)
}

/// Verify an ML-DSA-65 signature.
#[must_use] 
pub fn verify_65(sig: &[u8], msg: &[u8], pk: &[u8]) -> bool {
    if pk.len() != ml_dsa_65::PK_BYTES || sig.len() != ml_dsa_65::SIG_BYTES {
        return false;
    }
    dilithium::sign::verify(dilithium::params::ML_DSA_65, sig, msg, &[], pk)
}

/// Verify an ML-DSA-87 signature.
#[must_use] 
pub fn verify_87(sig: &[u8], msg: &[u8], pk: &[u8]) -> bool {
    if pk.len() != ml_dsa_87::PK_BYTES || sig.len() != ml_dsa_87::SIG_BYTES {
        return false;
    }
    dilithium::sign::verify(dilithium::params::ML_DSA_87, sig, msg, &[], pk)
}

/// Generate an ML-DSA-44 key pair (for testing / key generation).
///
/// Returns (pk, sk) as raw byte vectors.
#[must_use] 
pub fn keygen(seed: &[u8; 32]) -> (Vec<u8>, Zeroizing<Vec<u8>>) {
    let (pk, sk) = dilithium::sign::keypair(dilithium::params::ML_DSA_44, seed);
    (pk, Zeroizing::new(sk))
}

/// Sign a message using standard (non-threshold) ML-DSA-44.
///
/// This is used for reference comparison in tests.
pub fn sign_standard(msg: &[u8], sk: &[u8], rng_seed: &[u8; 32]) -> Result<Vec<u8>, Error> {
    sign_with_context(dilithium::params::ML_DSA_44, SK_BYTES, msg, &[], sk, rng_seed)
}

/// Generate an ML-DSA-65 key pair.
#[must_use]
pub fn keygen_65(seed: &[u8; 32]) -> (Vec<u8>, Zeroizing<Vec<u8>>) {
    let (pk, sk) = dilithium::sign::keypair(dilithium::params::ML_DSA_65, seed);
    (pk, Zeroizing::new(sk))
}

/// Sign a message using standard ML-DSA-65.
pub fn sign_standard_65(msg: &[u8], sk: &[u8], rng_seed: &[u8; 32]) -> Result<Vec<u8>, Error> {
    sign_with_context(dilithium::params::ML_DSA_65, ml_dsa_65::SK_BYTES, msg, &[], sk, rng_seed)
}

/// Generate an ML-DSA-87 key pair.
#[must_use]
pub fn keygen_87(seed: &[u8; 32]) -> (Vec<u8>, Zeroizing<Vec<u8>>) {
    let (pk, sk) = dilithium::sign::keypair(dilithium::params::ML_DSA_87, seed);
    (pk, Zeroizing::new(sk))
}

/// Sign a message using standard ML-DSA-87.
pub fn sign_standard_87(msg: &[u8], sk: &[u8], rng_seed: &[u8; 32]) -> Result<Vec<u8>, Error> {
    sign_with_context(dilithium::params::ML_DSA_87, ml_dsa_87::SK_BYTES, msg, &[], sk, rng_seed)
}

/// FIPS 204 `ML-DSA.Sign` for `mode` with context string `ctx`.
fn sign_with_context(
    mode: dilithium::params::DilithiumMode,
    sk_len: usize,
    msg: &[u8],
    ctx: &[u8],
    sk: &[u8],
    rng_seed: &[u8; 32],
) -> Result<Vec<u8>, Error> {
    if sk.len() != sk_len {
        return Err(Error::InvalidParameters);
    }
    let mut sig = vec![0u8; mode.signature_bytes()];
    // Fails only for a context longer than 255 bytes.
    if dilithium::sign::sign_signature(mode, &mut sig, msg, ctx, rng_seed, sk) != 0 {
        return Err(Error::InvalidParameters);
    }
    Ok(sig)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ml_dsa_65_keygen_sign_verify() {
        let seed = [42u8; 32];
        let (pk, sk) = keygen_65(&seed);
        assert_eq!(pk.len(), ml_dsa_65::PK_BYTES);
        assert_eq!(sk.len(), ml_dsa_65::SK_BYTES);

        let msg = b"ML-DSA-65 verification test";
        let sig = sign_standard_65(msg, &sk, &[99u8; 32]).unwrap();
        assert_eq!(sig.len(), ml_dsa_65::SIG_BYTES);
        assert!(verify_65(&sig, msg, &pk));
    }

    #[test]
    fn test_ml_dsa_87_keygen_sign_verify() {
        let seed = [42u8; 32];
        let (pk, sk) = keygen_87(&seed);
        assert_eq!(pk.len(), ml_dsa_87::PK_BYTES);
        assert_eq!(sk.len(), ml_dsa_87::SK_BYTES);

        let msg = b"ML-DSA-87 verification test";
        let sig = sign_standard_87(msg, &sk, &[99u8; 32]).unwrap();
        assert_eq!(sig.len(), ml_dsa_87::SIG_BYTES);
        assert!(verify_87(&sig, msg, &pk));
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;

    #[test]
    fn overlong_context_is_rejected() {
        let (_, sk) = keygen(&[1u8; 32]);
        let ctx = [0u8; 256];
        assert_eq!(
            sign_with_context(dilithium::params::ML_DSA_44, SK_BYTES, b"m", &ctx, &sk, &[0u8; 32]).err(),
            Some(Error::InvalidParameters)
        );
        assert!(sign_with_context(dilithium::params::ML_DSA_44, SK_BYTES, b"m", &ctx[..255], &sk, &[0u8; 32]).is_ok());
    }
}

#[cfg(test)]
mod wrapper_tests {
    use super::*;

    #[test]
    fn standard_signing_and_length_checks() {
        let (pk, sk) = keygen(&[20u8; 32]);
        let sig = sign_standard(b"m", &sk, &[0u8; 32]).unwrap();
        assert!(verify(&sig, b"m", &pk));
        assert!(!verify(&sig, b"m", &pk[1..]));
        assert!(!verify(&sig[1..], b"m", &pk));
        assert!(sign_standard(b"m", &sk[1..], &[0u8; 32]).is_err());

        let (pk65, sk65) = keygen_65(&[21u8; 32]);
        let sig65 = sign_standard_65(b"m", &sk65, &[0u8; 32]).unwrap();
        assert!(verify_65(&sig65, b"m", &pk65));
        assert!(!verify_65(&sig65, b"m", &pk65[1..]));

        let (pk87, sk87) = keygen_87(&[22u8; 32]);
        let sig87 = sign_standard_87(b"m", &sk87, &[0u8; 32]).unwrap();
        assert!(verify_87(&sig87, b"m", &pk87));
        assert!(!verify_87(&sig87, b"m", &pk87[1..]));
    }
}
