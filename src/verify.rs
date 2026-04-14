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

/// Verify an ML-DSA-44 signature against a public key and message.
///
/// This is a thin wrapper around `dilithium::sign::verify()` using
/// the ML-DSA-44 parameter set.
///
/// # Arguments
/// * `sig` — the encoded signature (SIG_BYTES bytes)
/// * `msg` — the signed message
/// * `pk` — the encoded public key (PK_BYTES bytes)
///
/// # Returns
/// `true` if the signature is valid, `false` otherwise.
pub fn verify(sig: &[u8], msg: &[u8], pk: &[u8]) -> bool {
    dilithium::sign::verify(dilithium::params::ML_DSA_44, sig, msg, &[], pk)
}

/// Verify an ML-DSA-65 signature.
pub fn verify_65(sig: &[u8], msg: &[u8], pk: &[u8]) -> bool {
    dilithium::sign::verify(dilithium::params::ML_DSA_65, sig, msg, &[], pk)
}

/// Verify an ML-DSA-87 signature.
pub fn verify_87(sig: &[u8], msg: &[u8], pk: &[u8]) -> bool {
    dilithium::sign::verify(dilithium::params::ML_DSA_87, sig, msg, &[], pk)
}

/// Generate an ML-DSA-44 key pair (for testing / key generation).
///
/// Returns (pk, sk) as raw byte vectors.
pub fn keygen(seed: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
    dilithium::sign::keypair(dilithium::params::ML_DSA_44, seed)
}

/// Sign a message using standard (non-threshold) ML-DSA-44.
///
/// This is used for reference comparison in tests.
pub fn sign_standard(msg: &[u8], sk: &[u8], rng_seed: &[u8; 32]) -> Vec<u8> {
    let mode = dilithium::params::ML_DSA_44;
    let mut sig = vec![0u8; mode.signature_bytes()];
    dilithium::sign::sign_signature(mode, &mut sig, msg, &[], rng_seed, sk);
    sig
}
