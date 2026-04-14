//! Error types for the threshold ML-DSA protocol.
//!
//! Uses manual `Display` impls instead of `thiserror` derive macros
//! to maintain `#![no_std]` compatibility.

use core::fmt;

/// Errors that can occur during threshold ML-DSA operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// A party's partial response z_i exceeded the hyperball L₂ bound.
    /// This is the expected "soft failure" in the Mithril protocol —
    /// the coordinator simply retries the round.
    LocalRejectionAbort,

    /// An RSS share failed validation (e.g., wrong subset, corrupted data).
    InvalidShare,

    /// Fewer than T parties provided valid responses after rejection sampling.
    InsufficientResponses,

    /// The aggregated or individual signature failed FIPS 204 verification.
    InvalidSignature,

    /// The (N, T) threshold configuration is invalid (e.g., T > N or N > MAX_PARTIES).
    InvalidParameters,

    /// The low-bits hint check failed (‖w − cs₂‖∞ ≥ γ₂ − β).
    HintCheckFailed,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::LocalRejectionAbort => {
                write!(f, "local hyperball rejection: ‖z_i‖₂² exceeded bound")
            }
            Error::InvalidShare => write!(f, "invalid RSS share"),
            Error::InsufficientResponses => {
                write!(f, "insufficient valid responses from threshold parties")
            }
            Error::InvalidSignature => write!(f, "signature verification failed"),
            Error::InvalidParameters => write!(f, "invalid threshold parameters (N, T)"),
            Error::HintCheckFailed => write!(f, "low-bits hint check failed"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}
