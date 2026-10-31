//! Error types for the threshold ML-DSA protocol.
//!
//! Uses manual `Display` impls instead of `thiserror` derive macros
//! to maintain `#![no_std]` compatibility.

use core::fmt;

/// Errors that can occur during threshold ML-DSA operations.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// A party's partial response `z_i` exceeded the hyperball L₂ bound.
    /// This is the expected "soft failure" in the Mithril protocol —
    /// the coordinator simply retries the round.
    LocalRejectionAbort,

    /// An RSS share failed validation (e.g., wrong subset, corrupted data).
    InvalidShare,

    /// Fewer than T parties provided valid responses after rejection sampling.
    InsufficientResponses,

    /// The aggregated or individual signature failed FIPS 204 verification.
    InvalidSignature,

    /// The (N, T) threshold configuration is invalid (e.g., T > N or N > `MAX_PARTIES`).
    InvalidParameters,

    /// The low-bits hint check failed (‖w − cs₂‖∞ ≥ γ₂ − β).
    HintCheckFailed,

    /// A [`crate::nonce::NonceGuard`] refused a Round 1 session ID that was
    /// already used with this key (or its counter is exhausted).
    SessionReplay,

    /// A DKG message is malformed, from an unexpected sender, or for another
    /// session.
    DkgMalformedMessage,

    /// A party's revealed `t_I` does not match its commitment.
    DkgCommitmentMismatch {
        /// The party whose reveal failed
        party: u8,
    },

    /// A received subset seed does not match the dealer's commitment. Either
    /// the dealer or the recipient is dishonest; abort the session.
    DkgSeedMismatch {
        /// Claimed dealer
        dealer: u8,
        /// Subset mask
        mask: u8,
    },

    /// Parties saw different DKG transcripts (equivocation).
    DkgInconsistentTranscript,
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
            Error::SessionReplay => write!(f, "round 1 session ID already used with this key"),
            Error::DkgMalformedMessage => write!(f, "malformed DKG message"),
            Error::DkgCommitmentMismatch { party } => {
                write!(f, "DKG reveal from party {party} does not match its commitment")
            }
            Error::DkgSeedMismatch { dealer, mask } => {
                write!(f, "DKG seed for subset {mask:#04x} from dealer {dealer} does not match its commitment")
            }
            Error::DkgInconsistentTranscript => write!(f, "DKG parties saw different transcripts"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}
