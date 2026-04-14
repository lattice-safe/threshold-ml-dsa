//! # Threshold ML-DSA
//!
//! A production-grade implementation of the "Efficient Threshold ML-DSA" scheme
//! (ePrint 2026/013, Celi et al. — the "Mithril" scheme).
//!
//! This crate implements a hardened 4-round distributed signing protocol for ML-DSA (FIPS 204)
//! using Replicated Secret Sharing (RSS) and hyperball-based local rejection sampling.
//! The resulting threshold signatures are bit-for-bit verifiable by any standard,
//! unmodified ML-DSA verifier.
//!
//! ## Key Innovations (from ePrint 2026/013)
//!
//! 1. **Short Replicated Secret Sharing** — avoids Lagrange multiplier blow-up
//! 2. **Hyperball Local Rejection** — L₂-norm rejection per party, no global aborts
//! 3. **Fail-Closed Aggregation** — returned signatures are end-to-end verified
//! 4. **SDK Layer** — high-level API in [`sdk`] for easier integration
//!
//! ## Usage
//!
//! ```rust,no_run
//! use threshold_ml_dsa::{coordinator, params, rss, sdk};
//! ```

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(not(feature = "std"))]
extern crate alloc;

pub mod coordinator;
pub mod error;
pub mod fvec;
pub mod params;
pub mod partition;
pub mod poly;
pub mod rss;
pub mod sdk;
pub mod sign;
pub mod verify;
