//! Nonce-reuse protection for Round 1.
//!
//! If a party's RNG ever repeats output (VM snapshot restore, process fork,
//! fixed seed), the hedged nonce derivation alone does not help: the
//! coordinator controls `session_id`, the active set and the message, so it
//! can replay them, obtain the same nonce `y` twice, answer it under two
//! different challenges and solve `z − z′ = (c − c′)·s₁_I` for the partial
//! secret.
//!
//! A [`NonceGuard`] closes this gap as long as its state is not rolled back
//! together with the RNG: [`crate::sign::round1`] reserves the session ID
//! before drawing randomness, refuses a session ID it has already seen, and
//! mixes a strictly increasing per-key counter into the nonce seed. Two
//! Round 1 calls therefore never derive the same nonce, even if the RNG
//! returns identical bytes.
//!
//! What remains unprotected:
//! - rolling back the RNG and the guard's storage together (VM snapshot,
//!   restoring a backup);
//! - running the same party key on two replicas with separate guards;
//! - losing the guard's state (counter restarts at its stored high-water
//!   mark) at the same time as an RNG repeat.

#[cfg(not(feature = "std"))]
use alloc::collections::{BTreeMap, BTreeSet, VecDeque};
#[cfg(feature = "std")]
use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::error::Error;

/// Records Round 1 sessions so a nonce is never derived twice.
pub trait NonceGuard {
    /// Reserve `session_id` for the key identified by `key_id`.
    ///
    /// Implementations MUST:
    /// - return `Err(Error::SessionReplay)` if `session_id` was already
    ///   reserved for `key_id` (within the retained history);
    /// - return strictly increasing counters per `key_id` over the life of
    ///   the persisted state;
    /// - make the reservation durable before returning `Ok` if the guard is
    ///   meant to survive restarts.
    fn reserve(&mut self, key_id: &[u8; 32], session_id: &[u8; 32]) -> Result<u64, Error>;
}

impl<G: NonceGuard + ?Sized> NonceGuard for &mut G {
    fn reserve(&mut self, key_id: &[u8; 32], session_id: &[u8; 32]) -> Result<u64, Error> {
        (**self).reserve(key_id, session_id)
    }
}

/// In-memory [`NonceGuard`]: per-key counters plus a bounded history of
/// session IDs.
///
/// **Not durable.** It protects a single long-lived process. To survive
/// restarts, persist [`InMemoryNonceGuard::counter`] and restore it with
/// [`InMemoryNonceGuard::with_start_counter`], or implement [`NonceGuard`]
/// over durable storage.
pub struct InMemoryNonceGuard {
    capacity: usize,
    keys: BTreeMap<[u8; 32], KeyHistory>,
}

#[derive(Default)]
struct KeyHistory {
    next: u64,
    seen: BTreeSet<[u8; 32]>,
    order: VecDeque<[u8; 32]>,
}

impl InMemoryNonceGuard {
    /// Default number of session IDs remembered per key.
    pub const DEFAULT_CAPACITY: usize = 4096;

    /// Create a guard remembering [`Self::DEFAULT_CAPACITY`] sessions per key.
    #[must_use]
    pub fn new() -> Self {
        Self::with_capacity(Self::DEFAULT_CAPACITY)
    }

    /// Create a guard remembering up to `capacity` session IDs per key.
    ///
    /// Older IDs are evicted first; the counter still guarantees distinct
    /// nonce seeds after eviction.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            keys: BTreeMap::new(),
        }
    }

    /// Restore a persisted high-water mark: the next counter handed out for
    /// `key_id` will be at least `next`.
    #[must_use]
    pub fn with_start_counter(mut self, key_id: [u8; 32], next: u64) -> Self {
        let h = self.keys.entry(key_id).or_default();
        h.next = h.next.max(next);
        self
    }

    /// The next counter value for `key_id` (persist this to resume later).
    #[must_use]
    pub fn counter(&self, key_id: &[u8; 32]) -> u64 {
        self.keys.get(key_id).map_or(0, |h| h.next)
    }
}

impl Default for InMemoryNonceGuard {
    fn default() -> Self {
        Self::new()
    }
}

impl NonceGuard for InMemoryNonceGuard {
    fn reserve(&mut self, key_id: &[u8; 32], session_id: &[u8; 32]) -> Result<u64, Error> {
        let h = self.keys.entry(*key_id).or_default();
        if h.seen.contains(session_id) {
            return Err(Error::SessionReplay);
        }
        let ctr = h.next;
        h.next = ctr.checked_add(1).ok_or(Error::SessionReplay)?;
        if h.order.len() >= self.capacity {
            if let Some(old) = h.order.pop_front() {
                h.seen.remove(&old);
            }
        }
        h.seen.insert(*session_id);
        h.order.push_back(*session_id);
        Ok(ctr)
    }
}
