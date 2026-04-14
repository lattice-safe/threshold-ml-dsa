//! Replicated Secret Sharing (RSS) for threshold ML-DSA key distribution.
//!
//! ## Overview
//!
//! Standard Shamir secret sharing introduces large Lagrange interpolation
//! multipliers that blow up lattice coefficient sizes, making the resulting
//! shares incompatible with the short-vector requirements of ML-DSA.
//!
//! Replicated Secret Sharing (RSS) avoids this entirely:
//! - The secret is split into additive shares indexed by all subsets of size `N-T+1`.
//! - Each party receives all shares for subsets that contain it.
//! - Any T qualifying parties can reconstruct the secret by summing exactly one
//!   copy of each indexed subset share.
//! - Crucially, **no Lagrange multipliers** are needed — the shares are short
//!   by construction.
//!
//! ## Threshold Semantics
//!
//! For (N, T)-threshold:
//! - N = total number of parties
//! - T = minimum number of parties needed to sign
//! - Share-index subsets = all (N-T+1)-subsets of {0, ..., N-1}
//! - Each such subset gets one additive share piece
//! - A party i holds share pieces for all those subsets containing i

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use crate::error::Error;
use crate::params::*;
use crate::poly::{PolyVecK, PolyVecL};
use rand_core::{CryptoRng, RngCore};
use zeroize::Zeroize;

/// A single party's key share: the collection of RSS share pieces
/// for all (N-T+1)-subsets containing this party.
#[derive(Clone, Zeroize)]
#[zeroize(drop)]
pub struct PartyKeyShare {
    /// Party index in [0, N).
    #[zeroize(skip)]
    pub party_id: usize,

    /// Total number of parties in this sharing instance.
    #[zeroize(skip)]
    pub n: usize,

    /// Threshold for this sharing instance.
    #[zeroize(skip)]
    pub t: usize,

    /// Share pieces for s₁ (the L-vector part of the secret key).
    /// One PolyVecL per (N-T+1)-subset containing this party.
    pub s1_shares: Vec<PolyVecL>,

    /// Share pieces for s₂ (the K-vector part of the secret key).
    /// One PolyVecK per (N-T+1)-subset containing this party.
    pub s2_shares: Vec<PolyVecK>,

    /// The indices of the (N-T+1)-subsets that this party participates in.
    /// Each entry is an index into the global subset list.
    #[zeroize(skip)]
    pub subset_indices: Vec<usize>,
}

impl core::fmt::Debug for PartyKeyShare {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "PartyKeyShare {{ party_id: {}, n: {}, t: {}, shares: {} }}",
            self.party_id,
            self.n,
            self.t,
            self.s1_shares.len()
        )
    }
}

impl PartyKeyShare {
    /// Compute this party's "effective" s₁ share by summing all share pieces.
    /// This is what the party uses in the signing protocol.
    pub fn effective_s1(&self) -> PolyVecL {
        let mut result = PolyVecL::zero();
        for share in &self.s1_shares {
            result.add_assign(share);
        }
        result
    }

    /// Compute this party's "effective" s₂ share by summing all share pieces.
    pub fn effective_s2(&self) -> PolyVecK {
        let mut result = PolyVecK::zero();
        for share in &self.s2_shares {
            result.add_assign(share);
        }
        result
    }
}

/// Enumerate all (N-T+1)-subsets of {0, 1, ..., N-1} used as RSS piece indices.
///
/// Returns a list of subsets, where each subset is a sorted Vec of party indices.
/// For example:
/// - N=3, T=2 → M=2: [{0,1}, {0,2}, {1,2}]
/// - N=4, T=2 → M=3: [{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}]
pub fn enumerate_subsets(n: usize, t: usize) -> Vec<Vec<usize>> {
    if n == 0 || t == 0 || t > n {
        return Vec::new();
    }
    let subset_size = n - t + 1;
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(subset_size);
    enumerate_subsets_helper(n, subset_size, 0, &mut current, &mut result);
    result
}

fn enumerate_subsets_helper(
    n: usize,
    t: usize,
    start: usize,
    current: &mut Vec<usize>,
    result: &mut Vec<Vec<usize>>,
) {
    if current.len() == t {
        result.push(current.clone());
        return;
    }
    let remaining = t - current.len();
    for i in start..=(n - remaining) {
        current.push(i);
        enumerate_subsets_helper(n, t, i + 1, current, result);
        current.pop();
    }
}

/// Distribute an ML-DSA secret key (s₁, s₂) among N parties using
/// Replicated Secret Sharing with threshold T.
///
/// # Arguments
/// * `s1` — the secret key vector s₁ ∈ ℤ_q^L
/// * `s2` — the secret key vector s₂ ∈ ℤ_q^K
/// * `n` — total number of parties
/// * `t` — threshold (minimum parties to sign)
/// * `rng` — cryptographically secure RNG
///
/// # Returns
/// A vector of N `PartyKeyShare`s, one per party. Each share is zeroized on drop.
///
/// # Errors
/// Returns `Error::InvalidParameters` if N > MAX_PARTIES, T > N, or T < 2.
pub fn distribute_key<R: RngCore + CryptoRng>(
    s1: &PolyVecL,
    s2: &PolyVecK,
    n: usize,
    t: usize,
    rng: &mut R,
) -> Result<Vec<PartyKeyShare>, Error> {
    // Validate parameters
    if n > MAX_PARTIES || t > n || t < 2 || n < 2 {
        return Err(Error::InvalidParameters);
    }

    let subsets = enumerate_subsets(n, t);
    let num_subsets = subsets.len();

    // Generate random additive shares for each subset.
    // For the last subset, the share is computed as: secret - (sum of all other shares).
    // This ensures that the sum of all shares equals the secret.
    let mut s1_pieces: Vec<PolyVecL> = Vec::with_capacity(num_subsets);
    let mut s2_pieces: Vec<PolyVecK> = Vec::with_capacity(num_subsets);

    // Running sum of shares (to compute the last one)
    let mut s1_sum = PolyVecL::zero();
    let mut s2_sum = PolyVecK::zero();

    for idx in 0..num_subsets {
        if idx < num_subsets - 1 {
            // Sample SHORT random share pieces, bounded by η.
            // This is the critical RSS property: shares stay short so each party's
            // effective contribution remains bounded after aggregating the pieces it holds.
            let mut s1_share = PolyVecL::zero();
            let mut s2_share = PolyVecK::zero();

            for l in 0..L {
                for c in 0..N {
                    let mut buf = [0u8; 4];
                    loop {
                        rng.fill_bytes(&mut buf);
                        let raw = u32::from_le_bytes(buf);
                        // Rejection sampling: reject values that cause modular bias
                        // Accept only if raw < largest multiple of (2η+1) ≤ 2^32
                        let range = 2 * ETA + 1; // 5 for η=2
                        let limit = u32::MAX - (u32::MAX % range);
                        if raw < limit {
                            let val = (raw % range) as i32 - ETA as i32;
                            s1_share.polys[l].coeffs[c] = val;
                            break;
                        }
                    }
                }
            }
            for k in 0..K {
                for c in 0..N {
                    let mut buf = [0u8; 4];
                    loop {
                        rng.fill_bytes(&mut buf);
                        let raw = u32::from_le_bytes(buf);
                        let range = 2 * ETA + 1;
                        let limit = u32::MAX - (u32::MAX % range);
                        if raw < limit {
                            let val = (raw % range) as i32 - ETA as i32;
                            s2_share.polys[k].coeffs[c] = val;
                            break;
                        }
                    }
                }
            }

            s1_sum.add_assign(&s1_share);
            s2_sum.add_assign(&s2_share);

            s1_pieces.push(s1_share);
            s2_pieces.push(s2_share);
        } else {
            // Last share: secret - sum_of_others
            // Since original secret ∈ [-η, η] and each random piece ∈ [-η, η],
            // the last piece = s - Σ(random) remains bounded by the number of
            // indexed subsets times η, so it is still short for our small-N setting.
            let mut s1_last = PolyVecL::zero();
            let mut s2_last = PolyVecK::zero();
            PolyVecL::sub(&mut s1_last, s1, &s1_sum);
            PolyVecK::sub(&mut s2_last, s2, &s2_sum);
            s1_pieces.push(s1_last);
            s2_pieces.push(s2_last);
        }
    }

    // SECURITY: zeroize running sums (contain secret-correlated data)
    s1_sum.zeroize();
    s2_sum.zeroize();

    // Now distribute share pieces to parties:
    // Party i receives all pieces for subsets containing i.
    let mut party_shares: Vec<PartyKeyShare> = Vec::with_capacity(n);
    for party_id in 0..n {
        let mut s1_shares = Vec::new();
        let mut s2_shares = Vec::new();
        let mut subset_indices = Vec::new();

        for (subset_idx, subset) in subsets.iter().enumerate() {
            if subset.contains(&party_id) {
                s1_shares.push(s1_pieces[subset_idx].clone());
                s2_shares.push(s2_pieces[subset_idx].clone());
                subset_indices.push(subset_idx);
            }
        }

        party_shares.push(PartyKeyShare {
            party_id,
            n,
            t,
            s1_shares,
            s2_shares,
            subset_indices,
        });
    }

    // SECURITY: zeroize intermediate share piece vectors (M1)
    // These contain copies of all RSS share pieces — secret material.
    for piece in s1_pieces.iter_mut() {
        piece.zeroize();
    }
    for piece in s2_pieces.iter_mut() {
        piece.zeroize();
    }

    Ok(party_shares)
}

fn polyvecl_equal(a: &PolyVecL, b: &PolyVecL) -> bool {
    for l in 0..L {
        for c in 0..N {
            if a.polys[l].coeffs[c] != b.polys[l].coeffs[c] {
                return false;
            }
        }
    }
    true
}

fn polyveck_equal(a: &PolyVecK, b: &PolyVecK) -> bool {
    for k in 0..K {
        for c in 0..N {
            if a.polys[k].coeffs[c] != b.polys[k].coeffs[c] {
                return false;
            }
        }
    }
    true
}

/// Reconstruct the secret key (s₁, s₂) from a qualifying set of parties.
///
/// A qualifying set is any subset of size ≥ T. The reconstruction works
/// by summing exactly one copy of each RSS share piece: the union of
/// subsets held by the qualifying parties covers all (N-T+1)-subsets.
///
/// # Arguments
/// * `party_shares` — the shares of the qualifying parties
/// * `n` — total number of parties
/// * `t` — threshold
///
/// # Returns
/// The reconstructed (s₁, s₂).
pub fn reconstruct(
    party_shares: &[&PartyKeyShare],
    n: usize,
    t: usize,
) -> Result<(PolyVecL, PolyVecK), Error> {
    if n > MAX_PARTIES || t < 2 || t > n {
        return Err(Error::InvalidParameters);
    }
    if party_shares.len() < t {
        return Err(Error::InsufficientResponses);
    }

    // Count expected number of subsets C(n, n-t+1)
    let subset_size = n - t + 1;
    let num_subsets = {
        let mut result = 1usize;
        for i in 0..subset_size {
            result = result * (n - i) / (i + 1);
        }
        result
    };

    // Track which subsets we've already accumulated
    let mut seen = vec![false; num_subsets];
    let mut first_piece_refs: Vec<Option<(&PolyVecL, &PolyVecK)>> = vec![None; num_subsets];
    let mut s1_acc = PolyVecL::zero();
    let mut s2_acc = PolyVecK::zero();

    for ps in party_shares {
        if ps.n != n || ps.t != t {
            return Err(Error::InvalidParameters);
        }
        if ps.s1_shares.len() != ps.s2_shares.len() || ps.s1_shares.len() != ps.subset_indices.len()
        {
            return Err(Error::InvalidShare);
        }
        for (local_idx, &global_subset_idx) in ps.subset_indices.iter().enumerate() {
            if global_subset_idx >= num_subsets {
                return Err(Error::InvalidShare);
            }
            if local_idx >= ps.s1_shares.len() || local_idx >= ps.s2_shares.len() {
                return Err(Error::InvalidShare);
            }
            let s1_piece = &ps.s1_shares[local_idx];
            let s2_piece = &ps.s2_shares[local_idx];

            if !seen[global_subset_idx] {
                // This party holds a copy of this subset's share.
                // We only need one copy, so mark it as seen.
                s1_acc.add_assign(s1_piece);
                s2_acc.add_assign(s2_piece);
                seen[global_subset_idx] = true;
                first_piece_refs[global_subset_idx] = Some((s1_piece, s2_piece));
            } else if let Some((expected_s1, expected_s2)) = first_piece_refs[global_subset_idx] {
                // Replicas must be consistent regardless of input ordering.
                if !polyvecl_equal(expected_s1, s1_piece) || !polyveck_equal(expected_s2, s2_piece)
                {
                    return Err(Error::InvalidShare);
                }
            } else {
                return Err(Error::InvalidShare);
            }
        }
    }

    // Verify all subsets were covered
    if !seen.iter().all(|&s| s) {
        return Err(Error::InsufficientResponses);
    }

    Ok((s1_acc, s2_acc))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_enumerate_subsets() {
        // C(3, 2) = 3 subsets
        let subsets = enumerate_subsets(3, 2);
        assert_eq!(subsets.len(), 3);
        assert_eq!(subsets[0], vec![0, 1]);
        assert_eq!(subsets[1], vec![0, 2]);
        assert_eq!(subsets[2], vec![1, 2]);

        // C(4, 3) = 4 subsets (for threshold T=2, subset size is N-T+1=3)
        let subsets = enumerate_subsets(4, 2);
        assert_eq!(subsets.len(), 4);

        // C(5, 3) = 10 subsets
        let subsets = enumerate_subsets(5, 3);
        assert_eq!(subsets.len(), 10);
    }

    #[test]
    fn test_distribute_and_reconstruct_basic() {
        let mut rng = StdRng::seed_from_u64(42);

        // Create a known secret
        let mut s1 = PolyVecL::zero();
        let mut s2 = PolyVecK::zero();
        for l in 0..L {
            for c in 0..N {
                s1.polys[l].coeffs[c] = ((l * N + c) % Q as usize) as i32;
            }
        }
        for k in 0..K {
            for c in 0..N {
                s2.polys[k].coeffs[c] = ((k * N + c + 1000) % Q as usize) as i32;
            }
        }

        // Distribute with N=3, T=2
        let shares = distribute_key(&s1, &s2, 3, 2, &mut rng).unwrap();
        assert_eq!(shares.len(), 3);

        // Reconstruct with parties {0, 1}
        let qualifying = vec![&shares[0], &shares[1]];
        let (r_s1, _r_s2) = reconstruct(&qualifying, 3, 2).unwrap();

        // Reduce and compare
        for l in 0..L {
            for c in 0..N {
                let expected = s1.polys[l].coeffs[c] % Q;
                let mut got = r_s1.polys[l].coeffs[c] % Q;
                if got < 0 {
                    got += Q;
                }
                let mut exp = expected;
                if exp < 0 {
                    exp += Q;
                }
                assert_eq!(got, exp, "s1 mismatch at [{l}][{c}]");
            }
        }
    }
}
