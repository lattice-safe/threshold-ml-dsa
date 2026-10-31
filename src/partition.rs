//! Balanced partition of RSS shares among active signers.
//!
//! Implements `RSSRecover` from ePrint 2026/013, Algorithm 6.
//! For a given active signer set, computes a partition of all
//! C(N, N-T+1) subset secrets among T signers, minimizing the
//! maximum load per party.
//!
//! Subsets are encoded as bitmasks: bit `i` is set if party `i`
//! is a member of the subset. This matches the reference Go
//! implementation.

extern crate alloc;
use crate::error::Error;
use crate::params::MAX_PARTIES;
use alloc::vec::Vec;

/// Compute the balanced partition of shares for the given active signer set.
///
/// Returns a `Vec` of length `T`, where element `i` contains the list of
/// subset bitmasks assigned to the `i`-th active signer (`active[i]`).
///
/// # Arguments
/// * `active` — sorted list of active party IDs (0-indexed), length = T
/// * `n` — total number of parties N
/// * `t` — threshold T (number of active signers = `active.len()`)
///
/// Returns `Error::InvalidParameters` if the active set or `(t, n)` is invalid.
pub fn rss_recover(active: &[u8], n: u8, t: u8) -> Result<Vec<Vec<u8>>, Error> {
    if active.len() != t as usize || t < 2 || t > n || (n as usize) > MAX_PARTIES {
        return Err(Error::InvalidParameters);
    }
    let mut prev: Option<u8> = None;
    for &id in active {
        if id >= n {
            return Err(Error::InvalidParameters);
        }
        if let Some(p) = prev {
            // Require sorted, duplicate-free signer sets.
            if id <= p {
                return Err(Error::InvalidParameters);
            }
        }
        prev = Some(id);
    }

    // Base case: T == N → each party has exactly one share (its own singleton)
    if t == n {
        let mut result = Vec::with_capacity(t as usize);
        for &party in active {
            // When T == N, each subset has size 1 (N-T+1 = 1), so each
            // subset is a single party bitmask
            result.push(alloc::vec![1u8 << party]);
        }
        return Ok(result);
    }

    // Hardcoded optimal partitions for canonical active set {0, 1, ..., T-1}
    // from ePrint 2026/013, Algorithm 6, lines 3-22.
    //
    // Subset bitmasks: e.g., "01" = parties {0,1} = bitmask 0b11 = 3
    //                        "02" = parties {0,2} = bitmask 0b101 = 5
    // The Go reference uses these exact values.
    let sharing = canonical_partition(t, n).ok_or(Error::InvalidParameters)?;

    // Build the permutation φ: canonical party index → actual party ID
    // (Algorithm 6, lines 23-31)
    let mut perm = alloc::vec![0u8; n as usize];
    let mut i1 = 0usize; // next slot for active parties
    let mut i2 = t as usize; // next slot for inactive parties
    for j in 0..n {
        if active.contains(&j) {
            perm[i1] = j;
            i1 += 1;
        } else {
            perm[i2] = j;
            i2 += 1;
        }
    }

    // Apply the permutation to translate canonical subset bitmasks
    // to actual subset bitmasks (Algorithm 6, lines 32-33)
    let mut result = Vec::with_capacity(t as usize);
    for party_shares in &sharing {
        let mut translated = Vec::with_capacity(party_shares.len());
        for &canonical_mask in party_shares {
            let mut actual_mask: u8 = 0;
            for bit in 0..n {
                if canonical_mask & (1 << bit) != 0 {
                    actual_mask |= 1 << perm[bit as usize];
                }
            }
            translated.push(actual_mask);
        }
        result.push(translated);
    }

    Ok(result)
}

/// Hardcoded optimal partition for the canonical active set `{0, …, T−1}`
/// (Algorithm 6, lines 3-22), or `None` for an unsupported `(T, N)` with
/// `T < N`.
fn canonical_partition(t: u8, n: u8) -> Option<Vec<Vec<u8>>> {
    Some(match (t, n) {
        (2, 3) => alloc::vec![
            alloc::vec![3, 5], // party 0: subsets {0,1}=3, {0,2}=5
            alloc::vec![6],    // party 1: subset  {1,2}=6
        ],
        (2, 4) => alloc::vec![
            alloc::vec![11, 13], // {0,1,3}=11, {0,2,3}=13
            alloc::vec![7, 14],  // {0,1,2}=7,  {1,2,3}=14
        ],
        (3, 4) => alloc::vec![
            alloc::vec![3, 9],  // {0,1}=3, {0,3}=9
            alloc::vec![6, 10], // {1,2}=6, {1,3}=10
            alloc::vec![12, 5], // {2,3}=12, {0,2}=5
        ],
        (2, 5) => alloc::vec![
            alloc::vec![27, 29, 23], // {0,1,3,4}=27, {0,2,3,4}=29, {0,1,2,4}=23
            alloc::vec![30, 15],     // {1,2,3,4}=30, {0,1,2,3}=15
        ],
        (3, 5) => alloc::vec![
            alloc::vec![25, 11, 19, 13], // {0,3,4}=25, {0,1,3}=11, {0,1,4}=19, {0,2,3}=13
            alloc::vec![7, 14, 22, 26],  // {0,1,2}=7, {1,2,3}=14, {1,2,4}=22, {1,3,4}=26
            alloc::vec![28, 21],         // {2,3,4}=28, {0,2,4}=21
        ],
        (4, 5) => alloc::vec![
            alloc::vec![3, 9, 17],  // {0,1}=3, {0,3}=9, {0,4}=17
            alloc::vec![6, 10, 18], // {1,2}=6, {1,3}=10, {1,4}=18
            alloc::vec![12, 5, 20], // {2,3}=12, {0,2}=5, {2,4}=20
            alloc::vec![24],        // {3,4}=24
        ],
        (2, 6) => alloc::vec![
            alloc::vec![61, 47, 55], // {0,2,3,4,5}=61, {0,1,2,3,5}=47, {0,1,2,4,5}=55
            alloc::vec![62, 31, 59], // {1,2,3,4,5}=62, {0,1,2,3,4}=31, {0,1,3,4,5}=59
        ],
        (3, 6) => alloc::vec![
            alloc::vec![27, 23, 43, 57, 39],
            alloc::vec![51, 58, 46, 30, 54],
            alloc::vec![45, 53, 29, 15, 60],
        ],
        (4, 6) => alloc::vec![
            alloc::vec![19, 13, 35, 7, 49],
            alloc::vec![42, 26, 38, 50, 22],
            alloc::vec![52, 21, 44, 28, 37],
            alloc::vec![25, 11, 14, 56, 41],
        ],
        (5, 6) => alloc::vec![
            alloc::vec![3, 5, 33],   // {0,1}, {0,2}, {0,5}
            alloc::vec![6, 10, 34],  // {1,2}, {1,3}, {1,5}
            alloc::vec![12, 20, 36], // {2,3}, {2,4}, {2,5}
            alloc::vec![9, 24, 40],  // {0,3}, {3,4}, {3,5}
            alloc::vec![48, 17, 18], // {4,5}, {0,4}, {1,4}
        ],
        // N = 7, 8: max-flow balanced partitions (same procedure as the
        // reference `params/recover.py`; it reproduces the N ≤ 6 tables above).
        (2, 7) => alloc::vec![
            alloc::vec![125, 95, 111, 119],
            alloc::vec![63, 126, 123],
        ],
        (3, 7) => alloc::vec![
            alloc::vec![79, 87, 103, 59, 121, 47, 55],
            alloc::vec![91, 118, 122, 94, 110, 107, 115],
            alloc::vec![109, 117, 93, 61, 31, 62, 124],
        ],
        (4, 7) => alloc::vec![
            alloc::vec![27, 51, 75, 23, 113, 43, 71, 83, 39],
            alloc::vec![58, 114, 46, 99, 86, 102, 78, 90, 54],
            alloc::vec![45, 53, 116, 92, 101, 108, 60, 77, 85],
            alloc::vec![120, 105, 89, 30, 29, 15, 57, 106],
        ],
        (5, 7) => alloc::vec![
            alloc::vec![19, 21, 97, 13, 67, 35, 7],
            alloc::vec![42, 70, 14, 26, 98, 38, 74],
            alloc::vec![52, 76, 44, 69, 37, 84, 100],
            alloc::vec![88, 25, 104, 11, 73, 56, 41],
            alloc::vec![82, 28, 81, 112, 50, 49, 22],
        ],
        (6, 7) => alloc::vec![
            alloc::vec![3, 5, 17, 65],
            alloc::vec![66, 6, 10, 34],
            alloc::vec![12, 20, 36, 68],
            alloc::vec![9, 24, 40, 72],
            alloc::vec![48, 80, 18],
            alloc::vec![96, 33],
        ],
        (2, 8) => alloc::vec![
            alloc::vec![253, 223, 239, 191],
            alloc::vec![247, 127, 251, 254],
        ],
        (3, 8) => alloc::vec![
            alloc::vec![119, 159, 215, 175, 207, 95, 111, 249, 183, 231],
            alloc::vec![235, 250, 190, 222, 238, 246, 243, 187, 219, 123],
            alloc::vec![63, 245, 125, 221, 252, 237, 126, 189],
        ],
        (4, 8) => alloc::vec![
            alloc::vec![151, 167, 91, 79, 199, 87, 103, 59, 241, 47, 55, 107, 155, 143],
            alloc::vec![227, 182, 242, 158, 110, 174, 118, 203, 94, 179, 211, 62, 206, 115],
            alloc::vec![229, 244, 157, 109, 230, 205, 181, 213, 117, 93, 188, 214, 124, 173],
            alloc::vec![248, 236, 186, 234, 185, 218, 171, 121, 217, 122, 220, 61, 31, 233],
        ],
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::num_subsets;

    /// Verify that the partition covers all subsets exactly once
    fn verify_partition_covers_all(partition: &[Vec<u8>], n: u8, t: u8) {
        let mut all_subsets: Vec<u8> = partition.iter().flat_map(|v| v.iter().copied()).collect();
        all_subsets.sort();
        all_subsets.dedup();

        let expected_count = num_subsets(n, t);
        assert_eq!(all_subsets.len(), expected_count, "(T={t}, N={n}): duplicate or missing subsets");
    }

    /// Verify each party only receives subsets it belongs to
    fn verify_party_membership(active: &[u8], partition: &[Vec<u8>]) {
        for (i, shares) in partition.iter().enumerate() {
            let party = active[i];
            for &mask in shares {
                assert!(
                    mask & (1 << party) != 0,
                    "Party {} assigned subset 0b{:06b} which it doesn't belong to",
                    party,
                    mask
                );
            }
        }
    }

    #[test]
    fn test_partition_3_2() {
        let active = [0, 1];
        let p = rss_recover(&active, 3, 2).unwrap();
        assert_eq!(p.len(), 2);
        verify_partition_covers_all(&p, 3, 2);
        verify_party_membership(&active, &p);
    }

    #[test]
    fn test_partition_4_3() {
        let active = [0, 1, 2];
        let p = rss_recover(&active, 4, 3).unwrap();
        assert_eq!(p.len(), 3);
        verify_partition_covers_all(&p, 4, 3);
        verify_party_membership(&active, &p);
    }

    #[test]
    fn test_partition_permuted_active_set() {
        // Use a non-canonical active set {1, 3, 4} for (T=3, N=5)
        let active = [1, 3, 4];
        let p = rss_recover(&active, 5, 3).unwrap();
        assert_eq!(p.len(), 3);
        verify_partition_covers_all(&p, 5, 3);
        verify_party_membership(&active, &p);
    }

    #[test]
    fn test_partition_6_6_base_case() {
        let active = [0, 1, 2, 3, 4, 5];
        let p = rss_recover(&active, 6, 6).unwrap();
        assert_eq!(p.len(), 6);
        // T == N: each party gets exactly one singleton subset
        for (i, shares) in p.iter().enumerate() {
            assert_eq!(shares.len(), 1);
            assert_eq!(shares[0], 1u8 << active[i]);
        }
    }

    #[test]
    fn test_all_valid_configs() {
        let configs: &[(u8, u8)] = &[
            (2, 2),
            (2, 3),
            (3, 3),
            (2, 4),
            (3, 4),
            (4, 4),
            (2, 5),
            (3, 5),
            (4, 5),
            (5, 5),
            (2, 6),
            (3, 6),
            (4, 6),
            (5, 6),
            (6, 6),
            (2, 7),
            (3, 7),
            (4, 7),
            (5, 7),
            (6, 7),
            (7, 7),
            (2, 8),
            (3, 8),
            (4, 8),
            (8, 8),
        ];
        for &(t, n) in configs {
            // Use canonical active set {0, 1, ..., T-1}
            let active: Vec<u8> = (0..t).collect();
            let p = rss_recover(&active, n, t).unwrap();
            assert_eq!(
                p.len(),
                t as usize,
                "Wrong partition size for ({}, {})",
                t,
                n
            );
            verify_partition_covers_all(&p, n, t);
            verify_party_membership(&active, &p);
        }
    }

    #[test]
    fn test_invalid_active_set_rejected() {
        // Duplicate party id.
        assert_eq!(rss_recover(&[0, 0], 3, 2), Err(Error::InvalidParameters));
        // Unsorted.
        assert_eq!(rss_recover(&[2, 1], 3, 2), Err(Error::InvalidParameters));
        // Out of range.
        assert_eq!(rss_recover(&[0, 3], 3, 2), Err(Error::InvalidParameters));
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;

    #[test]
    fn invalid_inputs_are_rejected() {
        assert!(rss_recover(&[0], 3, 2).is_err()); // wrong length
        assert!(rss_recover(&[0, 1], 9, 2).is_err()); // N too large
        assert!(rss_recover(&[0, 3], 3, 2).is_err()); // id ≥ N
        assert!(rss_recover(&[1, 0], 3, 2).is_err()); // unsorted
        assert!(canonical_partition(2, 9).is_none());
        assert!(canonical_partition(5, 8).is_none()); // unsupported (impractical K)
        assert!(canonical_partition(3, 3).is_none()); // T = N handled by the caller
    }
}
