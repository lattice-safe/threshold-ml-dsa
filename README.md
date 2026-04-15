# threshold-ml-dsa

**Threshold ML-DSA (FIPS 204) — Paper-Faithful Implementation of ePrint 2026/013**

A `#![no_std]`-compatible Rust implementation of threshold ML-DSA-44 based on the Mithril scheme ([ePrint 2026/013](https://eprint.iacr.org/2026/013)). Threshold signatures are **bit-for-bit compatible** with standard FIPS 204 verifiers.

[![Crates.io](https://img.shields.io/crates/v/threshold-ml-dsa.svg)](https://crates.io/crates/threshold-ml-dsa)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org)

---

## Overview

Standard threshold signature schemes based on Shamir secret sharing introduce large Lagrange interpolation coefficients that blow up lattice coefficient sizes, breaking the short-vector requirements of ML-DSA. The Mithril scheme solves this with **Replicated Secret Sharing (RSS)** and **hyperball-based local rejection sampling**.

### Key Properties

| Property | Description |
|----------|-------------|
| **FIPS 204 Compatible** | Output signatures pass any unmodified ML-DSA-44 verifier |
| **Fail-Closed** | Invalid aggregates are never returned as successful signatures |
| **Paper-Faithful** | Exact parameter sets from ePrint 2026/013, Figure 8 |
| **K-Parallel Repetitions** | Multiple commitment slots per round for amortized acceptance |
| **Hyperball Rejection** | Rényi-divergence-safe L₂ norm rejection via FVec + Box-Muller |
| **Balanced Partitions** | Algorithm 6 (RSSRecover) for optimal share assignment |
| **Fresh Keygen** | Independent secrets per subset — no existing key decomposition |
| **Single-Use Nonces** | Round-3 consumes nonce state by value — prevents replay |
| **Zeroize-on-Drop** | FVec, StRound1, StRound2, ThresholdPrivateKey all wiped on drop |
| **Zero Unsafe** | 100% safe Rust (zeroize via `zeroize` crate, not `write_volatile`) |
| **`#![no_std]`** | Suitable for embedded and TEE environments |

## Architecture

```
┌──────────────┐    ┌──────────────┐    ┌──────────────┐
│  params.rs   │───▶│   poly.rs    │───▶│   rss.rs     │
│ ML-DSA-44 +  │    │ NTT, norms,  │    │ Fresh keygen │
│ Figure 8     │    │ SHAKE sample │    │ per subset   │
│ ThresholdP.  │    └──────┬───────┘    └──────┬───────┘
└──────────────┘           │                    │
                    ┌──────▼───────┐    ┌───────▼──────┐
                    │  sign.rs     │◀───│partition.rs  │
                    │  K-parallel  │    │ Algorithm 6  │
                    │  3-round     │    │ balanced RSS │
                    └──────┬───────┘    └──────────────┘
                           │
┌──────────────┐    ┌──────▼───────┐    ┌──────────────┐
│   fvec.rs    │───▶│coordinator.rs│───▶│  verify.rs   │
│ SampleHyper  │    │  K-parallel  │    │  dilithium-rs│
│ ball + L₂    │    │  Combine     │    │  FIPS 204    │
└──────────────┘    └──────────────┘    └──────────────┘
```

## Quick Start

```rust
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use rand::rngs::OsRng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = OsRng;
    let seed = [42u8; 32]; // use OsRng in production

    // 1. Create a 2-of-3 threshold SDK (fresh keygen from seed)
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 3, 100)?;

    // 2. Sign with parties 0 and 1 (any 2 of 3)
    let msg = b"Hello, threshold ML-DSA!";
    let active = [0u8, 1];
    let sig = sdk.threshold_sign(&active, msg, &mut rng)?;

    // 3. Verify with any standard ML-DSA-44 verifier
    assert!(sdk.verify(msg, &sig));

    Ok(())
}
```

## Protocol (3 Rounds + K-Parallel)

The v0.3 protocol follows ePrint 2026/013 exactly:

| Round | Party Action | Coordinator Action |
|-------|-------------|-------------------|
| **1 — Commit** | Sample K hyperball vectors via `SampleHyperball(r₁, ν)`, compute K commitments `w_{i,k} = A·r_k + e_k`, broadcast `H(tag ‖ tr ‖ id ‖ act ‖ session ‖ μ ‖ w_packed)` | Collect commitment hashes |
| **2 — Reveal** | Send full K commitment vectors | Verify reveal/hash binding, aggregate `w_k = Σ w_{i,k}` for each slot k |
| **3 — Respond** | Recover partial secret via `RSSRecover(active)`, compute `z_{i,k} = c·s_I + (r_k, e_k)`, apply `FVec.Excess(r, ν)` rejection | Aggregate `z_k = Σ z_{i,k}`, try each k: check `‖z_k‖∞ < γ₁-β`, compute δ, generate hint, pack FIPS 204 signature, verify end-to-end |

### Key Differences from v0.2

| v0.2 | v0.3 |
|------|------|
| Decompose existing ML-DSA key | **Fresh independent keygen** per subset |
| Single commitment per round | **K parallel** commitment slots |
| L₂ norm on integer coefficients | **FVec float hyperball** (Box-Muller, ν-scaling) |
| Ad-hoc share assignment | **Algorithm 6** balanced partition |
| No δ norm check | **δ = ‖Az-2^d·c·t₁-w‖∞ < γ₂** coordinator check |

## Supported Configurations

All 15 (T, N) parameter sets from ePrint 2026/013, Figure 8:

| (T, N) | K reps | Radius r | Radius r₁ | ν |
|--------|--------|----------|-----------|---|
| (2, 2) | 2 | 252,778 | 252,833 | 3 |
| (2, 3) | 3 | 310,060 | 310,138 | 3 |
| (3, 3) | 4 | 246,490 | 246,546 | 3 |
| (2, 4) | 3 | 305,919 | 305,997 | 3 |
| (3, 4) | 7 | 279,235 | 279,314 | 3 |
| (4, 4) | 8 | 243,463 | 243,519 | 3 |
| (2, 5) | 3 | 285,363 | 285,459 | 3 |
| (3, 5) | 14 | 282,800 | 282,912 | 3 |
| (4, 5) | 30 | 259,427 | 259,526 | 3 |
| (5, 5) | 16 | 239,924 | 239,981 | 3 |
| (2, 6) | 4 | 300,265 | 300,362 | 3 |
| (3, 6) | 19 | 277,014 | 277,139 | 3 |
| (4, 6) | 74 | 268,705 | 268,831 | 3 |
| (5, 6) | 100 | 250,590 | 250,686 | 3 |
| (6, 6) | 37 | 219,245 | 219,301 | 3 |

## Modules

| Module | Description |
|--------|-------------|
| `params` | ML-DSA-44 constants + `ThresholdParams` with Figure 8 values for all 15 (T,N) pairs |
| `poly` | Polynomial arithmetic: NTT, norms, SHAKE sampling, packing |
| `fvec` | Float vector + `SampleHyperball` via Box-Muller + `Excess` L₂ check |
| `partition` | Algorithm 6 (`RSSRecover`): balanced partition for arbitrary active sets |
| `rss` | `keygen_from_seed()`: fresh independent secrets per subset (Figure 4) |
| `sign` | 3-round party protocol: K-parallel commit, reveal, respond |
| `coordinator` | K-parallel `Combine`: δ check, hint generation, FIPS 204 packing |
| `sdk` | High-level API: `from_seed()`, `threshold_sign()`, `verify()` |
| `verify` | FIPS 204 verification via `dilithium-rs` |
| `error` | `no_std`-compatible error types |

## Parameters (ML-DSA-44)

| Parameter | Value | Description |
|-----------|-------|-------------|
| q | 8,380,417 | Prime modulus |
| N | 256 | Polynomial degree |
| (K, L) | (4, 4) | Matrix dimensions |
| η | 2 | Secret key coefficient bound |
| γ₁ | 2¹⁷ | Masking vector bound |
| γ₂ | (q-1)/88 = 95,232 | Decompose parameter |
| τ | 39 | Challenge weight |
| ν | 3 | Expansion factor |
| Max parties | 6 | RSS subset limit |

## Testing

```bash
# Run all tests
cargo test

# Run end-to-end signing tests
cargo test --test v03_tests
```

### Test Coverage

| Category | Tests | Coverage |
|----------|-------|----------|
| Keygen | 4 | Determinism, all (T,N) configs, subset counts, valid keys |
| Partitions | 5 | All 15 configs, balanced coverage, permuted active sets |
| Hyperball sampling | 2 | Norm bound, excess check |
| FVec roundtrip | 1 | Poly ↔ FVec centering and reconstruction |
| SDK | 3 | Creation, invalid params, duplicate/Sybil rejection |
| **End-to-end threshold sign** | **6** | **(2,2), (2,3), (3,3), (3,4), (4,6), (5,6): `threshold_sign` + FIPS 204 verify** |
| Sign unit tests | 5 | pack/unpack roundtrip, bitmask, mu determinism |
| Coordinator unit tests | 6 | aggregate commitments/responses (zero, additive, short-input rejection) |
| Poly arithmetic | 6 | Add, sub, center, norms, power2round |
| Params | 4 | Lookup, binomial, num_subsets, rejection |
| Doctest | 1 | Usage example compilation |
| **Total** | **107** | **All passing, clean `no_std`; pedantic clippy currently reports non-security style warnings (mostly tests/docs)** |

## Security Notes

- **No key decomposition**: Keys are generated fresh per subset — no existing ML-DSA secret is ever split.
- **Fail-closed signing**: Returned signatures are verified by the standard FIPS 204 verifier before return.
- **Single-use nonce state**: `round3()` takes `StRound1` by value — Rust ownership prevents replay of the same nonce randomness with different coordinator challenges.
- **Input validation**: `aggregate_commitments`, `aggregate_responses`, and `combine` validate input lengths and return `Result`, preventing panics on malformed vectors.
- **Sybil / duplicate-ID protection**: The SDK rejects duplicate or unsorted party IDs in the active set.
- **Zeroize-on-drop**: `FVec`, `StRound1`, `StRound2`, and `ThresholdPrivateKey` all wipe sensitive material on drop via the `zeroize` crate.
- **Zero unsafe code**: The entire crate is 100% safe Rust.
- **Hyperball rejection**: Uses Rényi-divergence-safe float L₂ norms; branchless final comparison and branchless Box-Muller clamp.
- **Timing model note**: Hyperball sampling uses floating-point (`f64`/`libm`) and is not strictly constant-time. Deployments should isolate party execution.
- **Balanced share assignment**: Algorithm 6 ensures each party receives an equal number of subset secrets.
- **Commitment binding**: Round-1 hashes bind `(tr, id, act, session, μ, w_packed)` to prevent cross-session replay.
- **SDK limitation**: `ThresholdMlDsa44Sdk` is an in-process orchestrator. For distributed deployments, use the low-level `sign::round*` APIs with per-peer `verify_round2_reveal()` checks and authenticated transport.

## Dependencies

| Crate | Purpose |
|-------|---------|
| `dilithium-rs` | FIPS 204 NTT, Montgomery reduction, and canonical verification |
| `sha3` | SHAKE-256 for deterministic sampling |
| `libm` | `no_std`-compatible floating-point (Box-Muller, `sqrt`, `log`) |
| `zeroize` | Secure erasure of sensitive key material |
| `subtle` | Constant-time operations |
| `rand_core` | CSPRNG trait for key generation |

## References

- [FIPS 204 — ML-DSA](https://csrc.nist.gov/pubs/fips/204/final) — Module-Lattice-Based Digital Signature Standard
- [ePrint 2026/013](https://eprint.iacr.org/2026/013) — Efficient Threshold ML-DSA (Mithril Scheme)
- [Threshold-ML-DSA (Go)](https://github.com/Threshold-ML-DSA/Threshold-ML-DSA) — Reference Go implementation
- [`dilithium-rs`](https://crates.io/crates/dilithium-rs) — Pure-Rust FIPS 204 implementation

## License

This project is licensed under the [MIT License](LICENSE).
