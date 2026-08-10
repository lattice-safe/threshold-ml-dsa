# threshold-ml-dsa

**Threshold ML-DSA (FIPS 204) — Paper-Faithful Implementation of ePrint 2026/013**

A `#![no_std]`-compatible Rust implementation of threshold ML-DSA based on the Mithril scheme ([ePrint 2026/013](https://eprint.iacr.org/2026/013)). Threshold signatures are **bit-for-bit compatible** with standard FIPS 204 verifiers.

[![Crates.io](https://img.shields.io/crates/v/threshold-ml-dsa.svg)](https://crates.io/crates/threshold-ml-dsa)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org)

---

## Overview

Standard threshold signature schemes based on Shamir secret sharing introduce large Lagrange interpolation coefficients that blow up lattice coefficient sizes, breaking the short-vector requirements of ML-DSA. The Mithril scheme solves this with **Replicated Secret Sharing (RSS)** and **hyperball-based local rejection sampling**.

### Key Features & Capabilities

| Feature | Description | Reference |
|---|---|---|
| **FIPS 204 Verifier Compatible** | Output signatures pass any standard, unmodified ML-DSA verifier | FIPS 204 |
| **Fresh Threshold Keygen** | Independent secrets per subset — dealer-based key generation | Figure 4 |
| **A Posteriori Key Sharing** | Convert existing standard ML-DSA secret keys into threshold shares | §4.2 |
| **Distributed Key Generation** | 2-round commitment-based DKG without relying on a trusted dealer | §4.3 |
| **K-Parallel Repetitions** | Amortized rejection sampling via K parallel commitment slots per round | §3.2 |
| **Hyperball Rejection** | Rényi-divergence-safe L₂ norm rejection via FVec + Box-Muller | §2.7 |
| **Balanced RSS Partition** | Algorithm 6 (`RSSRecover`) + dynamic greedy partition solver for $N \le 8$ | Algorithm 6 |
| **Multi-Security-Level Ready** | Parameter modules and standard verifiers for ML-DSA-44, 65, and 87 | FIPS 204 |
| **Fail-Closed Verification** | Returned signatures are validated by the standard verifier before return | SDK Layer |
| **Single-Use Nonces** | Round 3 consumes nonce state by value — prevents replay | Security Hardening |
| **Zeroize-on-Drop** | Sensitive floats, shares, and private keys wiped on drop | Security Hardening |
| **Zero Unsafe Code** | 100% safe Rust | Hardening |
| **`#![no_std]` Compatible** | Suitable for embedded, mobile, and TEE environments | Embedded |

---

## Architecture

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   params.rs     │───▶│    poly.rs      │───▶│    rss.rs       │
│ ML-DSA-44/65/87 │    │ NTT, norms,     │    │ Fresh keygen    │
│ ThresholdParams │    │ SHAKE sampling  │    │ per subset      │
└────────┬────────┘    └────────┬────────┘    └────────┬────────┘
         │                      │                      │
         │             ┌────────▼────────┐    ┌────────▼────────┐
         │             │  aposteriori.rs │    │    dkg.rs       │
         │             │ Existing key    │    │ 2-round DKG     │
         │             │ sharing (§4.2)  │    │ protocol (§4.3) │
         │             └────────┬────────┘    └────────┬────────┘
         │                      │                      │
         │             ┌────────▼────────┐    ┌────────▼────────┐
         └────────────▶│    sign.rs      │◀───│  partition.rs   │
                       │ K-parallel 3rnd │    │ Alg 6 + Dynamic │
                       └────────┬────────┘    └─────────────────┘
                                │
┌─────────────────┐    ┌────────▼────────┐    ┌─────────────────┐
│    fvec.rs      │───▶│ coordinator.rs  │───▶│   verify.rs     │
│ SampleHyperball │    │ K-parallel      │    │  dilithium-rs   │
│ + L₂ Excess     │    │ Combine         │    │  FIPS 204       │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

---

## Quick Start

### 1. Fresh Key Generation (Trusted Dealer)

```rust
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use rand::rngs::OsRng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = OsRng;
    let seed = [42u8; 32]; // Use secure seed in production

    // Create a 2-of-3 threshold SDK with fresh key generation
    let sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 3, 100)?;

    // Sign message with active parties 0 and 1
    let msg = b"Hello, threshold ML-DSA!";
    let active = [0u8, 1];
    let sig = sdk.threshold_sign(&active, msg, &mut rng)?;

    // Verify with standard ML-DSA-44 verifier
    assert!(sdk.verify(msg, &sig));

    Ok(())
}
```

### 2. Distributed Key Generation (DKG — No Trusted Dealer)

```rust
use threshold_ml_dsa::dkg;
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use rand::rngs::OsRng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = OsRng;
    let (t, n) = (2u8, 3u8);

    // Round 1: Each party generates a contribution seed and commitment hash
    let mut contribs = Vec::new();
    let mut commits = Vec::new();
    for i in 0..n {
        let (c, h) = dkg::dkg_round1(i, &mut rng)?;
        contribs.push((i, c));
        commits.push((i, h));
    }

    // Round 2: Verify commitments and initialize SDK from DKG
    let sdk = ThresholdMlDsa44Sdk::from_dkg(&contribs, &commits, t, n, 100)?;

    // Sign and verify
    let msg = b"Signed via trustless DKG";
    let sig = sdk.threshold_sign(&[0, 2], msg, &mut rng)?;
    assert!(sdk.verify(msg, &sig));

    Ok(())
}
```

### 3. A Posteriori Key Sharing (Splitting an Existing Key)

```rust
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use threshold_ml_dsa::verify;
use rand::rngs::OsRng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = OsRng;
    let seed = [99u8; 32];

    // Generate standard ML-DSA-44 key pair
    let (pk_bytes, sk_bytes) = verify::keygen(&seed);

    // Split existing secret key into 3-of-4 threshold shares (preserves public key)
    let sdk = ThresholdMlDsa44Sdk::from_existing_key(&sk_bytes, &[77u8; 32], 3, 4, 100)?;

    // Derived threshold public key is bit-for-bit identical to original public key
    assert_eq!(sdk.pk(), &pk_bytes[..]);

    // Threshold sign with parties 0, 1, 2
    let msg = b"Splitting an existing FIPS 204 key";
    let sig = sdk.threshold_sign(&[0, 1, 2], msg, &mut rng)?;
    assert!(verify::verify(&sig, msg, &pk_bytes));

    Ok(())
}
```

---

## Supported Threshold Configurations

The implementation includes paper-exact parameters for all 15 $(T, N)$ sets from ePrint 2026/013 Figure 8, and extends support up to $N \le 8$:

| $(T, N)$ | Parallel Reps $K$ | Target Radius $r$ | Randomness Radius $r_1$ | Expansion $\nu$ |
|---|---|---|---|---|
| (2, 2) | 2 | 252,778 | 252,833 | 3.0 |
| (2, 3) | 3 | 310,060 | 310,138 | 3.0 |
| (3, 3) | 4 | 246,490 | 246,546 | 3.0 |
| (2, 4) | 3 | 305,919 | 305,997 | 3.0 |
| (3, 4) | 7 | 279,235 | 279,314 | 3.0 |
| (4, 4) | 8 | 243,463 | 243,519 | 3.0 |
| (2, 5) | 3 | 285,363 | 285,459 | 3.0 |
| (3, 5) | 14 | 282,800 | 282,912 | 3.0 |
| (4, 5) | 30 | 259,427 | 259,526 | 3.0 |
| (5, 5) | 16 | 239,924 | 239,981 | 3.0 |
| (2, 6) | 4 | 300,265 | 300,362 | 3.0 |
| (3, 6) | 19 | 277,014 | 277,139 | 3.0 |
| (4, 6) | 74 | 268,705 | 268,831 | 3.0 |
| (5, 6) | 100 | 250,590 | 250,686 | 3.0 |
| (6, 6) | 37 | 219,245 | 219,301 | 3.0 |
| **(T, 7)** | 5 – 120 | Extended | Extended | 3.0 |
| **(T, 8)** | 6 – 150 | Extended | Extended | 3.0 |

---

## Performance & Benchmarks

Measured on Apple Silicon (macOS `aarch64`), compiled with `cargo bench` (`opt-level = 3`):

### 1. Threshold Signing Latency (`threshold_sign`)

| $(T, N)$ Configuration | Parallel Slots $K$ | Average Latency | Status |
|---|---|---|---|
| **(2, 2)** | 2 | **1.39 ms** / op | Sub-2ms |
| **(2, 3)** | 3 | **998.69 µs** / op | **Sub-millisecond (< 1 ms)** |
| **(3, 3)** | 4 | **1.87 ms** / op | Sub-2ms |
| **(2, 4)** | 3 | **1.00 ms** / op | Sub-millisecond |
| **(3, 4)** | 7 | **3.04 ms** / op | ~3 ms |
| **(4, 4)** | 8 | **4.74 ms** / op | ~4.7 ms |
| **(3, 5)** | 14 | **5.86 ms** / op | ~5.8 ms |
| **(5, 5)** | 16 | **11.98 ms** / op | ~12 ms |
| **(6, 6)** | 37 | **34.94 ms** / op | ~35 ms |
| **(4, 6)** | 74 | **42.04 ms** / op | ~42 ms |
| **(5, 6)** | 100 | **72.78 ms** / op | ~72 ms |

### 2. Key Generation & Setup Latency

| Operation | Setup Type | (T, N) | Latency |
|---|---|---|---|
| **Fresh Keygen** (`from_seed`) | Trusted Dealer | (2, 3) | **167.37 µs** / op |
| **Fresh Keygen** (`from_seed`) | Trusted Dealer | (3, 4) | **193.42 µs** / op |
| **Fresh Keygen** (`from_seed`) | Trusted Dealer | (5, 6) | **403.75 µs** / op |
| **A Posteriori Sharing** (`from_existing_key`) | Key Decomposition | (2, 3) | **113.55 µs** / op |
| **A Posteriori Sharing** (`from_existing_key`) | Key Decomposition | (3, 4) | **185.67 µs** / op |
| **A Posteriori Sharing** (`from_existing_key`) | Key Decomposition | (5, 6) | **397.62 µs** / op |
| **Distributed Key Gen** (`from_dkg`) | 2-Round Trustless | (2, 3) | **120.64 µs** / op |
| **Distributed Key Gen** (`from_dkg`) | 2-Round Trustless | (3, 4) | **192.71 µs** / op |
| **Distributed Key Gen** (`from_dkg`) | 2-Round Trustless | (5, 6) | **417.59 µs** / op |

### 3. FIPS 204 Signature Verification

| Security Level | Parameter Set | Public Key Size | Signature Size | Verification Time |
|---|---|---|---|---|
| **Category 2** | ML-DSA-44 (Threshold Output) | 1,312 B | 2,420 B | **50.14 µs** / op |
| **Category 3** | ML-DSA-65 (Standard Verifier) | 1,952 B | 3,309 B | **79.68 µs** / op |
| **Category 5** | ML-DSA-87 (Standard Verifier) | 2,592 B | 4,627 B | **132.00 µs** / op |

### Running Benchmarks Locally

```bash
cargo bench --bench benchmarks -- --nocapture
```

---

## Crate Modules

| Module | Description |
|---|---|
| `params` | FIPS 204 constants for ML-DSA-44, 65, 87 + `ThresholdParams` lookup table |
| `poly` | Polynomial ring arithmetic: NTT, INTT, Montgomery reduction, norms, SHAKE-256 sampling |
| `fvec` | Floating-point vector (`FVec`) + Box-Muller `SampleHyperball` + $L_2$ norm rejection |
| `partition` | Algorithm 6 (`RSSRecover`) + dynamic greedy partition solver for $N \le 8$ |
| `rss` | `keygen_from_seed()`: fresh independent secrets per subset (Figure 4) |
| `aposteriori` | `share_existing_key()`: split an existing standard ML-DSA-44 secret key (§4.2) |
| `dkg` | 2-round commitment-based Distributed Key Generation protocol (§4.3) |
| `sign` | 3-round party protocol: K-parallel commit, reveal, respond |
| `coordinator` | K-parallel `Combine`: $\delta$-norm check, hint generation, FIPS 204 packing |
| `sdk` | High-level API: `from_seed()`, `from_dkg()`, `from_existing_key()`, `threshold_sign()`, `verify()` |
| `verify` | Standard FIPS 204 verification via `dilithium-rs` (ML-DSA-44, ML-DSA-65, ML-DSA-87) |
| `error` | `no_std`-compatible error handling |

---

## Testing & Verification

```bash
# Run complete test suite (124 tests)
cargo test

# Run full-coverage integration suite (DKG + Sign + Verify)
cargo test --test dkg_sign_verify_full_coverage

# Run end-to-end threshold signing tests
cargo test --test v03_tests
```

### Test Coverage Summary

| Test Suite | Tests | Description | Result |
|---|---|---|---|
| **Full Coverage Integration** | 10 | DKG, A Posteriori, non-canonical active sets, tampering, multi-level ML-DSA | ✅ **Passed** |
| **Internal Unit Tests** | 44 | Poly arithmetic, RSS keygen, DKG, aposteriori, partitions, sign rounds | ✅ **Passed** |
| **Ported Coverage Tests** | 51 | NIST KAT alignment, poly norms, packing, challenge weight, NTT roundtrip | ✅ **Passed** |
| **v0.3 E2E Integration** | 12 | Threshold signing for (2,2), (2,3), (3,3), (3,4), (4,6), (5,6) | ✅ **Passed** |
| **Component Specific Tests** | 6 | `ct0`, `delta`, `fvec`, `ntt`, `pack_overflow` | ✅ **Passed** |
| **Doc Tests** | 1 | Inline documentation code block compilation | ✅ **Passed** |
| **Total** | **124** | **All 124 tests passing cleanly** | **100%** |

---

## Security Hardening

- **No Key Reuse**: Fresh keygen samples independent secrets per subset.
- **Single-Use Nonces**: `round3()` consumes `StRound1` by value — Rust ownership guarantees nonce randomness cannot be re-used across different challenges.
- **Cryptographic Reveal Witness**: `verify_all_round2_reveals()` returns a bound witness token that `round3()` requires, enforcing round ordering.
- **Hedged Nonce Generation**: Round 1 binds PRNG entropy + long-term key + transcript context + party ID.
- **Fail-Closed SDK**: Output signatures are verified by the pure-Rust `dilithium-rs` FIPS 204 verifier before returning to caller.
- **Zeroize-on-Drop**: `FVec`, `StRound1`, `StRound2`, `ThresholdPrivateKey`, `Share`, and `DkgContribution` implement `Zeroize` and `Drop`.
- **100% Safe Rust**: Zero `unsafe` blocks across the codebase.

---

## References

- [FIPS 204 — ML-DSA](https://csrc.nist.gov/pubs/fips/204/final) — Module-Lattice-Based Digital Signature Standard
- [ePrint 2026/013](https://eprint.iacr.org/2026/013) — Efficient Threshold ML-DSA (Mithril Scheme)
- [`dilithium-rs`](https://crates.io/crates/dilithium-rs) — Pure-Rust FIPS 204 implementation

## License

Licensed under the [MIT License](LICENSE).
