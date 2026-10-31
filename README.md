# threshold-ml-dsa

**Threshold ML-DSA (FIPS 204) — Paper-Faithful Implementation of ePrint 2026/013**

A `#![no_std]`-compatible Rust implementation of threshold ML-DSA based on the Mithril scheme ([ePrint 2026/013](https://eprint.iacr.org/2026/013)). Threshold signatures are **bit-for-bit compatible** with standard FIPS 204 verifiers.

[![Crates.io](https://img.shields.io/crates/v/threshold-ml-dsa.svg)](https://crates.io/crates/threshold-ml-dsa)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.79%2B-orange.svg)](https://www.rust-lang.org)

> **Security notice.** Versions before 0.5.0 are affected by
> [GHSA-hqr3-xgxm-59jr](https://github.com/lattice-safe/threshold-ml-dsa/security/advisories/GHSA-hqr3-xgxm-59jr)
> (critical: the old DKG gave every party the full key) and 0.3.1–0.4.0 also by
> [GHSA-22vh-j6fg-hm7g](https://github.com/lattice-safe/threshold-ml-dsa/security/advisories/GHSA-22vh-j6fg-hm7g)
> (rejection sampling disabled). Upgrade to 0.5.1 or later and replace affected
> keys; see [SECURITY.md](SECURITY.md) for which keys and how to report issues,
> and [CHANGELOG.md](CHANGELOG.md) for details. 0.3.1–0.4.1 are yanked.

---

## Overview

Standard threshold signature schemes based on Shamir secret sharing introduce large Lagrange interpolation coefficients that blow up lattice coefficient sizes, breaking the short-vector requirements of ML-DSA. The Mithril scheme solves this with **Replicated Secret Sharing (RSS)** and **hyperball-based local rejection sampling**.

### Key Features & Capabilities

| Feature | Description | Reference |
|---|---|---|
| **FIPS 204 Verifier Compatible** | Output signatures pass any standard, unmodified ML-DSA verifier | FIPS 204 |
| **Fresh Threshold Keygen** | Independent secrets per subset — dealer-based key generation | Fig. 5 (USENIX version) |
| **K-Parallel Repetitions** | Amortized rejection sampling via K parallel commitment slots per round | §3.2 |
| **Hyperball Rejection** | Rényi-divergence-safe L₂ norm rejection via FVec + Box-Muller | §2.7 |
| **Balanced RSS Partition** | Algorithm 6 (`RSSRecover`) for $N \le 6$ | Algorithm 6 |
| **Standard Verifiers** | Single-party verifiers for ML-DSA-44, 65, and 87 (threshold signing is ML-DSA-44 only) | FIPS 204 |
| **Fail-Closed Verification** | Returned signatures are validated by the standard verifier before return | SDK Layer |
| **Single-Use Nonces** | Rounds 2 and 3 consume their state by value — prevents replay | Security Hardening |
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
         │             │ Gaussian coset  │    │ 5-round DKG     │
         │             │ sharing (§3.4)  │    │ (slides 18–19)  │
         │             └────────┬────────┘    └────────┬────────┘
         │                      │                      │
         │             ┌────────▼────────┐    ┌────────▼────────┐
         └────────────▶│    sign.rs      │◀───│  partition.rs   │
                       │ K-parallel 3rnd │    │ Alg. 6 tables   │
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

```toml
[dependencies]
threshold-ml-dsa = "0.5"
```

### 1. Fresh Key Generation (Trusted Dealer)

```rust
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use rand::rngs::OsRng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = OsRng;
    let seed = [42u8; 32]; // Use secure seed in production

    // Create a 2-of-3 threshold SDK with fresh key generation
    let mut sdk = ThresholdMlDsa44Sdk::from_seed(&seed, 2, 3, 100)?;

    // Sign message with active parties 0 and 1
    let msg = b"Hello, threshold ML-DSA!";
    let active = [0u8, 1];
    let sig = sdk.threshold_sign(&active, msg, &mut rng)?;

    // Verify with standard ML-DSA-44 verifier
    assert!(sdk.verify(msg, &sig));

    Ok(())
}
```

### 2. Distributed Key Generation (no trusted dealer)

Each party runs its own `DkgParty`. Broadcast messages need an authenticated
channel; `DkgSeedMsg` is secret and must go over an authenticated, encrypted
point-to-point channel to its recipient. No party ever holds a subset secret
it is not a member of.

The design follows the outline in the Mithril NIST MPTS slides (18–19):
per-subset keys `K_I`, a jointly sampled coin, `s_I = H(coin, K_I)`,
commit-and-reveal of `t_I`, `t = Σ t_I`. The paper defers its DKG to the full
version of ePrint 2026/013, which we have not been able to check; the exact
messages and checks here are this crate's own.

```rust
use threshold_ml_dsa::dkg::{DkgConfig, DkgParty};
use threshold_ml_dsa::nonce::InMemoryNonceGuard;

// On party `me`, with `transport` your networking layer:
let (party, r0) = DkgParty::start(DkgConfig { t, n, me, context }, &mut rng)?;
let all_r0 = transport.broadcast(r0);
let (party, coin, seed_msgs) = party.round1(&all_r0)?;
let all_coins = transport.broadcast(coin);
let my_seeds = transport.send_private(seed_msgs); // encrypted, per recipient
let (party, commit) = party.round2(&all_coins, my_seeds, &mut rng)?;
let all_commits = transport.broadcast(commit);
let (party, reveal) = party.round3(&all_commits)?;
let all_reveals = transport.broadcast(reveal);
let (party, confirm) = party.round4(&all_reveals)?;
let all_confirms = transport.broadcast(confirm);
let (pk, my_key) = party.finalize(&all_confirms)?;
let stored = my_key.to_bytes()?; // protect at rest
```

`dkg::simulate` runs all parties in one process for tests.

### 3. Sharing an existing key (a posteriori)

```rust
// On the key owner's machine; the public key stays the same.
let mut sdk = ThresholdMlDsa44Sdk::from_existing_key(&sk_bytes, t, n, 100, &mut rng)?;
// or, to distribute keys: aposteriori::share_existing_key(&sk_bytes, &params, &mut rng)
```

This follows §3.4 of the paper: for every coefficient, the S subset shares
are drawn jointly from a discrete Gaussian (σ = √2, the variance of fresh
keygen) conditioned on their sum being the original coefficient. The coset
Gaussian is sampled exactly and in constant time by rejection, instead of the
paper's GPV sampler. Signing works unchanged.

A short sharing of a fixed key always leaks a hint to T−1 parties. The paper
estimates the loss at **up to 12 bits of security** for σ = √2 (heuristic; its
formal 7-bit bound needs wider shares and parameters from the full version).
Prefer fresh keys when the public key need not be kept, deliver each party's
key over an encrypted channel, and delete the original secret key.

**Share a key at most once.** Every sharing gives a sub-threshold coalition
another independent hint on the same secret, so the loss grows with each
re-sharing (to refresh shares, change (T, N) or replace a party, generate a
new key with the DKG instead).

### Nonce safety

`sign::round1` takes a `NonceGuard`. It refuses a repeated session ID per key
and mixes a strictly increasing counter into the nonce seed, so a repeated
RNG output alone cannot reuse a nonce (which would reveal the party's partial
secret). `InMemoryNonceGuard` is not durable; persist its counter or implement
`NonceGuard` over durable storage. Rolling back the RNG and the guard state
together (e.g. a VM snapshot) is not protected against.

---

## Supported Threshold Configurations

$N \le 6$ uses the 15 parameter sets of ePrint 2026/013 (full version, Fig. 9; values as in the Go reference). $N = 7, 8$ use parameters derived with the reference script `params/hyperball.sage` (ν = 3, per-attempt success ≈ 1/2, K at the −2σ bound of the Monte Carlo estimate) and max-flow balanced partitions as in `params/recover.py`; the paper gives no table for them. Measured per-attempt success: 0.43–0.53.

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
| (2, 7) | 4 | 313,939 | 314,051 | 3.0 |
| (3, 7) | 31 | 285,453 | 285,601 | 3.0 |
| (4, 7) | 322 | 261,941 | 262,110 | 3.0 |
| (5, 7) | 1,071 | 241,018 | 241,166 | 3.0 |
| (6, 7) | 1,008 | 220,521 | 220,633 | 3.0 |
| (7, 7) | 91 | 192,000 | 192,055 | 3.0 |
| (2, 8) | 4 | 313,939 | 314,051 | 3.0 |
| (3, 8) | 58 | 291,228 | 291,406 | 3.0 |
| (4, 8) | 1,195 | 267,837 | 268,048 | 3.0 |
| (8, 8) | 227 | 182,918 | 182,974 | 3.0 |

Cost per party and signing attempt is about 5.25 kB of traffic and 19.3 KB of
memory per unit of K: (4, 8) ≈ 6.3 MB / 23 MB, (5, 7) ≈ 5.6 MB / 21 MB.
(5, 8), (6, 8) and (7, 8) would need K ≈ 14,800, 38,500 and 5,200 (27–198 MB
of traffic, 98–728 MB of memory per attempt) and are not supported.

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
| `partition` | Algorithm 6 (`RSSRecover`): paper tables for $N \le 6$, max-flow balanced partitions for $N = 7, 8$ |
| `rss` | `keygen_from_seed()`: fresh independent secrets per subset (Fig. 5, USENIX version) |
| `aposteriori` | Split an existing ML-DSA-44 secret key with Gaussian shares conditioned on their sum (§3.4) |
| `dkg` | Distributed key generation: per-subset keys, joint coin, commit/reveal of `t_I`, confirm round |
| `nonce` | `NonceGuard` / `InMemoryNonceGuard`: Round 1 nonce-reuse protection |
| `sign` | 3-round party protocol: K-parallel commit, reveal, respond |
| `coordinator` | K-parallel `Combine`: $\delta$-norm check, hint generation, FIPS 204 packing |
| `sdk` | High-level API: `from_seed()`, `from_party_keys()`, `threshold_sign()`, `verify()`, `from_existing_key()` |
| `verify` | Standard FIPS 204 verification via `dilithium-rs` (ML-DSA-44, ML-DSA-65, ML-DSA-87) |
| `error` | `no_std`-compatible error handling |

---

## Testing & Verification

```bash
# Run complete test suite
cargo test --release

# Security regression tests
cargo test --release --test security_regressions

# Run end-to-end threshold signing tests
cargo test --test v03_tests

# N = 7, 8 (large-K configurations are slow; run them in release)
cargo test --release --test n78_tests -- --include-ignored

# Line coverage of the unit tests alone (100%)
cargo llvm-cov --lib --show-missing-lines
```

---

## Security Hardening

- **No Key Reuse**: Fresh keygen samples independent secrets per subset.
- **Single-Use Nonces**: `round2()` consumes `StRound1` and `round3()` consumes `StRound2` by value, so a party reveals each commitment once and each nonce answers at most one challenge.
- **Transcript Binding**: `round2()` checks that the party's own commitment is in its slot; `round3()` rejects non-canonical aggregates and witnesses from another session.
- **Cryptographic Reveal Witness**: `verify_all_round2_reveals()` returns a bound witness token that `round3()` requires, enforcing round ordering.
- **Hedged Nonce Generation**: Round 1 binds PRNG entropy + long-term key + transcript context + party ID + a `NonceGuard` counter, so a repeated RNG output alone does not repeat a nonce. Rolling back the RNG **and** the guard state together (VM snapshots, restored backups) is not protected against.
- **Fail-Closed SDK**: Output signatures are verified by the pure-Rust `dilithium-rs` FIPS 204 verifier before returning to caller.
- **Zeroize-on-Drop (best effort)**: `FVec`, `StRound1`, `StRound2`, `ThresholdPrivateKey`, and `Share` implement `Zeroize` and `Drop`; nonce and secret intermediates are wiped explicitly. `sha3` states are wiped on drop; copies left by moves of other stack values are not.
- **100% Safe Rust**: `#![forbid(unsafe_code)]`.
- **Side channels**: the nonce sampler uses branch-free sine/cosine and rounding instead of `libm`'s branching versions; hardware `f64` arithmetic and `sqrt` remain (see `fvec.rs`).

Report vulnerabilities privately; see [SECURITY.md](SECURITY.md).

---

## References

- [FIPS 204 — ML-DSA](https://csrc.nist.gov/pubs/fips/204/final) — Module-Lattice-Based Digital Signature Standard
- [ePrint 2026/013](https://eprint.iacr.org/2026/013) — Efficient Threshold ML-DSA (Mithril Scheme)
- [`dilithium-rs`](https://crates.io/crates/dilithium-rs) — Pure-Rust FIPS 204 implementation

## License

Licensed under the [MIT License](LICENSE).
