# threshold-ml-dsa

**Threshold ML-DSA (FIPS 204) via Replicated Secret Sharing**

A `#![no_std]`-compatible Rust implementation of threshold ML-DSA-44 based on the Mithril scheme ([ePrint 2026/013](https://eprint.iacr.org/2026/013)). On successful aggregation, produced signatures are compatible with standard FIPS 204 verifiers. The coordinator is fail-closed: if a valid aggregate is not obtained, signing returns an error.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org)

---

## Overview

Standard threshold signature schemes based on Shamir secret sharing introduce large Lagrange interpolation multipliers that blow up lattice coefficient sizes, breaking the short-vector requirements of ML-DSA. The Mithril scheme solves this with **Replicated Secret Sharing (RSS)**, which keeps all shares short by construction.

### Key Properties

- **FIPS 204 Compatible on Success** — Returned signatures pass `dilithium-rs` verification unchanged
- **Fail-Closed Coordinator** — Invalid aggregates are never returned as successful signatures
- **Local Rejection Sampling** — Each party independently checks an L₂ hyperball bound; the coordinator retries across qualifying signer sets
- **Short Shares** — RSS preserves η-bounded coefficients (no Lagrange blow-up)
- **`#![no_std]`** — Suitable for embedded and TEE environments
- **Memory Safe** — Sensitive key material is `zeroize`d on drop; constant-time comparisons via `subtle`

## Architecture

```
┌─────────────┐     ┌──────────────┐     ┌──────────────┐
│   params.rs  │────▶│   poly.rs    │────▶│   rss.rs     │
│ ML-DSA-44    │     │ NTT, norms,  │     │ Key distrib. │
│ constants    │     │ sampling     │     │ & reconstruct│
└─────────────┘     └──────┬───────┘     └──────┬───────┘
                           │                     │
                    ┌──────▼───────┐     ┌───────▼──────┐
                    │  sign.rs     │◀────│ Party state  │
                    │  Commit/Sign │     │ machine      │
                    └──────┬───────┘     └──────────────┘
                           │
                    ┌──────▼───────┐     ┌──────────────┐
                    │coordinator.rs│────▶│  verify.rs   │
                    │  Aggregation │     │  dilithium-rs│
                    └──────────────┘     └──────────────┘
```

## Protocol (4 Rounds — Hardened)

| Round | Party Action | Coordinator Action |
|-------|-------------|-------------------| 
| **0 — Pre-commit** | Compute H(w_i ‖ party_id) binding hash | Collect all binding hashes |
| **1 — Reveal** | Send full commitment w_i | Verify each w_i against its hash (ADV-1) |
| **2 — Challenge** | — | Aggregate w = Σ w_i, compute c̃ = H(μ ‖ w₁_packed) |
| **3 — Sign** | Compute z_i = c·s₁_i + y_i with hyperball check, attach session binding H(c̃ ‖ party_id) | Verify session bindings (ADV-6), deduplicate party_ids (ADV-2), aggregate z = Σ z_i, compute Hint, encode FIPS 204 signature, verify end-to-end before return |

## Quick Start

```rust
use threshold_ml_dsa::{rss, sign, coordinator, verify, poly::MatrixA, params::*};
use rand::rngs::OsRng;

// 1. Generate an ML-DSA-44 key pair
let seed = [0u8; 32]; // use OsRng in production
let (pk, _sk) = verify::keygen(&seed);

// 2. Parse s₁, s₂ from sk and distribute via RSS
// let shares = rss::distribute_key(&s1, &s2, 3, 2, &mut OsRng)?;

// 3. Create signing parties
// let rho = &pk[..SEEDBYTES];
// let a_hat = MatrixA::expand(rho.try_into().unwrap());
// let parties: Vec<_> = shares.iter()
//     .map(|s| sign::Party::new(s, a_hat.clone()))
//     .collect();

// 4. Run threshold signing protocol (4-round hardened)
// let sig = coordinator::threshold_sign(&mut parties, msg, &tr, &pk, 2, 50, &mut OsRng)?;

// 5. Verify with any standard ML-DSA-44 verifier
// assert!(verify::verify(&sig.to_bytes(), msg, &pk));
```

## Modules

| Module | Description |
|--------|-------------|
| `params` | ML-DSA-44 constants (q, N, K, L, γ₁, γ₂, η, τ, β) and threshold hyperball bound B² |
| `poly` | Polynomial arithmetic over ℤ_q[X]/(X²⁵⁶+1): NTT, norms, SHAKE sampling, packing |
| `rss` | Replicated Secret Sharing: key distribution and reconstruction for (N, T)-threshold |
| `sign` | Party state machine with precommit/reveal/commit/sign rounds, hedged nonces, session binding |
| `coordinator` | Signature aggregation with commitment verification, Sybil protection, and automatic retry |
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
| τ | 39 | Challenge weight |
| B² | L·N·γ₁² | Hyperball rejection bound |
| Max parties | 6 | RSS subset limit |

## Testing

```bash
cargo test
```

### Cross-Validation Against dilithium-rs

Every critical primitive is verified byte-exact against the reference implementation:

| Primitive | Test | Coverage |
|-----------|------|----------|
| `decompose` | `test_decompose_matches_reference` | 2,200+ values |
| `make_hint` | `test_make_hint_matches_reference` | 8,402 values |
| `use_hint` | `test_use_hint_matches_reference` | 16,804 pairs |
| `pack_w1` | `test_w1_packing_matches_reference` | 4 full polynomials |
| `pack_z` | `test_z_packing_matches_reference` | 14 edge cases + roundtrip |
| Full protocol | `test_full_threshold_signature_flow` | E2E sign → FIPS verify |

### Security Test Coverage

- **RSS reconstruction** — Any qualifying T-of-N party set reconstructs; malformed/inconsistent replicas are rejected
- **Hyperball rejection** — Per-party L₂ checks with safe coordinator retries
- **Commitment binding** — Precommit/reveal/verify round
- **Session binding** — Stale z_i rejection across sessions
- **Full signing flow** — 4-round threshold signing with FIPS 204 verification

## Security Notes

- The coordinator does **not** reconstruct a full ML-DSA secret key during signing.
- Signing is fail-closed: if retries do not produce a verifiable aggregate, an error is returned.
- RSS reconstruction validates subset-index bounds and replica consistency.

## Dependencies

| Crate | Purpose |
|-------|---------|
| `dilithium-rs` | FIPS 204 NTT, Montgomery reduction, and canonical verification |
| `sha3` | SHAKE-128/256 for deterministic sampling |
| `zeroize` | Secure erasure of sensitive key material |
| `subtle` | Constant-time operations |
| `rand_core` | CSPRNG trait for key generation |

## References

- [FIPS 204 — ML-DSA (Module-Lattice-Based Digital Signature Standard)](https://csrc.nist.gov/pubs/fips/204/final)
- [ePrint 2026/013 — Mithril: Threshold ML-DSA via Replicated Secret Sharing](https://eprint.iacr.org/2026/013)
- [`dilithium-rs`](https://crates.io/crates/dilithium-rs) — Pure-Rust FIPS 204 implementation

## License

This project is licensed under the [MIT License](LICENSE).
