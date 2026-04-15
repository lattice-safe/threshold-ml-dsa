# Threshold ML-DSA SDK Guide

## Overview

`threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk` provides a one-call interface for
threshold ML-DSA-44 key generation and signing. It implements the complete
3-round protocol from ePrint 2026/013 internally.

## Quick Start

```rust
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;
use rand::rngs::OsRng;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut rng = OsRng;
    let seed = [42u8; 32]; // use OsRng in production

    // Create a 2-of-3 threshold SDK
    let sdk = ThresholdMlDsa44Sdk::from_seed(
        &seed,
        2,   // T (threshold)
        3,   // N (total parties)
        100, // max full-protocol retries
    )?;

    // Sign — any 2 of 3 parties
    let msg = b"Hello, threshold ML-DSA!";
    let active = [0u8, 1]; // parties 0 and 1
    match sdk.threshold_sign(&active, msg, &mut rng) {
        Ok(sig) => {
            // Verify with any standard ML-DSA-44 verifier
            assert!(sdk.verify(msg, &sig));
        }
        // Fail-closed outcomes (no signature returned)
        Err(threshold_ml_dsa::error::Error::InsufficientResponses)
        | Err(threshold_ml_dsa::error::Error::InvalidSignature) => {}
        Err(e) => return Err(Box::new(e)),
    }
    Ok(())
}
```

## API

### `ThresholdMlDsa44Sdk::from_seed(seed, t, n, max_retries)`

Generate fresh threshold keys from a 32-byte seed. This follows the paper's
Figure 4 keygen: each RSS subset gets an independently sampled secret.

**Parameters:**
- `seed: &[u8; 32]` — deterministic seed (use `OsRng` in production)
- `t: u8` — threshold (minimum parties to sign, 2 ≤ T ≤ N)
- `n: u8` — total parties (2 ≤ N ≤ 6)
- `max_retries: usize` — max full-protocol retry attempts

**Returns:** `Result<Self, Error>`

### `sdk.threshold_sign(active, msg, rng)`

Perform a full 3-round threshold signing protocol with K parallel repetitions.

**Parameters:**
- `active: &[u8]` — sorted list of T active signer IDs
- `msg: &[u8]` — message to sign
- `rng: &mut R` — CSPRNG for commitment randomness

**Returns:** `Result<[u8; 2420], Error>`

- `Ok(sig)`: FIPS 204-valid signature bytes
- `Err(Error::InsufficientResponses | Error::InvalidSignature)`: fail-closed abort
- other `Err(...)`: configuration/protocol failure

### `sdk.verify(msg, sig)`

Verify a signature using the standard ML-DSA-44 verifier.

**Returns:** `bool`

## What Happens Internally

```
from_seed(seed, T, N)
  └── rss::keygen_from_seed()
       ├── Derive ρ, per-party keys from SHAKE-256(seed)
       ├── For each (N-T+1)-subset (Gosper's hack):
       │     sample s_I ← χ_s independently
       ├── Compute s_total = Σ s_I
       └── pk = Pack(ρ, Power2Round(A·s₁ + s₂))

threshold_sign(active, msg, rng)
  ├── Round 1: Each party → K hyperball commitments (SampleHyperball)
  ├── Round 2: Reveal commitments, compute μ = CRH(tr ‖ msg)
  ├── Round 3: Each party → K responses with FVec.Excess rejection
  │     └── Partial secret recovered via partition::rss_recover()
  ├── Combine: Try each of K slots
  │     ├── ‖z_k‖∞ < γ₁ - β
  │     ├── δ = ‖Az - 2^d·c·t₁ - w‖∞ < γ₂
  │     └── MakeHint + pack_sig
  └── Fail-closed: verify(pk, msg, sig) before returning
```

## Supported Configurations

All (T, N) pairs with 2 ≤ T ≤ N ≤ 6 are supported, using exact parameters
from ePrint 2026/013, Figure 8.

## Security Model

- **No key decomposition**: Keys are generated fresh, not split from an existing ML-DSA key.
- **Fail-closed**: `Ok(sig)` is always verified before return; unsuccessful attempts return an error instead of an invalid signature.
- **In-process orchestrator**: The SDK runs all parties in-process. For distributed
  deployment, use the low-level `sign` + `coordinator` modules and add authenticated
  transport between isolated party processes.

## Migration from v0.2

The v0.3 SDK API is a **breaking change**:

| v0.2 | v0.3 |
|------|------|
| `MlDsa44ThresholdMaterial::from_secret_key(pk, sk)` | **Removed** (insecure) |
| `ThresholdMlDsa44Sdk::new(material, n, t, retries, rng)` | `ThresholdMlDsa44Sdk::from_seed(seed, t, n, retries)` |
| `ThresholdMlDsa44Sdk::from_seed(seed, n, t, retries, rng)` | `ThresholdMlDsa44Sdk::from_seed(seed, t, n, retries)` |
| `sdk.sign(msg, rng)` | `sdk.threshold_sign(active, msg, rng)` |
