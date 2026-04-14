# Threshold ML-DSA SDK Guide

This project provides two layers:

- **Low-level protocol modules** (`rss`, `sign`, `coordinator`) for custom orchestration
- **High-level SDK layer** (`sdk`) for application integration

## What The SDK Gives You

`threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk` wraps:

1. key material extraction for thresholdization
2. RSS distribution and party construction
3. fail-closed threshold signing
4. FIPS 204 verification helper

## Minimal Example

```rust
use rand::rngs::OsRng;
use threshold_ml_dsa::sdk::ThresholdMlDsa44Sdk;

let mut rng = OsRng;
let seed = [42u8; 32];

let (mut sdk, _sk) = ThresholdMlDsa44Sdk::from_seed(&seed, 4, 3, 128, &mut rng)?;
let msg = b"sdk example";
let sig = sdk.sign(msg, &mut rng)?;
assert!(sdk.verify(msg, &sig));
```

## Security Model Notes

- Signing is **fail-closed**: no unverifiable signature is returned as success.
- SDK is an **in-process orchestrator**; it does not provide a network protocol.
- For real distributed deployment:
  - isolate parties into separate processes/devices
  - use authenticated/encrypted message transport
  - enforce replay protection at transport/session layer

See `tests/encrypted_transport_tests.rs` for an example transport harness using
ML-KEM (`lattice-kyber`) + AEAD.
