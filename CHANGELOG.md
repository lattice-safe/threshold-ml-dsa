# Changelog

## [0.5.0] - Unreleased

Security fixes from an internal review. This release has breaking API changes.

### Security

- **New DKG.** The old DKG derived the master seed from contributions that
  every party reveals, so every participant (and any observer of the
  broadcast) could derive the full secret key. `dkg` now follows the outline
  in the Mithril slides (NIST MPTS 2026, slides 18–19): each subset key `K_I`
  is chosen by one of its members and sent privately to the other members,
  all `K_I` and coin shares are committed before a joint coin is revealed,
  `s_I = H(sid, coin, I, K_I)`, `t_I` is committed and revealed, and a
  confirmation round catches equivocation. No party receives a subset secret
  it is not a member of, a dealer cannot grind `s_I`, and no party can bias
  `ρ` except by aborting. The paper's own DKG is in the full version, which
  has not been checked. `dkg_round1/2`,
  `dkg_finalize`, `DkgContribution` and `from_dkg` are removed. **Keys created
  with the old `from_dkg` must be treated as compromised.**
- **A posteriori sharing is behind the non-default `insecure-keygen` feature**
  and renamed `insecure_share_existing_key` / `insecure_from_existing_key`.
  It samples uniform shares plus a non-short remainder, which is not the
  paper's algorithm (§3.4 of the USENIX version: Gaussian shares conditioned
  on their sum, 7–12 bits of security loss), and coalitions smaller than T
  learn part of the original secret. **Keys shared with `from_existing_key` must be
  treated as compromised.** It now also rejects malformed secret keys
  (coefficients outside `[-η, η]`, mismatched `tr`).
- **Nonce-reuse protection.** `round1` takes a `NonceGuard` that refuses a
  repeated session ID per key and mixes a strictly increasing counter into the
  nonce seed. Previously a repeated RNG output let a coordinator replay a
  session and recover the party's partial secret from two responses.
- **N = 7, 8 removed** (`MAX_PARTIES` is now 6). Their parameters were not
  derived from the paper's formulas, and their radii are 0.46–0.90 of what the
  rejection-sampling bound requires.
- **`round2` consumes `StRound1` by value.** It previously borrowed it, so a
  party could be made to reveal again against a transcript whose peer
  commitments were chosen after its own reveal was seen.
- **`round2` checks that the party's own commitment is in its slot** and that
  the transcript has one commitment per active party. A coordinator could
  previously replace a party's commitment and control the challenge.
- **`round3` rejects non-canonical aggregate commitments.** `pack_w_single` is
  not injective outside `[0, q)`, so a coordinator could send an aggregate that
  passes the hash binding but decomposes to a different challenge.
- **Parameters are validated.** `ThresholdParams` fields are public; a radius
  of NaN or ∞ disabled rejection and released every response. `round1/2/3`,
  `combine` and a posteriori sharing now require the canonical parameters for
  (T, N), and keys record their (T, N) and refuse other parameters.
- `round3` checks that the reveal witness has the same session ID as Round 2,
  and fails if the partition assigns a share the party does not hold.
- `verify*` return `false` for a wrong-length public key instead of
  panicking; `sign_standard*` check the secret key length and return code.
- The coordinator rejects party commitments and responses outside `[0, q)`
  (previously i32 overflow), and `Signature::to_bytes` checks field lengths
  instead of panicking.
- Secrets are wiped in more places: Round 1 nonce seed and integer nonce,
  Gaussian buffers in `sample_hyperball`, `c·s_I` intermediates and the
  partial secret in Round 3, keygen totals and a posteriori intermediates.
  `ThresholdPrivateKey` is no longer `Clone`; `keygen*` return the secret key
  in `Zeroizing`. `Poly`'s `Debug` no longer prints a coefficient.

### Added

- `ThresholdPrivateKey::{to_bytes, from_bytes}`: versioned, validated key
  serialization so each party can hold only its own key.
- `ThresholdMlDsa44Sdk::from_party_keys` to build an SDK from DKG or parsed keys.
- `nonce::{NonceGuard, InMemoryNonceGuard}`, `ThresholdPrivateKey::guard_key_id`.
- `dkg::simulate` to run a whole DKG in one process (tests, benchmarks).
- `ThresholdParams::validate`, `ThresholdPrivateKey::threshold`.

### Fixed

- `keygen_from_seed` now absorbs `[K, L]` after the seed, as the Go reference
  and FIPS 204 do; keys match the Go reference byte for byte (test vectors in
  `tests/security_regressions.rs`). **Keys derived from a given seed differ
  from 0.4.x.**
- `Poly::uniform_eta` and `MatrixA::expand` now match FIPS 204 `ExpandS` and
  `ExpandA` (previously uniform but not FIPS-identical). `Poly::pack_z` centers
  coefficients first. These are not used by the signing protocol.
- Documentation no longer claims the DKG and a posteriori sharing implement
  ePrint 2026/013 §4.3/§4.2, or "production-grade" status. Paper references
  were checked against the USENIX Security '26 version: keygen is Fig. 5, a
  posteriori sharing is §3.4, and the DKG and per-(T, N) parameters are in the
  full version.

### Changed (breaking)

- `sign::round1(.., rng, guard)` takes a `NonceGuard`.
- `sign::round2(.., st_rd1: StRound1, ..)` takes the Round 1 state by value.
- `sign::round3(sk, wfinals, st_rd2: StRound2, params, verified)` takes the
  Round 2 state by value and no longer takes `StRound1`.
- `ThresholdMlDsa44Sdk::threshold_sign` takes `&mut self`.
- `get_threshold_params` returns `None` for N > 6.
- `Error` is `#[non_exhaustive]` and has new variants.
- `verify::sign_standard*` return `Result`; `verify::keygen*` return the
  secret key as `Zeroizing<Vec<u8>>`; `Signature::to_bytes` returns `Result`.
- Round 1 commitments and nonces use new domain tags (`-v2`).

## [0.4.1] - 2026-10-31

### Security

- Fix hyperball sampler norm (`src/fvec.rs`, `sample_hyperball`). The squared
  norm was accumulated after ν-scaling the first L·N coordinates, so samples had
  ν-metric norm ≈0.447·r₁ instead of ≈r₁. Round-3 hyperball rejection therefore
  almost never triggered, and accepted responses were not independent of the
  secret shares as the scheme's security argument requires. Affects 0.3.1–0.4.0
  (introduced in d429a0b); 0.3.0 is unaffected. All threshold-signing users
  should upgrade, and keys used for threshold signing with affected versions
  should be considered for rotation.
  Thanks to Md Alquamah for responsible disclosure.

### Changed

- Signing now rejects at the rate given in ePrint 2026/013, so expect more
  retries per signature than with 0.3.1–0.4.0.

### Tests

- Hyperball sampler tests now also assert a lower bound (≥ 0.99·radius) on the
  ν-metric norm.
- New test checking that per-attempt signing success at (2, 2) with
  `max_retries = 1` is well below 1.
