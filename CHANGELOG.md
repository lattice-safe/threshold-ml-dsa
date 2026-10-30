# Changelog

## [0.5.0] - Unreleased

Security fixes from an internal review. This release has breaking API changes.

### Security

- **DKG and a posteriori sharing are now behind the non-default
  `insecure-keygen` feature.** The DKG derived the master seed from
  contributions that every party reveals, so every participant (and any
  observer of the broadcast) could derive the full secret key. A posteriori
  sharing gave the last subset `s − Σ others`, so coalitions smaller than T
  learn `s − s_M` for a short `s_M`, which leaks the original secret.
  Keys created with `from_dkg` or `from_existing_key` should be treated as
  compromised and replaced with keys from `from_seed`.
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
- `round3` checks that the reveal witness has the same session ID as Round 2,
  and fails if the partition assigns a share the party does not hold.
- Secrets are wiped in more places: the Round 1 nonce seed, `c·s_I`
  intermediates and the partial secret in Round 3, and the full secret totals
  in keygen. `Poly`'s `Debug` no longer prints a coefficient.
- Documented that the hedged nonce derivation does not protect against RNG
  rollback or repeated output: a repeated nonce reveals the partial secret.

### Fixed

- `Poly::uniform_eta` and `MatrixA::expand` now match FIPS 204 `ExpandS` and
  `ExpandA` (previously uniform but not FIPS-identical). `Poly::pack_z` centers
  coefficients first. These are not used by the signing protocol.

### Changed (breaking)

- `sign::round2(.., st_rd1: StRound1, ..)` takes the Round 1 state by value.
- `sign::round3(sk, wfinals, st_rd2: StRound2, params, verified)` takes the
  Round 2 state by value and no longer takes `StRound1`.
- `get_threshold_params` returns `None` for N > 6.
- `dkg`, `aposteriori`, `from_dkg` and `from_existing_key` require the
  `insecure-keygen` feature.

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
