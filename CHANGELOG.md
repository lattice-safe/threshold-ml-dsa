# Changelog

## [0.4.1] - Unreleased

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
