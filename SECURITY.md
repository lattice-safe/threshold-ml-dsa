# Security Policy

## Supported versions

| Version | Supported |
|---|---|
| 0.5.x | Yes |
| < 0.5.0 | No: affected by the advisories below; 0.3.1–0.4.1 are yanked |

## Reporting a vulnerability

Please report vulnerabilities privately, not in public issues:

- GitHub: **Security → Report a vulnerability** on this repository, or
- e-mail: suleyman@lattice-safe.com

Include the affected version, a description, and if possible a reproduction.
We aim to acknowledge reports within a few days, agree on a disclosure date
with the reporter, and credit reporters in the advisory unless they prefer
otherwise.

## Published advisories

| Advisory | Severity | Affected | Fixed | Summary |
|---|---|---|---|---|
| [GHSA-hqr3-xgxm-59jr](https://github.com/lattice-safe/threshold-ml-dsa/security/advisories/GHSA-hqr3-xgxm-59jr) | Critical | < 0.5.0 | 0.5.0 | DKG revealed the full secret key; insecure a posteriori sharing; N = 7, 8 parameters; nonce reuse and transcript-binding issues |
| [GHSA-22vh-j6fg-hm7g](https://github.com/lattice-safe/threshold-ml-dsa/security/advisories/GHSA-22vh-j6fg-hm7g) | High | 0.3.1 – 0.4.0 | 0.4.1 | Hyperball sampler norm bug disabled rejection sampling |

## Key rotation after upgrading

Replace keys that were:

- created with `from_dkg` (0.4.x) or shared with `from_existing_key` /
  `from_secret_key` before 0.5.0;
- used for threshold signing with 0.3.1 – 0.4.0, or with N = 7 or 8 before
  0.5.0.

## Scope and assumptions

See the "Security Hardening" section of the [README](README.md) and
[PROTOCOL.md](PROTOCOL.md) for the threat model, including the remaining
assumptions: authenticated broadcast and encrypted point-to-point channels for
the DKG, a durable `NonceGuard` across restarts, and the 7–12-bit security loss
of a posteriori sharing (share a key at most once).
