# Threshold ML-DSA Protocol Specification

> **Based on:** FIPS 204 (ML-DSA) + ePrint 2026/013 (Mithril Scheme)
>
> This document describes the complete cryptographic protocol implemented in
> `threshold-ml-dsa`, covering key generation, key distribution via Replicated
> Secret Sharing, the hardened 4-round threshold signing protocol, adversarial
> defenses, and verification.

---

## Table of Contents

1. [Notation and Parameters](#1-notation-and-parameters)
2. [Algebraic Foundation](#2-algebraic-foundation)
3. [Standard ML-DSA (FIPS 204)](#3-standard-ml-dsa-fips-204)
4. [Why Threshold Lattice Signing Is Hard](#4-why-threshold-lattice-signing-is-hard)
5. [Replicated Secret Sharing (RSS)](#5-replicated-secret-sharing-rss)
6. [Threshold Key Generation](#6-threshold-key-generation)
7. [Threshold Signing Protocol (4 Rounds)](#7-threshold-signing-protocol-4-rounds)
8. [Hyperball Rejection Sampling](#8-hyperball-rejection-sampling)
9. [Signature Aggregation and Encoding](#9-signature-aggregation-and-encoding)
10. [Verification](#10-verification)
11. [Adversarial Security Hardening](#11-adversarial-security-hardening)
12. [Sequence Diagrams](#12-sequence-diagrams)
13. [Activity Diagrams](#13-activity-diagrams)
14. [Security Properties](#14-security-properties)

---

## 1. Notation and Parameters

### Symbols

| Symbol | Meaning |
|--------|---------|
| $\mathbb{Z}_q$ | Integers modulo $q$ |
| $R_q$ | Ring $\mathbb{Z}_q[X]/(X^{256}+1)$ |
| $\mathbf{s}_1 \in R_q^\ell$ | Secret key vector (short, $\eta$-bounded) |
| $\mathbf{s}_2 \in R_q^k$ | Secret key vector (short, $\eta$-bounded) |
| $\mathbf{A} \in R_q^{k \times \ell}$ | Public matrix (expanded from seed $\rho$) |
| $\mathbf{t} = \mathbf{A}\mathbf{s}_1 + \mathbf{s}_2$ | Public key vector |
| $c \in R_q$ | Challenge polynomial ($\tau$ non-zero $\pm 1$ coefficients) |
| $\mathbf{y} \in R_q^\ell$ | Masking vector ($\gamma_1$-bounded) |
| $\mathbf{z} = c \cdot \mathbf{s}_1 + \mathbf{y}$ | Signature response vector |
| $\lVert \cdot \rVert_2$ | Euclidean (L₂) norm |
| $\lVert \cdot \rVert_\infty$ | Infinity (L∞) norm |

### ML-DSA-44 Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Modulus | $q$ | $8{,}380{,}417 = 2^{23} - 2^{13} + 1$ | Prime for $\mathbb{Z}_q$ |
| Ring degree | $n$ | $256$ | Polynomial degree |
| Matrix rows | $k$ | $4$ | Rows in **A** |
| Matrix cols | $\ell$ | $4$ | Columns in **A** |
| Secret bound | $\eta$ | $2$ | $\mathbf{s}_1, \mathbf{s}_2 \in [-\eta, \eta]^n$ |
| Challenge weight | $\tau$ | $39$ | Non-zero entries in $c$ |
| Masking bound | $\gamma_1$ | $2^{17} = 131{,}072$ | $\mathbf{y} \in [-\gamma_1+1, \gamma_1]^n$ |
| Decompose param | $\gamma_2$ | $(q-1)/88 = 95{,}232$ | For HighBits/LowBits |
| Norm bound | $\beta$ | $\tau \cdot \eta = 78$ | $\lVert c \cdot \mathbf{s} \rVert_\infty \leq \beta$ |
| Rounding bits | $d$ | $13$ | Power2Round parameter |
| Hint budget | $\omega$ | $80$ | Max hint ones |

### Threshold-Specific Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Max parties | $N_{\max}$ | $6$ | Upper bound on $N$ |
| Hyperball bound | $B^2$ | $\ell \cdot n \cdot \gamma_1^2$ | Per-party L₂ rejection |
|  |  | $= 17{,}592{,}186{,}044{,}416$ | $\approx 1.76 \times 10^{13}$ |

---

## 2. Algebraic Foundation

### The Ring $R_q$

All arithmetic takes place in the polynomial ring:

$$R_q = \mathbb{Z}_q[X] / (X^{256} + 1)$$

Elements are polynomials of degree $< 256$ with coefficients in $\mathbb{Z}_q$. Multiplication is polynomial multiplication modulo $X^{256} + 1$.

### Number Theoretic Transform (NTT)

Polynomial multiplication is accelerated via the NTT. For $a, b \in R_q$:

$$a \cdot b = \text{NTT}^{-1}(\text{NTT}(a) \circ \text{NTT}(b))$$

where $\circ$ denotes coefficient-wise (pointwise) multiplication in the NTT domain using Montgomery reduction:

$$\text{PointwiseMont}(\hat{a}_i, \hat{b}_i) = \hat{a}_i \cdot \hat{b}_i \cdot R^{-1} \bmod q$$

with Montgomery constant $R = 2^{32}$.

### Norms

For a polynomial $p \in R_q$ with centered coefficients $p_0, \ldots, p_{255}$:

$$\lVert p \rVert_2^2 = \sum_{i=0}^{255} p_i^2 \qquad \lVert p \rVert_\infty = \max_{i} |p_i|$$

For a vector $\mathbf{v} = (v_0, \ldots, v_{\ell-1}) \in R_q^\ell$:

$$\lVert \mathbf{v} \rVert_2^2 = \sum_{j=0}^{\ell-1} \lVert v_j \rVert_2^2 \qquad \lVert \mathbf{v} \rVert_\infty = \max_{j} \lVert v_j \rVert_\infty$$

### Centering

Coefficients are stored in $[0, q)$ but centered to $[-(q-1)/2, (q-1)/2]$ for norm computations:

$$\text{center}(a) = \begin{cases} a & \text{if } a \leq (q-1)/2 \\ a - q & \text{if } a > (q-1)/2 \end{cases}$$

---

## 3. Standard ML-DSA (FIPS 204)

### 3.1 Key Generation

```
KeyGen(seed):
    (ρ, ρ', K) ← SHAKE-256(seed)
    A ← ExpandA(ρ)                        // A ∈ R_q^{k×ℓ}
    s₁ ← SampleShort(ρ', η, 0..ℓ)        // s₁ ∈ R_q^ℓ,  coeffs ∈ [-η, η]
    s₂ ← SampleShort(ρ', η, ℓ..ℓ+k)     // s₂ ∈ R_q^k,  coeffs ∈ [-η, η]
    t = A·s₁ + s₂                         // t ∈ R_q^k
    (t₁, t₀) = Power2Round(t, d)
    pk = (ρ, t₁)
    tr = SHAKE-256(pk)
    sk = (ρ, K, tr, s₁, s₂, t₀)
    return (pk, sk)
```

**Math:**

$$\mathbf{t} = \mathbf{A} \cdot \mathbf{s}_1 + \mathbf{s}_2 \in R_q^k$$

$$\mathbf{t} = \mathbf{t}_1 \cdot 2^d + \mathbf{t}_0 \quad \text{(Power2Round)}$$

### 3.2 Standard Signing

```
Sign(sk, msg):
    μ = SHAKE-256(tr ‖ msg)
    ρ' = SHAKE-256(K ‖ rnd ‖ μ)
    κ ← 0
    loop:
        y ← ExpandMask(ρ', κ)             // y ∈ S_{γ₁}^ℓ
        w = A · NTT(y)
        w₁ = HighBits(w)
        c̃ = SHAKE-256(μ ‖ w₁_packed)
        c = SampleInBall(c̃)               // τ non-zero ±1 coefficients
        z = y + c · s₁
        if ‖z‖∞ ≥ γ₁ - β: κ += ℓ; continue
        r₀ = LowBits(w - c·s₂)
        if ‖r₀‖∞ ≥ γ₂ - β: κ += ℓ; continue
        h = MakeHint(...)
        if |h| > ω: κ += ℓ; continue
        return σ = (c̃, z, h)
```

**Core equation (Fiat-Shamir with aborts):**

$$\mathbf{z} = \mathbf{y} + c \cdot \mathbf{s}_1$$

**Rejection conditions:**

$$\lVert \mathbf{z} \rVert_\infty < \gamma_1 - \beta$$

### 3.3 Verification

```
Verify(pk, msg, σ = (c̃, z, h)):
    if ‖z‖∞ ≥ γ₁ - β: return ⊥
    A ← ExpandA(ρ)
    μ = SHAKE-256(tr ‖ msg)
    c = SampleInBall(c̃)
    w₁' = UseHint(h, A·z - c·t₁·2^d)
    c̃' = SHAKE-256(μ ‖ w₁'_packed)
    return c̃ == c̃'    (constant-time)
```

**Verification equation:**

$$\mathbf{A} \cdot \mathbf{z} - c \cdot \mathbf{t}_1 \cdot 2^d \equiv \mathbf{A}\mathbf{y} + c(\mathbf{A}\mathbf{s}_1 + \mathbf{s}_2) - c \cdot \mathbf{t} \cdot 2^d$$

---

## 4. Why Threshold Lattice Signing Is Hard

### The Shamir Problem

In Shamir $(T, N)$-sharing, reconstruction uses Lagrange coefficients $\lambda_i$:

$$s = \sum_{i \in \mathcal{Q}} \lambda_i \cdot s_i$$

For lattice schemes, this is catastrophic because $\lambda_i$ can be large integers.
If $\lVert s \rVert_\infty \leq \eta$ and $|\lambda_i| \gg 1$, then:

$$\lVert \lambda_i \cdot s_i \rVert_\infty = |\lambda_i| \cdot \lVert s_i \rVert_\infty \gg \eta$$

This **blows up** the coefficient magnitudes, making the response $z_i = y_i + c \cdot \lambda_i \cdot s_i$ far too large to pass the FIPS 204 norm checks.

### The Mithril Solution

**Replicated Secret Sharing (RSS)** avoids Lagrange multipliers entirely:

- Secret is split into $\binom{N}{N-T+1}$ additive pieces — one per $(N-T+1)$-subset
- Each piece is **short by construction** (bounded by $\eta$)
- Reconstruction is **pure addition** — no multipliers needed
- Each party's local holding is bounded by $\binom{N-1}{N-T} \cdot \eta$

---

## 5. Replicated Secret Sharing (RSS)

### 5.1 Subset Structure

For $(N, T)$-threshold, define $M = N-T+1$ and enumerate all $M$-subsets of $\{0, \ldots, N-1\}$:

$$\mathcal{S}_M = \{ S \subseteq \{0,\ldots,N-1\} : |S| = M \}$$

$$|\mathcal{S}_M| = \binom{N}{M} = \binom{N}{N-T+1}$$

**Example** $(N=3, T=2)$:

$$M=2,\ \mathcal{S}_2 = \big\{ \{0,1\},\ \{0,2\},\ \{1,2\} \big\} \qquad |\mathcal{S}_2| = 3$$

### 5.2 Share Distribution

Given secret $\mathbf{s} \in R_q^\ell$ with $\lVert \mathbf{s} \rVert_\infty \leq \eta$:

$$\text{For each } S \in \mathcal{S}_M \text{ except the last: sample } \mathbf{s}^{(S)} \xleftarrow{\$} [-\eta, \eta]^{n \cdot \ell}$$

$$\text{Last share: } \mathbf{s}^{(S_{\text{last}})} = \mathbf{s} - \sum_{S \neq S_{\text{last}}} \mathbf{s}^{(S)}$$

**Additive reconstruction property:**

$$\sum_{S \in \mathcal{S}_M} \mathbf{s}^{(S)} = \mathbf{s}$$

**Party $i$ receives:** all shares $\mathbf{s}^{(S)}$ for subsets $S \ni i$.

### 5.3 Per-Party Effective Share

Party $i$'s effective secret share is:

$$\mathbf{s}_{1,i}^{\text{hold}} = \sum_{\substack{S \in \mathcal{S}_M \\ i \in S}} \mathbf{s}_1^{(S)}$$

Party $i$ participates in $\binom{N-1}{M-1}=\binom{N-1}{N-T}$ subsets, so:

$$\lVert \mathbf{s}_{1,i}^{\text{hold}} \rVert_\infty \leq \binom{N-1}{N-T} \cdot \eta$$

**Example** $(N=3, T=2)$: Each party holds $\binom{2}{1} = 2$ shares, so $\lVert \mathbf{s}^{\text{hold}}_{1,i} \rVert_\infty \leq 2 \cdot 2 = 4$.

### 5.4 Reconstruction

Any $T$ qualifying parties collectively cover **all** $\binom{N}{N-T+1}$ subsets. Each share piece is replicated across all members of its subset, so any one member can provide it:

$$\mathbf{s}_1 = \sum_{S \in \mathcal{S}_M} \mathbf{s}_1^{(S)} = \sum_{\text{one copy per subset from qualifying parties}}$$

```mermaid
graph LR
    subgraph "Secret s₁"
        S["s₁"]
    end
    subgraph "RSS Shares (N=3, T=2)"
        S01["s⁽⁰'¹⁾"]
        S02["s⁽⁰'²⁾"]
        S12["s⁽¹'²⁾"]
    end
    subgraph "Party Holdings"
        P0["Party 0: s⁽⁰'¹⁾ + s⁽⁰'²⁾"]
        P1["Party 1: s⁽⁰'¹⁾ + s⁽¹'²⁾"]
        P2["Party 2: s⁽⁰'²⁾ + s⁽¹'²⁾"]
    end
    S --> S01 & S02 & S12
    S01 --> P0 & P1
    S02 --> P0 & P2
    S12 --> P1 & P2
```

---

## 6. Threshold Key Generation

```mermaid
flowchart TD
    A["Trusted Dealer generates seed"] --> B["KeyGen(seed) → pk, sk"]
    B --> C["Parse sk → ρ, K, tr, s₁, s₂, t₀"]
    C --> D["RSS.distribute(s₁, s₂, N, T)"]
    D --> E1["Party 0: share₀"]
    D --> E2["Party 1: share₁"]
    D --> E3["Party 2: share₂"]
    E1 & E2 & E3 --> F["Each party stores replicated share pieces"]
    F --> G["Per-sign session derives assigned signing share"]
    
    style A fill:#ff6b6b,color:#fff
    style B fill:#ffa94d,color:#fff
    style D fill:#51cf66,color:#fff
    style G fill:#339af0,color:#fff
```

### Step-by-Step

**Step 1 — Standard Key Generation:**

$$(\rho, \rho', K) \leftarrow \text{SHAKE-256}(\text{seed})$$
$$\mathbf{A} \leftarrow \text{ExpandA}(\rho)$$
$$\mathbf{s}_1 \leftarrow \text{CBD}_\eta(\rho'), \quad \mathbf{s}_2 \leftarrow \text{CBD}_\eta(\rho')$$
$$\mathbf{t} = \mathbf{A}\mathbf{s}_1 + \mathbf{s}_2$$
$$pk = (\rho, \mathbf{t}_1), \quad tr = H(pk)$$

**Step 2 — RSS Distribution:**

$$\text{For each } S \in \mathcal{S}_M \setminus \{S_{\text{last}}\}: \quad \mathbf{s}_1^{(S)} \xleftarrow{\$} [-\eta, \eta]^{n \cdot \ell}$$
$$\mathbf{s}_1^{(S_\text{last})} = \mathbf{s}_1 - \sum_{S \neq S_\text{last}} \mathbf{s}_1^{(S)}$$

$$\text{Identically for } \mathbf{s}_2$$

**Step 3 — Party Initialization:**

$$\mathbf{s}^{\text{hold}}_{1,i} = \sum_{\substack{S \ni i}} \mathbf{s}_1^{(S)}, \qquad \mathbf{s}^{\text{hold}}_{2,i} = \sum_{\substack{S \ni i}} \mathbf{s}_2^{(S)}$$
$$\hat{\mathbf{A}} = \text{NTT}(\mathbf{A})$$

---

## 7. Threshold Signing Protocol (4 Rounds)

The protocol runs in **4 rounds** between a coordinator and an active signer set $\mathcal{A}$ of size $T$ (chosen from available parties). The implementation retries over qualifying signer sets and fails closed if no verifiable aggregate is found.

### 7.0 Round 0 — Pre-Commit (ADV-1: Commitment Binding)

Each party $P_i$ independently:

**Hedged nonce derivation (ADV-4):**

$$\text{seed}_i = \text{SHAKE-256}(\text{rng\_entropy} \| \mathbf{s}^{\text{sign}}_{1,i} \| \text{session\_context} \| i \| \text{counter}_i)$$

This binds the masking seed to fresh entropy, signing-share material, session context, and a local counter to mitigate RNG rollback/reuse.

**Sample masking vector:**

$$\mathbf{y}_i \leftarrow \text{ExpandMask}(\text{seed}_i) \in S_{\gamma_1}^\ell$$

**Compute partial commitment:**

$$\hat{\mathbf{y}}_i = \text{NTT}(\mathbf{y}_i)$$
$$\mathbf{w}_i = \text{NTT}^{-1}(\hat{\mathbf{A}} \circ \hat{\mathbf{y}}_i) \in R_q^k$$

**Compute binding hash:**

$$h_i = \text{SHAKE-256}(\mathbf{w}_i \| i)$$

**Send** $h_i$ to coordinator. **Store** $\mathbf{y}_i$ and $\mathbf{w}_i$ locally.

### 7.1 Round 1 — Reveal Commitments

After all binding hashes are collected:

**Send** full commitment $\mathbf{w}_i$ to coordinator.

**Coordinator verifies** each commitment against its pre-committed hash:

$$\text{SHAKE-256}(\mathbf{w}_i \| i) \stackrel{?}{=} h_i \quad \text{(constant-time)}$$

This prevents the coordinator from selectively including/excluding parties to grind for a favorable challenge.

### 7.2 Round 2 — Challenge (Coordinator)

The coordinator computes:

**Aggregate commitment:**

$$\mathbf{w} = \sum_{i \in \mathcal{A}} \mathbf{w}_i$$

**Extract high bits and pack (reference-faithful layout):**

$$\mathbf{w}_1 = \text{HighBits}(\mathbf{w}, 2\gamma_2)$$

$$\text{w1\_packed} = \text{polyw1\_pack}(\mathbf{w}_1)$$

**Fiat-Shamir challenge:**

$$\mu = \text{SHAKE-256}(tr \| \text{msg})$$

$$\tilde{c} = \text{SHAKE-256}(\mu \| \text{w1\_packed})$$

$$c = \text{SampleInBall}(\tilde{c}) \in R_q \quad \text{(exactly } \tau \text{ coefficients} \in \{-1, +1\}\text{)}$$

**Broadcast** $\tilde{c}$ to all parties.

### 7.3 Round 3 — Sign (with Local Rejection + Session Binding)

Each party $P_i \in \mathcal{A}$:

**Decode challenge:**

$$c = \text{SampleInBall}(\tilde{c})$$

**Compute partial response:**

$$\hat{c} = \text{NTT}(c)$$
$$\hat{\mathbf{s}}^{\text{sign}}_{1,i} = \text{NTT}(\mathbf{s}^{\text{sign}}_{1,i})$$
$$\text{For } j = 0, \ldots, \ell-1:$$
$$\quad (c \cdot \mathbf{s}^{\text{sign}}_{1,i})_j = \text{NTT}^{-1}(\hat{c} \circ \hat{\mathbf{s}}^{\text{sign}}_{1,i,j})$$
$$\quad \mathbf{z}_{i,j} = (c \cdot \mathbf{s}^{\text{sign}}_{1,i})_j + \mathbf{y}_{i,j}$$

**Hyperball rejection check:**

$$\lVert \mathbf{z}_i \rVert_2^2 = \sum_{j=0}^{\ell-1} \sum_{m=0}^{n-1} \text{center}(\mathbf{z}_{i,j,m})^2$$

$$\text{If } \lVert \mathbf{z}_i \rVert_2^2 > B^2: \quad \text{REJECT (abort locally)}$$

**Compute session binding (ADV-6):**

$$\text{binding}_i = \text{SHAKE-256}(\tilde{c} \| i)$$

**If accepted:** Send $(\mathbf{z}_i, \text{binding}_i)$ to coordinator.

> [!IMPORTANT]
> The L∞ check $\lVert \mathbf{z}_i \rVert_\infty < \gamma_1 - \beta$ is **not** applied per-party.
> It applies only to the **aggregated** $\mathbf{z}$.
>
> Implementation note: each RSS subset-piece is deterministically assigned to one
> active signer (smallest active party id in that subset), producing
> $\mathbf{s}^{\text{sign}}_{1,i}$ and preventing replicated-piece double counting.

---

## 8. Hyperball Rejection Sampling

### The Exponential Degradation Problem

Standard ML-DSA uses **hypercube** rejection: accept iff $\lVert \mathbf{z} \rVert_\infty < \gamma_1 - \beta$.

In a threshold setting with $N$ parties, each performing independent L∞-rejection, the combined acceptance probability is:

$$p_{\text{cube}} = \left(\frac{\gamma_1 - \beta}{\gamma_1}\right)^{N \cdot \ell \cdot n}$$

This decays **exponentially** in $N$ — for $N = 3$, $p_{\text{cube}}$ is negligible.

### The Hyperball Solution

Replace the hypercube with a **hyperball** (L₂-norm ball):

$$\text{Accept iff } \lVert \mathbf{z}_i \rVert_2^2 \leq B^2$$

where $B^2 = \ell \cdot n \cdot \gamma_1^2$.

**Expected value** of $\lVert \mathbf{z}_i \rVert_2^2$ (dominated by uniform $\mathbf{y}_i$):

$$\mathbb{E}[\lVert \mathbf{z}_i \rVert_2^2] \approx \ell \cdot n \cdot \frac{\gamma_1^2}{3} = \frac{B^2}{3}$$

The bound is $3\times$ the expected value, improving acceptance versus naive per-party
L∞ checks. Actual protocol success still depends on aggregate signer behavior and retry policy.

### Numerical Values (ML-DSA-44)

| Quantity | Value |
|----------|-------|
| $\ell \cdot n$ | $4 \times 256 = 1{,}024$ |
| $\gamma_1^2$ | $131{,}072^2 = 1.72 \times 10^{10}$ |
| $B^2 = \ell \cdot n \cdot \gamma_1^2$ | $1.76 \times 10^{13}$ |
| $\mathbb{E}[\lVert \mathbf{z}_i \rVert_2^2]$ | $\approx 5.86 \times 10^{12}$ |
| Ratio $\mathbb{E}/B^2$ | $\approx 0.33$ |

---

## 9. Signature Aggregation and Encoding

### Aggregation (Coordinator)

**Step 1 — Validate partials:**

- **ADV-2 (Sybil):** Deduplicate by party_id — reject duplicate submissions
- **ADV-6 (Replay):** Verify each session binding: $\text{SHAKE-256}(\tilde{c} \| i) \stackrel{?}{=} \text{binding}_i$

**Step 2 — Sum:**

$$\mathbf{z} = \sum_{i \in \mathcal{Q}} \mathbf{z}_i = \sum_i (\mathbf{y}_i + c \cdot \mathbf{s}_{1,i}) = \sum_i \mathbf{y}_i + c \cdot \sum_i \mathbf{s}_{1,i}$$

Since $\sum_i \mathbf{s}_{1,i} = \mathbf{s}_1$ (RSS reconstruction via signing):

$$\mathbf{z} = \Big(\sum_i \mathbf{y}_i\Big) + c \cdot \mathbf{s}_1$$

This is **exactly** the standard ML-DSA response with the aggregate masking $\mathbf{y} = \sum_i \mathbf{y}_i$.

> [!WARNING]
> After summation, $\mathbf{z}$ must remain in **centered form** (via `reduce()` only, no `caddq()`). The `pack_z` function computes $\gamma_1 - z_{\text{coeff}}$ which requires centered coefficients in $[-(\gamma_1-1), \gamma_1]$. Mapping to $[0, q)$ via `caddq()` produces garbage.

**Step 3 — L∞ check on aggregate (constant-time):**

$$\lVert \mathbf{z} \rVert_\infty < \gamma_1 - \beta$$

If this fails, the entire protocol is retried.

### Hint Computation ($\mathbf{h}$)

FIPS 204 mandates the encoding of a hint vector $\mathbf{h}$. Because threshold signers do not distribute $\mathbf{s}_2$ or $\mathbf{t}_0$, the coordinator derives the approximation exactly as the verifier does:

$$\mathbf{w}_{\text{approx}} = \mathbf{A} \cdot \mathbf{z} - c \cdot \mathbf{t}_1 \cdot 2^d$$

> [!IMPORTANT]
> The subtraction is performed **in the NTT domain** with a single INTT at the end, matching the reference verification path exactly.
> The shift $\mathbf{t}_1 \cdot 2^d$ uses plain `<<=` without modular reduction, matching the reference `polyveck_shiftl`.

The hint bit is set where $\text{HighBits}(\mathbf{w}) \neq \text{HighBits}(\mathbf{w}_{\text{approx}})$. The `use_hint` function in the verifier corrects the discrepancy using the sign of $\text{LowBits}(\mathbf{w}_{\text{approx}})$.

### FIPS 204 Encoding

The final signature is encoded as:

$$\sigma = \tilde{c}\ \|\ \text{pack}_z(\mathbf{z})\ \|\ \text{pack}_h(\mathbf{h})$$

| Field | Size (bytes) | Encoding |
|-------|-------------|----------|
| $\tilde{c}$ | 32 | Raw challenge hash |
| $\mathbf{z}$ | $\ell \times 576 = 2{,}304$ | 18 bits/coefficient (centered, $\gamma_1$-bounded) |
| $\mathbf{h}$ | $\omega + k = 84$ | Sparse hint: ascending indices per polynomial |
| **Total** | **2,420** | Identical to standard ML-DSA-44 |

---

## 10. Verification

The threshold signature is verified using an **unmodified** FIPS 204 verifier:

$$\text{Verify}(pk, \text{msg}, \sigma) \to \{0, 1\}$$

**Verification equation:**

$$\mathbf{w}_1' = \text{UseHint}\big(\mathbf{h},\ \mathbf{A}\mathbf{z} - c \cdot \mathbf{t}_1 \cdot 2^d\big)$$

$$\tilde{c}' = \text{SHAKE-256}(\mu \| \text{pack}(\mathbf{w}_1'))$$

$$\text{Accept iff } \tilde{c} = \tilde{c}' \quad \text{(constant-time)}$$

**Why it works:** The aggregated $\mathbf{z}$ satisfies the same algebraic relation as standard ML-DSA:

$$\mathbf{A}\mathbf{z} = \mathbf{A}(\mathbf{y} + c \cdot \mathbf{s}_1) = \mathbf{A}\mathbf{y} + c \cdot (\mathbf{t} - \mathbf{s}_2) = \mathbf{w} + c \cdot \mathbf{t} - c \cdot \mathbf{s}_2$$

The HighBits extraction and hint recover $\mathbf{w}_1$, completing the Fiat-Shamir check.

---

## 11. Adversarial Security Hardening

The implementation defends against the following adversarial capabilities:

| ID | Attack | Defense | Mechanism |
|----|--------|---------|----------|
| ADV-1 | Coordinator grinds challenges by selectively including/excluding parties | Commitment binding | 4-round pre-commit/reveal with $H(\mathbf{w}_i \| i)$ |
| ADV-2 | Corrupt party submits multiple $\mathbf{z}_i$ to bias aggregate | Sybil protection | `seen_ids[MAX_PARTIES]` deduplication |
| ADV-3 | Coordinator observes individual $\mathbf{z}_i$ | **Inherent** | Mitigated by ADV-4 (no nonce reuse → no key extraction) |
| ADV-4 | RNG failure (VM snapshot, weak PRNG) → nonce reuse | Hedged nonces | $\text{seed} = H(\text{rng} \| \mathbf{s}_{1,i})$ |
| ADV-5 | Corrupt party submits garbage $\mathbf{w}_i$ | Binding verification | ADV-1 hash covers full $\mathbf{w}_i$ content |
| ADV-6 | Stale $\mathbf{z}_i$ replayed across sessions | Session binding | $H(\tilde{c} \| i)$ tag verified by coordinator |
| ADV-7 | Coordinator leaks acceptance statistics via iteration pattern | Fixed iteration | Always iterates ALL parties unconditionally |
| ADV-8 | Timing leak in $\mathbf{t}_1 \cdot 2^d \bmod q$ | Barrett reduction | Constant-time `reduce_i32()` |

All constant-time comparisons use `subtle::ConstantTimeEq`. All norm checks (`chknorm`) iterate all coefficients without early exit.

---

## 12. Sequence Diagrams

### 11.1 Key Generation and Distribution

```mermaid
sequenceDiagram
    participant D as Trusted Dealer
    participant P0 as Party 0
    participant P1 as Party 1
    participant P2 as Party 2

    Note over D: KeyGen
    D->>D: seed ← random(32)
    D->>D: (ρ, ρ', K) ← SHAKE-256(seed)
    D->>D: A ← ExpandA(ρ)
    D->>D: s₁ ←$ CBD_η, s₂ ←$ CBD_η
    D->>D: t = A·s₁ + s₂
    D->>D: pk = (ρ, t₁), sk = (ρ, K, tr, s₁, s₂, t₀)

    Note over D: RSS Distribution (N=3, T=2)
    D->>D: Enumerate subsets: {0,1}, {0,2}, {1,2}
    D->>D: s₁⁽⁰¹⁾ ←$ [-η, η]ⁿˡ
    D->>D: s₁⁽⁰²⁾ ←$ [-η, η]ⁿˡ
    D->>D: s₁⁽¹²⁾ = s₁ − s₁⁽⁰¹⁾ − s₁⁽⁰²⁾

    D->>P0: {s₁⁽⁰¹⁾, s₁⁽⁰²⁾}, {s₂⁽⁰¹⁾, s₂⁽⁰²⁾}, pk
    D->>P1: {s₁⁽⁰¹⁾, s₁⁽¹²⁾}, {s₂⁽⁰¹⁾, s₂⁽¹²⁾}, pk
    D->>P2: {s₁⁽⁰²⁾, s₁⁽¹²⁾}, {s₂⁽⁰²⁾, s₂⁽¹²⁾}, pk

    P0->>P0: s₁₀ = s₁⁽⁰¹⁾ + s₁⁽⁰²⁾
    P1->>P1: s₁₁ = s₁⁽⁰¹⁾ + s₁⁽¹²⁾
    P2->>P2: s₁₂ = s₁⁽⁰²⁾ + s₁⁽¹²⁾

    D->>D: Zeroize s₁, s₂, seed
```

### 12.2 Threshold Signing Protocol (4 Rounds — Hardened)

```mermaid
sequenceDiagram
    participant P0 as Party 0
    participant P1 as Party 1
    participant P2 as Party 2
    participant C as Coordinator

    Note over P0,C: ═══ ROUND 0: PRE-COMMIT (ADV-1) ═══

    P0->>P0: seed₀ = H(rng ‖ s₁₀)
    P0->>P0: y₀ ← ExpandMask(seed₀)
    P0->>P0: w₀ = INTT(A_hat · NTT(y₀))
    P0->>P0: h₀ = H(w₀ ‖ 0)
    P1->>P1: seed₁ = H(rng ‖ s₁₁)
    P1->>P1: y₁, w₁, h₁
    P2->>P2: seed₂ = H(rng ‖ s₁₂)
    P2->>P2: y₂, w₂, h₂

    P0->>C: h₀
    P1->>C: h₁
    P2->>C: h₂

    Note over P0,C: ═══ ROUND 1: REVEAL ═══

    P0->>C: w₀
    P1->>C: w₁
    P2->>C: w₂

    C->>C: Verify H(w₀‖0) = h₀ ✓
    C->>C: Verify H(w₁‖1) = h₁ ✓
    C->>C: Verify H(w₂‖2) = h₂ ✓

    Note over P0,C: ═══ ROUND 2: CHALLENGE ═══

    C->>C: w = w₀ + w₁ + w₂
    C->>C: μ = H(tr ‖ msg)
    C->>C: c̃ = H(μ ‖ pack_w1(HighBits(w)))

    C->>P0: c̃
    C->>P1: c̃
    C->>P2: c̃

    Note over P0,C: ═══ ROUND 3: SIGN (ADV-6) ═══

    P0->>P0: z₀ = c · s₁₀ + y₀, check ‖z₀‖₂² ≤ B²
    P0->>P0: binding₀ = H(c̃ ‖ 0)
    P1->>P1: z₁ = c · s₁₁ + y₁, check ‖z₁‖₂² ≤ B²
    P1->>P1: binding₁ = H(c̃ ‖ 1)
    P2->>P2: z₂ = c · s₁₂ + y₂, check ‖z₂‖₂² ≤ B²
    P2->>P2: binding₂ = H(c̃ ‖ 2)

    P0->>C: (z₀, binding₀) ✓
    P1->>C: (z₁, binding₁) ✓
    P2->>C: (z₂, binding₂) ✓

    Note over P0,C: ═══ AGGREGATION ═══

    C->>C: Dedup party_ids (ADV-2)
    C->>C: Verify all bindings (ADV-6)
    C->>C: z = z₀ + z₁ + z₂ (centered form)
    C->>C: Check ‖z‖∞ < γ₁ − β (constant-time)
    C->>C: w_approx = Az - c·t₁·2ᵈ (single-INTT)
    C->>C: h = hint where HighBits(w) ≠ HighBits(w_approx)
    C->>C: σ = (c̃, pack_z(z), pack_h(h))
    C-->>C: Output σ
```

### 11.3 Signing with Rejection and Retry

```mermaid
sequenceDiagram
    participant P0 as Party 0
    participant P1 as Party 1
    participant P2 as Party 2
    participant C as Coordinator

    Note over P0,C: Attempt 1

    P0->>C: w₀
    P1->>C: w₁
    P2->>C: w₂
    C->>P0: c̃
    C->>P1: c̃
    C->>P2: c̃

    P0->>P0: z₀: ‖z₀‖₂² > B² → REJECT
    P1->>C: z₁ ✓
    P2->>C: z₂ ✓

    C->>C: Only 2/3 passed (≥ T=2) → aggregate
    C->>C: z = z₁ + z₂
    C->>C: ‖z‖∞ ≥ γ₁−β → FAIL

    Note over P0,C: Attempt 2 (full retry)

    P0->>C: w₀'
    P1->>C: w₁'
    P2->>C: w₂'
    C->>P0: c̃'
    C->>P1: c̃'
    C->>P2: c̃'

    P0->>C: z₀' ✓
    P1->>C: z₁' ✓
    P2->>C: z₂' ✓

    C->>C: z = z₀' + z₁' + z₂'
    C->>C: ‖z‖∞ < γ₁−β → ✓
    C-->>C: σ = (c̃', z, h)
```

---

## 12. Activity Diagrams

### 12.1 Overall Protocol Flow

```mermaid
flowchart TD
    START([Start]) --> KG[Key Generation]
    KG --> RSS[RSS Distribution to N parties]
    RSS --> INIT[Parties compute effective shares]

    INIT --> R1["ROUND 1: All parties sample yᵢ, compute wᵢ"]
    R1 --> AGG_W["Coordinator: w = Σ wᵢ"]
    AGG_W --> CHAL["Coordinator: c̃ = H(tr ‖ w₁ ‖ msg)"]
    CHAL --> BCAST["Broadcast c̃ to all parties"]

    BCAST --> R3["ROUND 3: Each party computes zᵢ = c·s₁ᵢ + yᵢ"]
    R3 --> REJ{"Each party:\n‖zᵢ‖₂² ≤ B²?"}
    REJ -->|Accept| SEND["Send zᵢ to coordinator"]
    REJ -->|Reject| ABORT["Party aborts locally"]

    SEND --> COUNT{"Coordinator:\n≥ T responses?"}
    ABORT --> COUNT
    COUNT -->|No| RETRY{"Retries left?"}
    RETRY -->|Yes| R1
    RETRY -->|No| FAIL([Protocol Failed])

    COUNT -->|Yes| AGG_Z["Coordinator: z = Σ zᵢ"]
    AGG_Z --> NORM{"‖z‖∞ < γ₁ − β?"}
    NORM -->|No| RETRY
    NORM -->|Yes| ENCODE["Compute h = MakeHint(...)<br>Encode σ = (c̃, z, h)"]
    ENCODE --> VERIFY["Verify via FIPS 204"]
    VERIFY --> DONE([Valid Signature σ])

    style START fill:#339af0,color:#fff
    style DONE fill:#51cf66,color:#fff
    style FAIL fill:#ff6b6b,color:#fff
    style REJ fill:#ffa94d,color:#000
    style NORM fill:#ffa94d,color:#000
```

### 12.2 Per-Party Sign Operation (Round 3)

```mermaid
flowchart TD
    START([Receive c̃]) --> DECODE["c = SampleInBall(c̃)"]
    DECODE --> NTT_C["ĉ = NTT(c)"]
    NTT_C --> NTT_S["ŝ₁ᵢ = NTT(s₁ᵢ)"]

    NTT_S --> LOOP["For j = 0 to ℓ−1:"]
    LOOP --> PW["cs₁ⱼ = INTT(ĉ ⊙ ŝ₁ᵢⱼ)"]
    PW --> RED1["reduce(cs₁ⱼ), caddq(cs₁ⱼ)"]
    RED1 --> ADD["zⱼ = cs₁ⱼ + yⱼ"]
    ADD --> RED2["reduce(zⱼ)"]
    RED2 --> NEXT{j < ℓ−1?}
    NEXT -->|Yes| LOOP
    NEXT -->|No| NORM

    NORM["Compute ‖z‖₂² = Σⱼ Σₘ center(zⱼₘ)²"]
    NORM --> CHECK{"‖z‖₂² ≤ B²?"}
    CHECK -->|Yes| ACCEPT(["Return zᵢ ✓"])
    CHECK -->|No| REJECT(["Return LocalRejectionAbort ✗"])

    style START fill:#339af0,color:#fff
    style ACCEPT fill:#51cf66,color:#fff
    style REJECT fill:#ff6b6b,color:#fff
    style CHECK fill:#ffa94d,color:#000
```

### 12.3 RSS Key Distribution

```mermaid
flowchart TD
    START([Input: s₁, s₂, N, T]) --> VALIDATE{"N ≤ 6?\nT ≤ N?\nT ≥ 2?"}
    VALIDATE -->|No| ERROR([InvalidParameters])
    VALIDATE -->|Yes| ENUM["Enumerate all C(N,N−T+1) subsets"]

    ENUM --> LOOP["For each subset S (except last):"]
    LOOP --> SAMPLE["Sample s₁⁽ˢ⁾ ←$ [-η, η]ⁿˡ"]
    SAMPLE --> SUM["sum += s₁⁽ˢ⁾"]
    SUM --> NEXT{More subsets?}
    NEXT -->|Yes| LOOP
    NEXT -->|No| LAST

    LAST["s₁⁽ˡᵃˢᵗ⁾ = s₁ − sum"] --> DIST["Distribute to parties"]
    DIST --> PARTY["Party i receives all s⁽ˢ⁾ where i ∈ S"]
    PARTY --> EFF["Party has local holding s₁,holdᵢ = Σ_{S∋i} s₁⁽ˢ⁾"]
    EFF --> DONE([Party Key Shares Ready])

    style START fill:#339af0,color:#fff
    style DONE fill:#51cf66,color:#fff
    style ERROR fill:#ff6b6b,color:#fff
```

---

## 14. Security Properties

### Correctness

The threshold signature is **identical** in format to a standard ML-DSA signature. The verification equation holds because:

$$\mathbf{z} = \sum_{i \in \mathcal{A}} \mathbf{z}_i
= \sum_{i \in \mathcal{A}} (\mathbf{y}_i + c \cdot \mathbf{s}^{\text{sign}}_{1,i})
= \underbrace{\sum_{i \in \mathcal{A}} \mathbf{y}_i}_{\text{aggregate mask}} + c \cdot \underbrace{\sum_{i \in \mathcal{A}} \mathbf{s}^{\text{sign}}_{1,i}}_{\mathbf{s}_1}$$

### Unforgeability

The threshold scheme inherits the _Module-LWE/Module-SIS_ hardness assumptions from ML-DSA. An adversary controlling fewer than $T$ parties cannot reconstruct $\mathbf{s}_1$ because:

- They are missing at least one RSS share piece (for a subset they don't fully cover)
- Each missing piece is uniformly random in $[-\eta, \eta]^{n \cdot \ell}$, information-theoretically hiding the secret

### Privacy

- Share pieces are $\eta$-bounded random values — statistically independent of the secret
- The last share (computed as $\mathbf{s} - \sum \text{random}$) has larger but still bounded magnitude
- Sensitive material ($\mathbf{s}_{1,i}$, $\mathbf{y}_i$, seed) is zeroized on drop via `zeroize`
- Challenge and session binding comparisons use constant-time equality (`subtle::ConstantTimeEq`)
- All norm checks (`chknorm`) iterate all coefficients/polynomials without early exit

### Adversarial Resilience

- **Malicious coordinator** cannot grind challenges (commitment binding, ADV-1)
- **Malicious party** cannot bias the aggregate signature (deduplication, ADV-2)
- **RNG failure** cannot cause nonce reuse across parties (hedged derivation, ADV-4)
- **Replay attacks** are detected via session binding (ADV-6)
- **Timing attacks** are mitigated by constant-time norms, comparisons, and fixed iteration (ALG-4, ADV-7)

### Liveness

Liveness depends on concrete parameterization, active-set choice, and rejection outcomes.
The implementation maintains safety by retrying with qualifying signer sets and returning
an error if no verifiable aggregate is produced within the configured retry budget.

---

> **References:**
>
> - NIST FIPS 204 — _Module-Lattice-Based Digital Signature Standard_ (2024)
> - ePrint 2026/013 — _Mithril: Threshold ML-DSA via Replicated Secret Sharing_
> - [`dilithium-rs`](https://crates.io/crates/dilithium-rs) — Pure-Rust FIPS 204 implementation (used for NTT, Montgomery reduction, and canonical verification)
