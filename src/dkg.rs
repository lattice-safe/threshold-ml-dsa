//! Distributed Key Generation (DKG) for threshold ML-DSA-44.
//!
//! Produces the same replicated-secret-sharing key structure as
//! [`crate::rss::keygen_from_seed`] without any party learning the full
//! secret: every (N−T+1)-subset secret `s_I` is dealt by one of its own
//! members, `d(I) = min(I)`, and sent only to the other members of `I`.
//!
//! This construction is this crate's own design. It has not been checked
//! against ePrint 2026/013 §4.3, and the Go reference has no DKG.
//!
//! ## Rounds
//!
//! 0. Every party broadcasts a fresh 32-byte nonce. The session ID is
//!    `sid = H(tag ‖ t ‖ n ‖ context ‖ nonce_0 ‖ … ‖ nonce_{N−1})` and the
//!    matrix seed is `ρ = H(tag ‖ sid)`.
//! 1. For each subset it deals, a party samples `seed_I`, derives the short
//!    secret `s_I` from it, computes `t_I = A·s₁_I + s₂_I` and broadcasts the
//!    commitment `C_I = H(tag ‖ sid ‖ dealer ‖ I ‖ t_I)`. It sends `seed_I`
//!    privately to every other member of `I`.
//! 2. Members re-derive `s_I` from the received seed (so `s_I` is short by
//!    construction even for a malicious dealer) and check it against `C_I`.
//!    Then every party reveals its `t_I`.
//! 3. Everyone checks each `t_I` against `C_I`, computes `t = Σ t_I`, the
//!    public key and a confirmation digest over `(sid, pk, all C_I)`.
//! 4. Parties exchange digests; the key is output only if all agree.
//!
//! ## Security notes
//!
//! - A coalition of at most T−1 parties misses at least one subset `M`
//!   (`M ∩ C = ∅`); its dealer is honest and `seed_M` never leaves `M`. The
//!   coalition sees `t_M`, a fresh MLWE sample under the same `A` as an
//!   ML-DSA public key.
//! - Commit-then-reveal stops a party that holds a whole subset from choosing
//!   its `t_I` after seeing the others (rogue-key attack).
//! - The confirm round detects equivocation in the broadcast.
//! - A party can still abort after seeing the reveals (attributable over an
//!   authenticated broadcast); restart without it and with a fresh session.
//! - A seed complaint is not attributable (the dealer or the member may be
//!   lying); abort the session.
//!
//! ## Transport requirements
//!
//! - Broadcast messages need an authenticated channel; they are public.
//! - [`DkgSeedMsg`] carries a secret and MUST be sent over an authenticated
//!   and encrypted point-to-point channel to `to`.

extern crate alloc;
use alloc::{collections::BTreeMap, vec::Vec};

use crate::error::Error;
use crate::params::{get_threshold_params, PK_BYTES, Q};
use crate::poly::PolyVecK;
use crate::rss::{
    compute_t, pk_from_t, share_from_seed, subset_masks, tr_from_pk, Share, ThresholdPrivateKey,
};
use crate::sign::{pack_w_single, pack_w_single_size, unpack_w_single};
use rand_core::{CryptoRng, RngCore};
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use subtle::ConstantTimeEq;
use zeroize::{Zeroize, Zeroizing};

/// Configuration of one party in a DKG session.
#[derive(Clone, Copy, Debug)]
pub struct DkgConfig {
    /// Threshold T
    pub t: u8,
    /// Number of parties N
    pub n: u8,
    /// This party's ID ∈ [0, N)
    pub me: u8,
    /// Application context bound into the session ID, e.g. a hash of the
    /// participants' long-term identity keys.
    pub context: [u8; 32],
}

/// Round 0 broadcast: a fresh session nonce.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgRound0Msg {
    /// Sender
    pub from: u8,
    /// Fresh random nonce
    pub nonce: [u8; 32],
}

/// Round 1 broadcast: commitments to `t_I` for each subset the sender deals.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgCommitMsg {
    /// Sender
    pub from: u8,
    /// Session ID
    pub sid: [u8; 32],
    /// `(subset mask, C_I)` in ascending mask order
    pub commits: Vec<(u8, [u8; 32])>,
}

/// Round 2 broadcast: the committed `t_I` (23-bit packed, coefficients in `[0, q)`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgRevealMsg {
    /// Sender
    pub from: u8,
    /// Session ID
    pub sid: [u8; 32],
    /// `(subset mask, packed t_I)` in ascending mask order
    pub t_parts: Vec<(u8, Vec<u8>)>,
}

/// Round 3 broadcast: digest of the sender's view of the transcript.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgConfirmMsg {
    /// Sender
    pub from: u8,
    /// `H(tag ‖ sid ‖ pk ‖ (mask ‖ C_I)*)`
    pub digest: [u8; 32],
}

/// Round 1 private message: the seed of subset `mask`, from its dealer to
/// one other member.
///
/// **Secret.** Send it only over an authenticated, encrypted channel to `to`.
/// Not `Clone`/`Debug`; the seed is wiped on drop.
pub struct DkgSeedMsg {
    from: u8,
    to: u8,
    sid: [u8; 32],
    mask: u8,
    seed: Zeroizing<[u8; 64]>,
}

const SEED_MSG_VERSION: u8 = 1;
const SEED_MSG_BYTES: usize = 1 + 1 + 1 + 32 + 1 + 64;

impl DkgSeedMsg {
    /// Sender (the subset's dealer).
    #[must_use]
    pub fn from(&self) -> u8 {
        self.from
    }

    /// Recipient.
    #[must_use]
    pub fn to(&self) -> u8 {
        self.to
    }

    /// Subset mask.
    #[must_use]
    pub fn mask(&self) -> u8 {
        self.mask
    }

    /// Encode for transport: `1 ‖ from ‖ to ‖ sid ‖ mask ‖ seed`.
    /// The output is secret and wiped on drop.
    #[must_use]
    pub fn to_bytes(&self) -> Zeroizing<Vec<u8>> {
        let mut out = Zeroizing::new(Vec::with_capacity(SEED_MSG_BYTES));
        out.extend_from_slice(&[SEED_MSG_VERSION, self.from, self.to]);
        out.extend_from_slice(&self.sid);
        out.push(self.mask);
        out.extend_from_slice(&self.seed[..]);
        out
    }

    /// Decode a message produced by [`DkgSeedMsg::to_bytes`].
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error> {
        if bytes.len() != SEED_MSG_BYTES || bytes[0] != SEED_MSG_VERSION {
            return Err(Error::DkgMalformedMessage);
        }
        let mut seed = Zeroizing::new([0u8; 64]);
        seed.copy_from_slice(&bytes[36..100]);
        Ok(Self {
            from: bytes[1],
            to: bytes[2],
            sid: bytes[3..35].try_into().map_err(|_| Error::DkgMalformedMessage)?,
            mask: bytes[35],
            seed,
        })
    }
}

/// A party after Round 0, waiting for every party's nonce.
pub struct DkgParty {
    cfg: DkgConfig,
    nonce: [u8; 32],
}

/// A party after Round 1, waiting for commitments and its seeds.
pub struct DkgAwaitingCommits {
    cfg: DkgConfig,
    sid: [u8; 32],
    rho: [u8; 32],
    key: Zeroizing<[u8; 32]>,
    /// Subsets this party deals: mask → (share, t_I)
    dealt: BTreeMap<u8, (Share, PolyVecK)>,
    own_commits: Vec<(u8, [u8; 32])>,
}

/// A party after Round 2, waiting for reveals.
pub struct DkgAwaitingReveals {
    cfg: DkgConfig,
    sid: [u8; 32],
    rho: [u8; 32],
    key: Zeroizing<[u8; 32]>,
    shares: BTreeMap<u8, Share>,
    /// All commitments, mask → C_I
    commits: BTreeMap<u8, [u8; 32]>,
}

/// A party after Round 3, waiting for confirmations.
pub struct DkgAwaitingConfirm {
    cfg: DkgConfig,
    rho: [u8; 32],
    key: Zeroizing<[u8; 32]>,
    shares: BTreeMap<u8, Share>,
    pk: [u8; PK_BYTES],
    digest: [u8; 32],
}

/// Dealer of a subset: its lowest member.
fn dealer(mask: u8) -> u8 {
    mask.trailing_zeros() as u8
}

fn hash32(parts: &[&[u8]]) -> [u8; 32] {
    let mut h = Shake256::default();
    for p in parts {
        h.update(p);
    }
    let mut out = [0u8; 32];
    h.finalize_xof().read(&mut out);
    out
}

fn commitment(sid: &[u8; 32], mask: u8, t_packed: &[u8]) -> [u8; 32] {
    hash32(&[b"th-ml-dsa-dkg-commit-v2", sid, &[dealer(mask), mask], t_packed])
}

/// Unpack a `t_I`, rejecting wrong lengths and coefficients ≥ q.
fn unpack_t_canonical(bytes: &[u8]) -> Result<PolyVecK, Error> {
    if bytes.len() != pack_w_single_size() {
        return Err(Error::DkgMalformedMessage);
    }
    let t = unpack_w_single(bytes);
    if t.polys.iter().any(|p| p.coeffs.iter().any(|&c| c >= Q)) {
        return Err(Error::DkgMalformedMessage);
    }
    Ok(t)
}

/// Check that `froms` is exactly `0..n`, in any order; return the order.
fn check_senders(froms: impl Iterator<Item = u8>, n: u8) -> Result<(), Error> {
    let mut seen = 0u16;
    let mut count = 0;
    for f in froms {
        if f >= n || seen & (1 << f) != 0 {
            return Err(Error::DkgMalformedMessage);
        }
        seen |= 1 << f;
        count += 1;
    }
    if count != n as usize {
        return Err(Error::DkgMalformedMessage);
    }
    Ok(())
}

impl DkgParty {
    /// Round 0: validate the configuration and draw this party's nonce.
    pub fn start<R: RngCore + CryptoRng>(cfg: DkgConfig, rng: &mut R) -> Result<(Self, DkgRound0Msg), Error> {
        get_threshold_params(cfg.t, cfg.n).ok_or(Error::InvalidParameters)?;
        if cfg.me >= cfg.n {
            return Err(Error::InvalidParameters);
        }
        let mut nonce = [0u8; 32];
        rng.fill_bytes(&mut nonce);
        Ok((Self { cfg, nonce }, DkgRound0Msg { from: cfg.me, nonce }))
    }

    /// Round 1: derive the session, deal this party's subsets.
    ///
    /// Returns the commitment broadcast and the private seed messages, one per
    /// (subset, other member).
    pub fn round1<R: RngCore + CryptoRng>(
        self,
        r0: &[DkgRound0Msg],
        rng: &mut R,
    ) -> Result<(DkgAwaitingCommits, DkgCommitMsg, Vec<DkgSeedMsg>), Error> {
        let cfg = self.cfg;
        check_senders(r0.iter().map(|m| m.from), cfg.n)?;
        let mut nonces = alloc::vec![[0u8; 32]; cfg.n as usize];
        for m in r0 {
            nonces[m.from as usize] = m.nonce;
        }
        if nonces[cfg.me as usize] != self.nonce {
            return Err(Error::DkgInconsistentTranscript);
        }

        let mut h = Shake256::default();
        h.update(b"th-ml-dsa-dkg-sid-v2");
        h.update(&[cfg.t, cfg.n]);
        h.update(&cfg.context);
        for nonce in &nonces {
            h.update(nonce);
        }
        let mut sid = [0u8; 32];
        h.finalize_xof().read(&mut sid);
        let rho = hash32(&[b"th-ml-dsa-dkg-rho-v2", &sid]);

        let mut key = Zeroizing::new([0u8; 32]);
        rng.fill_bytes(&mut key[..]);

        let mut dealt = BTreeMap::new();
        let mut own_commits = Vec::new();
        let mut seed_msgs = Vec::new();
        for mask in subset_masks(cfg.n, cfg.t) {
            if dealer(mask) != cfg.me {
                continue;
            }
            let mut seed = Zeroizing::new([0u8; 64]);
            rng.fill_bytes(&mut seed[..]);
            let share = share_from_seed(&seed);
            let t_i = compute_t(&rho, &share.s1h, &share.s2);
            own_commits.push((mask, commitment(&sid, mask, &pack_w_single(&t_i))));
            for to in 0..cfg.n {
                if to != cfg.me && mask & (1 << to) != 0 {
                    seed_msgs.push(DkgSeedMsg { from: cfg.me, to, sid, mask, seed: seed.clone() });
                }
            }
            dealt.insert(mask, (share, t_i));
        }

        let msg = DkgCommitMsg { from: cfg.me, sid, commits: own_commits.clone() };
        Ok((DkgAwaitingCommits { cfg, sid, rho, key, dealt, own_commits }, msg, seed_msgs))
    }
}

impl DkgAwaitingCommits {
    /// Round 2: check all commitments and the seeds received from other
    /// dealers, then reveal this party's `t_I`.
    pub fn round2(
        self,
        commits: &[DkgCommitMsg],
        seeds_for_me: Vec<DkgSeedMsg>,
    ) -> Result<(DkgAwaitingReveals, DkgRevealMsg), Error> {
        let cfg = self.cfg;
        check_senders(commits.iter().map(|m| m.from), cfg.n)?;
        let masks = subset_masks(cfg.n, cfg.t);

        let mut all = BTreeMap::new();
        for m in commits {
            if !bool::from(m.sid.ct_eq(&self.sid)) {
                return Err(Error::DkgMalformedMessage);
            }
            let expected: Vec<u8> = masks.iter().copied().filter(|&k| dealer(k) == m.from).collect();
            let got: Vec<u8> = m.commits.iter().map(|c| c.0).collect();
            if got != expected {
                return Err(Error::DkgMalformedMessage);
            }
            if m.from == cfg.me && m.commits != self.own_commits {
                return Err(Error::DkgInconsistentTranscript);
            }
            for &(mask, c) in &m.commits {
                all.insert(mask, c);
            }
        }

        // Exactly one seed for each subset containing me that I do not deal.
        let mut shares = BTreeMap::new();
        let wanted: Vec<u8> = masks
            .iter()
            .copied()
            .filter(|&k| k & (1 << cfg.me) != 0 && dealer(k) != cfg.me)
            .collect();
        if seeds_for_me.len() != wanted.len() {
            return Err(Error::DkgMalformedMessage);
        }
        for msg in seeds_for_me {
            if msg.to != cfg.me
                || !wanted.contains(&msg.mask)
                || msg.from != dealer(msg.mask)
                || !bool::from(msg.sid.ct_eq(&self.sid))
                || shares.contains_key(&msg.mask)
            {
                return Err(Error::DkgMalformedMessage);
            }
            let share = share_from_seed(&msg.seed);
            let t_i = compute_t(&self.rho, &share.s1h, &share.s2);
            let c = commitment(&self.sid, msg.mask, &pack_w_single(&t_i));
            if !bool::from(c.ct_eq(&all[&msg.mask])) {
                return Err(Error::DkgSeedMismatch { dealer: msg.from, mask: msg.mask });
            }
            shares.insert(msg.mask, share);
        }

        let mut t_parts = Vec::new();
        for (mask, (share, t_i)) in self.dealt {
            t_parts.push((mask, pack_w_single(&t_i)));
            if mask & (1 << cfg.me) != 0 {
                shares.insert(mask, share);
            }
        }

        let reveal = DkgRevealMsg { from: cfg.me, sid: self.sid, t_parts };
        Ok((
            DkgAwaitingReveals { cfg, sid: self.sid, rho: self.rho, key: self.key, shares, commits: all },
            reveal,
        ))
    }
}

impl DkgAwaitingReveals {
    /// Round 3: check every reveal against its commitment, compute the
    /// public key and this party's confirmation digest.
    pub fn round3(self, reveals: &[DkgRevealMsg]) -> Result<(DkgAwaitingConfirm, DkgConfirmMsg), Error> {
        let cfg = self.cfg;
        check_senders(reveals.iter().map(|m| m.from), cfg.n)?;

        let mut t = PolyVecK::zero();
        for m in reveals {
            if !bool::from(m.sid.ct_eq(&self.sid)) {
                return Err(Error::DkgMalformedMessage);
            }
            let expected: Vec<u8> = self.commits.keys().copied().filter(|&k| dealer(k) == m.from).collect();
            let got: Vec<u8> = m.t_parts.iter().map(|p| p.0).collect();
            if got != expected {
                return Err(Error::DkgMalformedMessage);
            }
            for (mask, packed) in &m.t_parts {
                let t_i = unpack_t_canonical(packed)?;
                let c = commitment(&self.sid, *mask, packed);
                if !bool::from(c.ct_eq(&self.commits[mask])) {
                    return Err(Error::DkgCommitmentMismatch { party: m.from });
                }
                t.add_assign(&t_i);
                t.reduce();
            }
        }
        // reduce() leaves (-q, q); bring t into [0, q) for Power2Round.
        for p in &mut t.polys {
            for c in &mut p.coeffs {
                *c = c.rem_euclid(Q);
            }
        }
        let pk = pk_from_t(&self.rho, &t);

        let mut h = Shake256::default();
        h.update(b"th-ml-dsa-dkg-confirm-v2");
        h.update(&self.sid);
        h.update(&pk);
        for (mask, c) in &self.commits {
            h.update(&[*mask]);
            h.update(c);
        }
        let mut digest = [0u8; 32];
        h.finalize_xof().read(&mut digest);

        Ok((
            DkgAwaitingConfirm { cfg, rho: self.rho, key: self.key, shares: self.shares, pk, digest },
            DkgConfirmMsg { from: cfg.me, digest },
        ))
    }
}

impl DkgAwaitingConfirm {
    /// Final step: output the key only if every party saw the same transcript.
    pub fn finalize(mut self, confirms: &[DkgConfirmMsg]) -> Result<([u8; PK_BYTES], ThresholdPrivateKey), Error> {
        check_senders(confirms.iter().map(|m| m.from), self.cfg.n)?;
        for m in confirms {
            if !bool::from(m.digest.ct_eq(&self.digest)) {
                return Err(Error::DkgInconsistentTranscript);
            }
        }
        let shares = core::mem::take(&mut self.shares);
        let sk = ThresholdPrivateKey {
            t: self.cfg.t,
            n: self.cfg.n,
            id: self.cfg.me,
            rho: self.rho,
            key: *self.key,
            tr: tr_from_pk(&self.pk),
            shares,
        };
        self.key.zeroize();
        Ok((self.pk, sk))
    }
}

/// All messages of a simulated DKG run (see [`simulate`]).
pub struct DkgTranscript {
    /// Round 0 broadcasts
    pub round0: Vec<DkgRound0Msg>,
    /// Round 1 broadcasts
    pub commits: Vec<DkgCommitMsg>,
    /// Round 1 private messages
    pub seeds: Vec<DkgSeedMsg>,
    /// Round 2 broadcasts
    pub reveals: Vec<DkgRevealMsg>,
    /// Round 3 broadcasts
    pub confirms: Vec<DkgConfirmMsg>,
}

/// Run all N parties of a DKG in one process (tests and benchmarks only).
///
/// Holding every party's output in one place defeats the purpose of a DKG;
/// real deployments run one [`DkgParty`] per participant.
pub fn simulate<R: RngCore + CryptoRng>(
    t: u8,
    n: u8,
    context: [u8; 32],
    rng: &mut R,
) -> Result<([u8; PK_BYTES], Vec<ThresholdPrivateKey>, DkgTranscript), Error> {
    let mut parties = Vec::new();
    let mut round0 = Vec::new();
    for me in 0..n {
        let (p, m) = DkgParty::start(DkgConfig { t, n, me, context }, rng)?;
        parties.push(p);
        round0.push(m);
    }
    let mut st1 = Vec::new();
    let mut commits = Vec::new();
    let mut seeds = Vec::new();
    for p in parties {
        let (s, c, sm) = p.round1(&round0, rng)?;
        st1.push(s);
        commits.push(c);
        seeds.extend(sm);
    }
    let mut st2 = Vec::new();
    let mut reveals = Vec::new();
    for (me, s) in st1.into_iter().enumerate() {
        let mine: Vec<DkgSeedMsg> = seeds
            .iter()
            .filter(|m| m.to as usize == me)
            .map(|m| DkgSeedMsg { from: m.from, to: m.to, sid: m.sid, mask: m.mask, seed: m.seed.clone() })
            .collect();
        let (s2, r) = s.round2(&commits, mine)?;
        st2.push(s2);
        reveals.push(r);
    }
    let mut st3 = Vec::new();
    let mut confirms = Vec::new();
    for s in st2 {
        let (s3, c) = s.round3(&reveals)?;
        st3.push(s3);
        confirms.push(c);
    }
    let mut pk = [0u8; PK_BYTES];
    let mut sks = Vec::new();
    for s in st3 {
        let (p, sk) = s.finalize(&confirms)?;
        pk = p;
        sks.push(sk);
    }
    Ok((pk, sks, DkgTranscript { round0, commits, seeds, reveals, confirms }))
}
