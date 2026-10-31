//! Distributed Key Generation (DKG) for threshold ML-DSA-44.
//!
//! Produces the same replicated-secret-sharing key structure as
//! [`crate::rss::keygen_from_seed`] without any party learning the full
//! secret. It follows the outline in the Mithril slides (NIST MPTS 2026,
//! slides 18–19): exchange a key `K_I` within each (N−T+1)-subset `I`,
//! jointly sample a coin, derive `s_I = H(coin, K_I)`, commit-and-reveal
//! `t_I = A·s₁_I + s₂_I`, and set `t = Σ t_I`. The USENIX version of
//! ePrint 2026/013 defers the DKG to the full version, which has not been
//! checked; message formats and checks here are this crate's own.
//!
//! ## Rounds
//!
//! 0. Every party samples a coin share `c_i` and, for each subset it deals
//!    (`d(I) = min(I)`), a key `K_I`. It broadcasts `H(c_i)` and every
//!    `H(K_I)`. The session ID is the hash of all round-0 messages.
//! 1. Dealers send `K_I` privately to the other members of `I`. Every party
//!    reveals `c_i`.
//! 2. Members check `K_I` against its commitment and every `c_i` against
//!    its commitment. Then `coin = H(sid ‖ c_0 ‖ … ‖ c_{N−1})`,
//!    `ρ = H(coin)` and `s_I = H(sid ‖ coin ‖ I ‖ K_I)`, which is short by
//!    construction. Dealers broadcast `C_I = H(sid ‖ I ‖ t_I)`.
//! 3. Members check `C_I` against their own `t_I`. Dealers reveal `t_I`.
//! 4. Everyone checks each `t_I` against `C_I`, computes `t = Σ t_I`, the
//!    public key and a digest of the transcript.
//! 5. Parties exchange digests; the key is output only if all agree.
//!
//! ## Security notes
//!
//! - A coalition of at most T−1 parties misses at least one subset `M`
//!   (`M ∩ C = ∅`); its dealer is honest and `K_M` never leaves `M`. The
//!   coalition sees `t_M`, a fresh MLWE sample under the same `A` as an
//!   ML-DSA public key (the paper also publishes every `vk_I`).
//! - `K_I` and `c_i` are committed before the coin is revealed, so a
//!   dealer cannot grind `s_I` for a subset that has an honest member, and
//!   no party can bias `ρ` other than by aborting.
//! - Commit-then-reveal of `t_I` stops a coalition that holds a whole
//!   subset from choosing its `t_I` after seeing the others (rogue key).
//! - The confirm round detects equivocation in the broadcast.
//! - A party can still abort after seeing the coin or the reveals
//!   (attributable over an authenticated broadcast); restart without it.
//! - A key complaint is not attributable (the dealer or the member may be
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

/// Round 0 broadcast: commitments to the sender's coin share and subset keys.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgRound0Msg {
    /// Sender
    pub from: u8,
    /// `H(tag ‖ t ‖ n ‖ context ‖ from ‖ c_i)`
    pub coin_commit: [u8; 32],
    /// `(subset mask, H(tag ‖ t ‖ n ‖ context ‖ dealer ‖ mask ‖ K_I))` for
    /// each subset the sender deals, in ascending mask order
    pub key_commits: Vec<(u8, [u8; 32])>,
}

/// Round 1 broadcast: the sender's coin share.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgCoinMsg {
    /// Sender
    pub from: u8,
    /// Session ID
    pub sid: [u8; 32],
    /// Coin share `c_i`
    pub coin: [u8; 32],
}

/// Round 2 broadcast: commitments to `t_I` for each subset the sender deals.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgCommitMsg {
    /// Sender
    pub from: u8,
    /// Session ID
    pub sid: [u8; 32],
    /// `(subset mask, C_I)` in ascending mask order
    pub commits: Vec<(u8, [u8; 32])>,
}

/// Round 3 broadcast: the committed `t_I` (23-bit packed, coefficients in `[0, q)`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgRevealMsg {
    /// Sender
    pub from: u8,
    /// Session ID
    pub sid: [u8; 32],
    /// `(subset mask, packed t_I)` in ascending mask order
    pub t_parts: Vec<(u8, Vec<u8>)>,
}

/// Round 4 broadcast: digest of the sender's view of the transcript.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DkgConfirmMsg {
    /// Sender
    pub from: u8,
    /// `H(tag ‖ sid ‖ coin ‖ pk ‖ (mask ‖ C_I)*)`
    pub digest: [u8; 32],
}

/// Round 1 private message: the key `K_I` of subset `mask`, from its dealer
/// to one other member.
///
/// **Secret.** Send it only over an authenticated, encrypted channel to `to`.
/// Not `Clone`/`Debug`; the key is wiped on drop.
pub struct DkgSeedMsg {
    from: u8,
    to: u8,
    sid: [u8; 32],
    mask: u8,
    seed: Zeroizing<[u8; 32]>,
}

const SEED_MSG_VERSION: u8 = 2;
const SEED_MSG_BYTES: usize = 1 + 1 + 1 + 32 + 1 + 32;

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

    /// Encode for transport: `2 ‖ from ‖ to ‖ sid ‖ mask ‖ K_I`.
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
        let mut seed = Zeroizing::new([0u8; 32]);
        seed.copy_from_slice(&bytes[36..68]);
        Ok(Self {
            from: bytes[1],
            to: bytes[2],
            sid: bytes[3..35].try_into().map_err(|_| Error::DkgMalformedMessage)?,
            mask: bytes[35],
            seed,
        })
    }
}

/// A party after Round 0, waiting for every party's commitments.
pub struct DkgParty {
    cfg: DkgConfig,
    coin: [u8; 32],
    /// Subsets this party deals: mask → K_I
    keys: BTreeMap<u8, Zeroizing<[u8; 32]>>,
    msg: DkgRound0Msg,
}

/// A party after Round 1, waiting for coin shares and its subset keys.
pub struct DkgAwaitingCoins {
    cfg: DkgConfig,
    sid: [u8; 32],
    /// All key commitments, mask → H(K_I)
    key_commits: BTreeMap<u8, [u8; 32]>,
    /// Coin commitments by sender
    coin_commits: Vec<[u8; 32]>,
    my_coin: [u8; 32],
    keys: BTreeMap<u8, Zeroizing<[u8; 32]>>,
}

/// A party after Round 2, waiting for `t_I` commitments.
pub struct DkgAwaitingCommits {
    cfg: DkgConfig,
    sid: [u8; 32],
    coin: [u8; 32],
    rho: [u8; 32],
    key: Zeroizing<[u8; 32]>,
    shares: BTreeMap<u8, Share>,
    /// `t_I` of every subset this party belongs to or deals
    known_t: BTreeMap<u8, PolyVecK>,
    own_commits: Vec<(u8, [u8; 32])>,
}

/// A party after Round 3, waiting for reveals.
pub struct DkgAwaitingReveals {
    cfg: DkgConfig,
    sid: [u8; 32],
    coin: [u8; 32],
    rho: [u8; 32],
    key: Zeroizing<[u8; 32]>,
    shares: BTreeMap<u8, Share>,
    /// All commitments, mask → C_I
    commits: BTreeMap<u8, [u8; 32]>,
}

/// A party after Round 4, waiting for confirmations.
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

fn coin_commitment(cfg: &DkgConfig, from: u8, coin: &[u8; 32]) -> [u8; 32] {
    hash32(&[b"th-ml-dsa-dkg-coin-v3", &[cfg.t, cfg.n], &cfg.context, &[from], coin])
}

fn key_commitment(cfg: &DkgConfig, mask: u8, key: &[u8; 32]) -> [u8; 32] {
    hash32(&[b"th-ml-dsa-dkg-key-v3", &[cfg.t, cfg.n], &cfg.context, &[dealer(mask), mask], key])
}

fn commitment(sid: &[u8; 32], mask: u8, t_packed: &[u8]) -> [u8; 32] {
    hash32(&[b"th-ml-dsa-dkg-commit-v3", sid, &[dealer(mask), mask], t_packed])
}

/// `s_I` from the joint coin and the subset key.
fn subset_share(sid: &[u8; 32], coin: &[u8; 32], mask: u8, key: &[u8; 32]) -> Share {
    let mut seed = Zeroizing::new([0u8; 64]);
    let mut h = Shake256::default();
    h.update(b"th-ml-dsa-dkg-subset-v3");
    h.update(sid);
    h.update(coin);
    h.update(&[mask]);
    h.update(key);
    h.finalize_xof().read(&mut seed[..]);
    share_from_seed(&seed)
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

/// Check that `froms` is exactly `0..n`, in any order.
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

/// Subsets dealt by `from`, in ascending order.
fn dealt_by(cfg: &DkgConfig, from: u8) -> Vec<u8> {
    subset_masks(cfg.n, cfg.t).into_iter().filter(|&k| dealer(k) == from).collect()
}

impl DkgParty {
    /// Round 0: validate the configuration, draw this party's coin share and
    /// subset keys, and commit to them.
    pub fn start<R: RngCore + CryptoRng>(cfg: DkgConfig, rng: &mut R) -> Result<(Self, DkgRound0Msg), Error> {
        get_threshold_params(cfg.t, cfg.n).ok_or(Error::InvalidParameters)?;
        if cfg.me >= cfg.n {
            return Err(Error::InvalidParameters);
        }
        let mut coin = [0u8; 32];
        rng.fill_bytes(&mut coin);
        let mut keys = BTreeMap::new();
        let mut key_commits = Vec::new();
        for mask in dealt_by(&cfg, cfg.me) {
            let mut k = Zeroizing::new([0u8; 32]);
            rng.fill_bytes(&mut k[..]);
            key_commits.push((mask, key_commitment(&cfg, mask, &k)));
            keys.insert(mask, k);
        }
        let msg = DkgRound0Msg { from: cfg.me, coin_commit: coin_commitment(&cfg, cfg.me, &coin), key_commits };
        Ok((Self { cfg, coin, keys, msg: msg.clone() }, msg))
    }

    /// Round 1: derive the session ID, reveal the coin share and send each
    /// subset key to the other members of the subset.
    ///
    /// Returns the coin broadcast and the private key messages, one per
    /// (subset, other member).
    pub fn round1(self, r0: &[DkgRound0Msg]) -> Result<(DkgAwaitingCoins, DkgCoinMsg, Vec<DkgSeedMsg>), Error> {
        let cfg = self.cfg;
        check_senders(r0.iter().map(|m| m.from), cfg.n)?;
        let mut sorted: Vec<&DkgRound0Msg> = r0.iter().collect();
        sorted.sort_by_key(|m| m.from);

        let mut h = Shake256::default();
        h.update(b"th-ml-dsa-dkg-sid-v3");
        h.update(&[cfg.t, cfg.n]);
        h.update(&cfg.context);
        let mut key_commits = BTreeMap::new();
        let mut coin_commits = Vec::new();
        for m in &sorted {
            if m.from == cfg.me && **m != self.msg {
                return Err(Error::DkgInconsistentTranscript);
            }
            let got: Vec<u8> = m.key_commits.iter().map(|c| c.0).collect();
            if got != dealt_by(&cfg, m.from) {
                return Err(Error::DkgMalformedMessage);
            }
            h.update(&[m.from]);
            h.update(&m.coin_commit);
            coin_commits.push(m.coin_commit);
            for (mask, c) in &m.key_commits {
                h.update(&[*mask]);
                h.update(c);
                key_commits.insert(*mask, *c);
            }
        }
        let mut sid = [0u8; 32];
        h.finalize_xof().read(&mut sid);

        let mut seed_msgs = Vec::new();
        for (&mask, k) in &self.keys {
            for to in 0..cfg.n {
                if to != cfg.me && mask & (1 << to) != 0 {
                    seed_msgs.push(DkgSeedMsg { from: cfg.me, to, sid, mask, seed: k.clone() });
                }
            }
        }

        let msg = DkgCoinMsg { from: cfg.me, sid, coin: self.coin };
        Ok((
            DkgAwaitingCoins { cfg, sid, key_commits, coin_commits, my_coin: self.coin, keys: self.keys },
            msg,
            seed_msgs,
        ))
    }
}

impl DkgAwaitingCoins {
    /// Round 2: check the coin shares and the keys received from other
    /// dealers, derive the subset secrets and commit to this party's `t_I`.
    pub fn round2<R: RngCore + CryptoRng>(
        mut self,
        coins: &[DkgCoinMsg],
        seeds_for_me: Vec<DkgSeedMsg>,
        rng: &mut R,
    ) -> Result<(DkgAwaitingCommits, DkgCommitMsg), Error> {
        let cfg = self.cfg;
        check_senders(coins.iter().map(|m| m.from), cfg.n)?;
        let mut shares_in = alloc::vec![[0u8; 32]; cfg.n as usize];
        for m in coins {
            if !bool::from(m.sid.ct_eq(&self.sid)) {
                return Err(Error::DkgMalformedMessage);
            }
            if m.from == cfg.me && m.coin != self.my_coin {
                return Err(Error::DkgInconsistentTranscript);
            }
            let c = coin_commitment(&cfg, m.from, &m.coin);
            if !bool::from(c.ct_eq(&self.coin_commits[m.from as usize])) {
                return Err(Error::DkgCommitmentMismatch { party: m.from });
            }
            shares_in[m.from as usize] = m.coin;
        }
        let mut h = Shake256::default();
        h.update(b"th-ml-dsa-dkg-coin-v3");
        h.update(&self.sid);
        for c in &shares_in {
            h.update(c);
        }
        let mut coin = [0u8; 32];
        h.finalize_xof().read(&mut coin);
        let rho = hash32(&[b"th-ml-dsa-dkg-rho-v3", &coin]);

        // Exactly one key for each subset containing me that I do not deal.
        let masks = subset_masks(cfg.n, cfg.t);
        let wanted: Vec<u8> = masks
            .iter()
            .copied()
            .filter(|&k| k & (1 << cfg.me) != 0 && dealer(k) != cfg.me)
            .collect();
        if seeds_for_me.len() != wanted.len() {
            return Err(Error::DkgMalformedMessage);
        }
        let mut shares = BTreeMap::new();
        let mut known_t = BTreeMap::new();
        for msg in seeds_for_me {
            if msg.to != cfg.me
                || !wanted.contains(&msg.mask)
                || msg.from != dealer(msg.mask)
                || !bool::from(msg.sid.ct_eq(&self.sid))
                || shares.contains_key(&msg.mask)
            {
                return Err(Error::DkgMalformedMessage);
            }
            let c = key_commitment(&cfg, msg.mask, &msg.seed);
            if !bool::from(c.ct_eq(&self.key_commits[&msg.mask])) {
                return Err(Error::DkgSeedMismatch { dealer: msg.from, mask: msg.mask });
            }
            let share = subset_share(&self.sid, &coin, msg.mask, &msg.seed);
            known_t.insert(msg.mask, compute_t(&rho, &share.s1h, &share.s2));
            shares.insert(msg.mask, share);
        }

        let mut own_commits = Vec::new();
        for (mask, k) in core::mem::take(&mut self.keys) {
            let share = subset_share(&self.sid, &coin, mask, &k);
            let t_i = compute_t(&rho, &share.s1h, &share.s2);
            own_commits.push((mask, commitment(&self.sid, mask, &pack_w_single(&t_i))));
            known_t.insert(mask, t_i);
            if mask & (1 << cfg.me) != 0 {
                shares.insert(mask, share);
            }
        }

        let mut key = Zeroizing::new([0u8; 32]);
        rng.fill_bytes(&mut key[..]);

        let msg = DkgCommitMsg { from: cfg.me, sid: self.sid, commits: own_commits.clone() };
        Ok((
            DkgAwaitingCommits { cfg, sid: self.sid, coin, rho, key, shares, known_t, own_commits },
            msg,
        ))
    }
}

impl DkgAwaitingCommits {
    /// Round 3: check all commitments (against this party's own `t_I` where
    /// it is a member), then reveal this party's `t_I`.
    pub fn round3(self, commits: &[DkgCommitMsg]) -> Result<(DkgAwaitingReveals, DkgRevealMsg), Error> {
        let cfg = self.cfg;
        check_senders(commits.iter().map(|m| m.from), cfg.n)?;

        let mut all = BTreeMap::new();
        for m in commits {
            if !bool::from(m.sid.ct_eq(&self.sid)) {
                return Err(Error::DkgMalformedMessage);
            }
            let got: Vec<u8> = m.commits.iter().map(|c| c.0).collect();
            if got != dealt_by(&cfg, m.from) {
                return Err(Error::DkgMalformedMessage);
            }
            if m.from == cfg.me && m.commits != self.own_commits {
                return Err(Error::DkgInconsistentTranscript);
            }
            for &(mask, c) in &m.commits {
                if let Some(t_i) = self.known_t.get(&mask) {
                    let mine = commitment(&self.sid, mask, &pack_w_single(t_i));
                    if !bool::from(mine.ct_eq(&c)) {
                        return Err(Error::DkgCommitmentMismatch { party: m.from });
                    }
                }
                all.insert(mask, c);
            }
        }

        let t_parts = self
            .known_t
            .iter()
            .filter(|(&mask, _)| dealer(mask) == cfg.me)
            .map(|(&mask, t_i)| (mask, pack_w_single(t_i)))
            .collect();

        let reveal = DkgRevealMsg { from: cfg.me, sid: self.sid, t_parts };
        Ok((
            DkgAwaitingReveals {
                cfg,
                sid: self.sid,
                coin: self.coin,
                rho: self.rho,
                key: self.key,
                shares: self.shares,
                commits: all,
            },
            reveal,
        ))
    }
}

impl DkgAwaitingReveals {
    /// Round 4: check every reveal against its commitment, compute the
    /// public key and this party's confirmation digest.
    pub fn round4(self, reveals: &[DkgRevealMsg]) -> Result<(DkgAwaitingConfirm, DkgConfirmMsg), Error> {
        let cfg = self.cfg;
        check_senders(reveals.iter().map(|m| m.from), cfg.n)?;

        let mut t = PolyVecK::zero();
        for m in reveals {
            if !bool::from(m.sid.ct_eq(&self.sid)) {
                return Err(Error::DkgMalformedMessage);
            }
            let got: Vec<u8> = m.t_parts.iter().map(|p| p.0).collect();
            if got != dealt_by(&cfg, m.from) {
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
        h.update(b"th-ml-dsa-dkg-confirm-v3");
        h.update(&self.sid);
        h.update(&self.coin);
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
    pub coins: Vec<DkgCoinMsg>,
    /// Round 1 private messages
    pub seeds: Vec<DkgSeedMsg>,
    /// Round 2 broadcasts
    pub commits: Vec<DkgCommitMsg>,
    /// Round 3 broadcasts
    pub reveals: Vec<DkgRevealMsg>,
    /// Round 4 broadcasts
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
    let mut coins = Vec::new();
    let mut seeds = Vec::new();
    for p in parties {
        let (s, c, sm) = p.round1(&round0)?;
        st1.push(s);
        coins.push(c);
        seeds.extend(sm);
    }
    let mut st2 = Vec::new();
    let mut commits = Vec::new();
    for (me, s) in st1.into_iter().enumerate() {
        let mine: Vec<DkgSeedMsg> = seeds
            .iter()
            .filter(|m| m.to as usize == me)
            .map(|m| DkgSeedMsg { from: m.from, to: m.to, sid: m.sid, mask: m.mask, seed: m.seed.clone() })
            .collect();
        let (s2, c) = s.round2(&coins, mine, rng)?;
        st2.push(s2);
        commits.push(c);
    }
    let mut st3 = Vec::new();
    let mut reveals = Vec::new();
    for s in st2 {
        let (s3, r) = s.round3(&commits)?;
        st3.push(s3);
        reveals.push(r);
    }
    let mut st4 = Vec::new();
    let mut confirms = Vec::new();
    for s in st3 {
        let (s4, c) = s.round4(&reveals)?;
        st4.push(s4);
        confirms.push(c);
    }
    let mut pk = [0u8; PK_BYTES];
    let mut sks = Vec::new();
    for s in st4 {
        let (p, sk) = s.finalize(&confirms)?;
        pk = p;
        sks.push(sk);
    }
    Ok((pk, sks, DkgTranscript { round0, coins, seeds, commits, reveals, confirms }))
}
