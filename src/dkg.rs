//! Distributed Key Generation (DKG) for threshold ML-DSA-44.
//!
//! Produces the same replicated-secret-sharing key structure as
//! [`crate::rss::keygen_from_seed`] without any party learning the full
//! secret. It follows the structure of the paper's DKG (ePrint 2026/013 full
//! version, App. D, Fig. 13): the lowest member of each (N−T+1)-subset `I`
//! sends a key `K_I` privately to the other members, the parties
//! commit-and-reveal a joint coin, `s_I = H(I, K_I, coin)`, each
//! `t_I = A·s₁_I + s₂_I` is committed and then revealed, and
//! `t = Σ t_I` goes through Power2Round.
//!
//! Differences from Fig. 13 (none weakens it; encodings are this crate's own):
//! - `K_I` is also committed in round 0. The paper binds it only by
//!   delivery timing; the commitment catches a bad key one round earlier.
//! - A session ID (hash of all round-0 messages) and the application
//!   `context` are bound into every later message and hash.
//! - Transcript agreement uses a separate round of unsigned digests over an
//!   authenticated broadcast. The paper instead signs the transcript with
//!   long-term keys (a PKI); to get that, sign
//!   [`DkgAwaitingConfirm::transcript_digest`] with your identity keys.
//! - `tr` is the 64-byte FIPS 204 `H(pk)` (the paper writes a 256-bit
//!   `H(ρ‖t₁)`), which standard ML-DSA verification requires.
//! - Five message rounds instead of four.
//!
//! The paper's DKG security argument is informal (App. D.3); its
//! unforgeability theorem assumes the trusted dealer of Fig. 5.
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
//!   ML-DSA public key. Fig. 13 also publishes every `t_I`; the proof's
//!   key-replacement step handles one honest `t_I`, and a standard
//!   multi-sample MLWE hybrid covers the rest.
//! - `K_I` and `c_i` are committed before the coin is revealed, so a
//!   dealer cannot grind `s_I` for a subset that has an honest member, and
//!   no party can bias `ρ` other than by aborting.
//! - Commit-then-reveal of `t_I` stops a coalition that holds a whole
//!   subset from choosing its `t_I` after seeing the others (rogue key). It
//!   does not force such a `t_I` to be well-formed: the coalition can commit
//!   to any value, which only makes the key unusable (a denial of service it
//!   already has), since the all-honest subset's `t_M` still hides the key.
//! - The confirm round detects equivocation in the broadcast.
//! - A party can still abort after seeing the coin or the reveals
//!   (attributable over an authenticated broadcast). Each abort-and-restart
//!   lets it pick among candidate keys (a log₂(restarts)-bit choice, which
//!   does not weaken ML-DSA); restarting without the party means a new
//!   party set and a new (T, N). (The paper's "cannot bias" claim does not
//!   consider abort-and-restart.)
//! - Termination is not unanimous: a party can send a correct digest to some
//!   honest parties and a wrong one to others, so some finalize and some
//!   abort. Do not use the public key until every party has confirmed
//!   success out of band.
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
        let mut sid = [0u8; 32];
        sid.copy_from_slice(&bytes[3..35]);
        Ok(Self {
            from: bytes[1],
            to: bytes[2],
            sid,
            mask: bytes[35],
            seed,
        })
    }
}

/// A party after Round 0, waiting for every party's commitments.
pub struct DkgParty {
    cfg: DkgConfig,
    /// Secret until revealed in round 1 (an early leak lets others bias ρ).
    coin: Zeroizing<[u8; 32]>,
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
fn check_senders(froms: &mut dyn Iterator<Item = u8>, n: u8) -> Result<(), Error> {
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
        let mut coin = Zeroizing::new([0u8; 32]);
        rng.fill_bytes(&mut coin[..]);
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
        check_senders(&mut r0.iter().map(|m| m.from), cfg.n)?;
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

        let msg = DkgCoinMsg { from: cfg.me, sid, coin: *self.coin };
        Ok((
            DkgAwaitingCoins { cfg, sid, key_commits, coin_commits, my_coin: *self.coin, keys: self.keys },
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
        check_senders(&mut coins.iter().map(|m| m.from), cfg.n)?;
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
        h.update(b"th-ml-dsa-dkg-coin-derive-v3");
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
        check_senders(&mut commits.iter().map(|m| m.from), cfg.n)?;

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
        check_senders(&mut reveals.iter().map(|m| m.from), cfg.n)?;

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
    /// Digest of this party's view of the transcript (session ID, coin,
    /// public key, all `t_I` commitments), the value exchanged in the
    /// confirm round. Applications with long-term identity keys can sign it
    /// before broadcasting, as the paper's Fig. 13 does.
    #[must_use]
    pub fn transcript_digest(&self) -> [u8; 32] {
        self.digest
    }

    /// Final step: output the key only if every party saw the same transcript.
    pub fn finalize(mut self, confirms: &[DkgConfirmMsg]) -> Result<([u8; PK_BYTES], ThresholdPrivateKey), Error> {
        check_senders(&mut confirms.iter().map(|m| m.from), self.cfg.n)?;
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

#[cfg(test)]
mod coverage_tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    const T: u8 = 2;
    const N: u8 = 3;

    fn copy(m: &DkgSeedMsg) -> DkgSeedMsg {
        DkgSeedMsg { from: m.from, to: m.to, sid: m.sid, mask: m.mask, seed: m.seed.clone() }
    }

    fn until_round1(rng: &mut StdRng) -> (Vec<DkgAwaitingCoins>, Vec<DkgCoinMsg>, Vec<DkgSeedMsg>) {
        let mut parties = Vec::new();
        let mut r0 = Vec::new();
        for me in 0..N {
            let (p, m) = DkgParty::start(DkgConfig { t: T, n: N, me, context: [0; 32] }, rng).unwrap();
            parties.push(p);
            r0.push(m);
        }
        let mut st = Vec::new();
        let mut coins = Vec::new();
        let mut seeds = Vec::new();
        for p in parties {
            let (s, c, sm) = p.round1(&r0).unwrap();
            st.push(s);
            coins.push(c);
            seeds.extend(sm);
        }
        (st, coins, seeds)
    }

    fn mine(seeds: &[DkgSeedMsg], me: u8) -> Vec<DkgSeedMsg> {
        seeds.iter().filter(|m| m.to == me).map(copy).collect()
    }

    fn until_round2(rng: &mut StdRng) -> (Vec<DkgAwaitingCommits>, Vec<DkgCommitMsg>) {
        let (st, coins, seeds) = until_round1(rng);
        let mut out = Vec::new();
        let mut commits = Vec::new();
        for (me, s) in st.into_iter().enumerate() {
            let (s2, c) = s.round2(&coins, mine(&seeds, me as u8), rng).unwrap();
            out.push(s2);
            commits.push(c);
        }
        (out, commits)
    }

    fn until_round3(rng: &mut StdRng) -> (Vec<DkgAwaitingReveals>, Vec<DkgRevealMsg>) {
        let (st, commits) = until_round2(rng);
        let mut out = Vec::new();
        let mut reveals = Vec::new();
        for s in st {
            let (s3, r) = s.round3(&commits).unwrap();
            out.push(s3);
            reveals.push(r);
        }
        (out, reveals)
    }

    #[test]
    fn seed_message_accessors_and_parsing() {
        let mut rng = StdRng::seed_from_u64(1);
        let (_, _, seeds) = until_round1(&mut rng);
        let m = &seeds[0];
        let back = DkgSeedMsg::from_bytes(&m.to_bytes()).unwrap();
        assert_eq!((back.from(), back.to(), back.mask()), (m.from(), m.to(), m.mask()));
        let mut bad = m.to_bytes().to_vec();
        bad[0] = 1;
        assert!(DkgSeedMsg::from_bytes(&bad).is_err());
    }

    #[test]
    fn helpers_reject_malformed_input() {
        assert_eq!(unpack_t_canonical(&[0u8; 5]).err(), Some(Error::DkgMalformedMessage));
        let mut t = PolyVecK::zero();
        t.polys[0].coeffs[0] = Q - 1;
        let mut packed = pack_w_single(&t);
        assert!(unpack_t_canonical(&packed).is_ok());
        packed.iter_mut().take(3).for_each(|b| *b = 0xFF); // first coefficient = 2²³ − 1 ≥ q
        assert_eq!(unpack_t_canonical(&packed).err(), Some(Error::DkgMalformedMessage));
        assert!(check_senders(&mut [0u8, 1].into_iter(), 3).is_err()); // missing sender
        assert!(check_senders(&mut [0u8, 3, 1].into_iter(), 3).is_err()); // out of range
        assert!(check_senders(&mut [0u8, 1, 2].into_iter(), 3).is_ok());
    }

    #[test]
    fn round2_rejects_duplicate_seed() {
        let mut rng = StdRng::seed_from_u64(2);
        let (st, coins, seeds) = until_round1(&mut rng);
        // Party 2 needs two keys; send the same one twice.
        let m = mine(&seeds, 2);
        let dup = vec![copy(&m[0]), copy(&m[0])];
        let s2 = st.into_iter().nth(2).unwrap();
        assert_eq!(s2.round2(&coins, dup, &mut rng).err(), Some(Error::DkgMalformedMessage));
    }

    #[test]
    fn round3_rejects_malformed_commits() {
        let mut rng = StdRng::seed_from_u64(3);
        let (st, commits) = until_round2(&mut rng);
        let mut st = st.into_iter();
        let (s0, s1) = (st.next().unwrap(), st.next().unwrap());

        let mut bad_sid = commits.clone();
        bad_sid[1].sid[0] ^= 1;
        assert_eq!(s0.round3(&bad_sid).err(), Some(Error::DkgMalformedMessage));

        let mut bad_masks = commits.clone();
        bad_masks[1].commits.clear();
        assert_eq!(s1.round3(&bad_masks).err(), Some(Error::DkgMalformedMessage));

        let (st, commits) = until_round2(&mut rng);
        let s0 = st.into_iter().next().unwrap();
        let mut own = commits.clone();
        own[0].commits[1].1[0] ^= 1; // party 0's own commitment, altered in transit
        assert_eq!(s0.round3(&own).err(), Some(Error::DkgInconsistentTranscript));
    }

    #[test]
    fn round4_rejects_malformed_reveals() {
        let mut rng = StdRng::seed_from_u64(4);
        let (st, reveals) = until_round3(&mut rng);
        let mut st = st.into_iter();
        let (s0, s1) = (st.next().unwrap(), st.next().unwrap());

        let mut bad_sid = reveals.clone();
        bad_sid[1].sid[0] ^= 1;
        assert_eq!(s0.round4(&bad_sid).err(), Some(Error::DkgMalformedMessage));

        let mut bad_masks = reveals.clone();
        bad_masks[0].t_parts.pop();
        assert_eq!(s1.round4(&bad_masks).err(), Some(Error::DkgMalformedMessage));
    }
}

#[cfg(test)]
mod protocol_tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn cfg(me: u8) -> DkgConfig {
        DkgConfig { t: 2, n: 3, me, context: [5; 32] }
    }

    fn copy(m: &DkgSeedMsg) -> DkgSeedMsg {
        DkgSeedMsg { from: m.from, to: m.to, sid: m.sid, mask: m.mask, seed: m.seed.clone() }
    }

    fn start_all(rng: &mut StdRng) -> (Vec<DkgParty>, Vec<DkgRound0Msg>) {
        (0..3).map(|me| DkgParty::start(cfg(me), rng).unwrap()).unzip()
    }

    fn round1_all(rng: &mut StdRng) -> (Vec<DkgAwaitingCoins>, Vec<DkgCoinMsg>, Vec<DkgSeedMsg>) {
        let (parties, r0) = start_all(rng);
        let mut st = Vec::new();
        let mut coins = Vec::new();
        let mut seeds = Vec::new();
        for p in parties {
            let (s, c, sm) = p.round1(&r0).unwrap();
            st.push(s);
            coins.push(c);
            seeds.extend(sm);
        }
        (st, coins, seeds)
    }

    fn mine(seeds: &[DkgSeedMsg], me: u8) -> Vec<DkgSeedMsg> {
        seeds.iter().filter(|m| m.to == me).map(copy).collect()
    }

    #[test]
    fn simulate_produces_consistent_keys() {
        let mut rng = StdRng::seed_from_u64(10);
        let (pk, sks, transcript) = simulate(2, 3, [0; 32], &mut rng).unwrap();
        assert_eq!(sks.len(), 3);
        assert_eq!(transcript.confirms.len(), 3);
        let tr = tr_from_pk(&pk);
        for (i, sk) in sks.iter().enumerate() {
            assert_eq!(sk.id as usize, i);
            assert_eq!(sk.tr, tr);
            assert_eq!(sk.rho[..], pk[..32]);
        }
        assert!(simulate(3, 2, [0; 32], &mut rng).is_err());
    }

    #[test]
    fn start_and_round1_reject_bad_input() {
        let mut rng = StdRng::seed_from_u64(11);
        assert_eq!(DkgParty::start(cfg(3), &mut rng).err(), Some(Error::InvalidParameters));

        let (parties, r0) = start_all(&mut rng);
        let mut parties = parties.into_iter();
        let p0 = parties.next().unwrap();
        let p1 = parties.next().unwrap();
        let mut own = r0.clone();
        own[0].coin_commit[0] ^= 1;
        assert_eq!(p0.round1(&own).err(), Some(Error::DkgInconsistentTranscript));
        let mut masks = r0.clone();
        masks[0].key_commits.pop();
        assert_eq!(p1.round1(&masks).err(), Some(Error::DkgMalformedMessage));
    }

    #[test]
    fn round2_rejects_bad_coins_and_seeds() {
        let mut rng = StdRng::seed_from_u64(12);
        // wrong session in a coin message
        let (st, coins, seeds) = round1_all(&mut rng);
        let mut bad = coins.clone();
        bad[1].sid[0] ^= 1;
        let s0 = st.into_iter().next().unwrap();
        assert_eq!(s0.round2(&bad, mine(&seeds, 0), &mut rng).err(), Some(Error::DkgMalformedMessage));

        // own coin changed in transit
        let (st, coins, seeds) = round1_all(&mut rng);
        let mut bad = coins.clone();
        bad[0].coin[0] ^= 1;
        let s0 = st.into_iter().next().unwrap();
        assert_eq!(s0.round2(&bad, mine(&seeds, 0), &mut rng).err(), Some(Error::DkgInconsistentTranscript));

        // another party's coin does not open its commitment
        let (st, coins, seeds) = round1_all(&mut rng);
        let mut bad = coins.clone();
        bad[1].coin[0] ^= 1;
        let s0 = st.into_iter().next().unwrap();
        assert_eq!(s0.round2(&bad, mine(&seeds, 0), &mut rng).err(), Some(Error::DkgCommitmentMismatch { party: 1 }));

        // missing seed, and a seed that does not open its key commitment
        let (st, coins, seeds) = round1_all(&mut rng);
        let mut st = st.into_iter();
        let s1 = st.nth(1).unwrap();
        let s2 = st.next().unwrap();
        assert_eq!(s1.round2(&coins, Vec::new(), &mut rng).err(), Some(Error::DkgMalformedMessage));
        let mut wrong = mine(&seeds, 2);
        wrong[0].seed[0] ^= 1;
        let mask = wrong[0].mask;
        assert_eq!(
            s2.round2(&coins, wrong, &mut rng).err(),
            Some(Error::DkgSeedMismatch { dealer: dealer(mask), mask })
        );
    }

    #[test]
    fn later_rounds_detect_bad_commitments() {
        let mut rng = StdRng::seed_from_u64(13);
        let (st, coins, seeds) = round1_all(&mut rng);
        let mut st2 = Vec::new();
        let mut commits = Vec::new();
        for (me, s) in st.into_iter().enumerate() {
            let (s2, c) = s.round2(&coins, mine(&seeds, me as u8), &mut rng).unwrap();
            st2.push(s2);
            commits.push(c);
        }
        // Party 0's commitment for {0,1} is wrong: member 1 catches it in round 3.
        let mut forged = commits.clone();
        forged[0].commits[0].1[0] ^= 1;
        let mut st2 = st2.into_iter();
        let s0 = st2.next().unwrap();
        let s1 = st2.next().unwrap();
        let s2 = st2.next().unwrap();
        assert_eq!(s1.round3(&forged).err(), Some(Error::DkgCommitmentMismatch { party: 0 }));

        // Non-member 2 accepts the forged commitment, then the honest reveal
        // does not open it.
        let (s2, _) = s2.round3(&forged).unwrap();
        let (_, reveal0) = s0.round3(&commits).unwrap();
        let reveals = vec![reveal0, DkgRevealMsg { from: 1, sid: s2.sid, t_parts: Vec::new() }, DkgRevealMsg {
            from: 2,
            sid: s2.sid,
            t_parts: Vec::new(),
        }];
        assert_eq!(s2.round4(&reveals).err(), Some(Error::DkgCommitmentMismatch { party: 0 }));
    }

    #[test]
    fn finalize_requires_matching_confirmations() {
        let mut rng = StdRng::seed_from_u64(14);
        let (st, coins, seeds) = round1_all(&mut rng);
        let mut st2 = Vec::new();
        let mut commits = Vec::new();
        for (me, s) in st.into_iter().enumerate() {
            let (s2, c) = s.round2(&coins, mine(&seeds, me as u8), &mut rng).unwrap();
            st2.push(s2);
            commits.push(c);
        }
        let mut st3 = Vec::new();
        let mut reveals = Vec::new();
        for s in st2 {
            let (s3, r) = s.round3(&commits).unwrap();
            st3.push(s3);
            reveals.push(r);
        }
        let mut st4 = Vec::new();
        let mut confirms = Vec::new();
        for s in st3 {
            let (s4, c) = s.round4(&reveals).unwrap();
            st4.push(s4);
            confirms.push(c);
        }
        assert!(confirms.iter().all(|c| c.digest == confirms[0].digest));
        assert!(st4.iter().all(|s| s.transcript_digest() == confirms[0].digest));
        confirms[2].digest[0] ^= 1;
        let s0 = st4.into_iter().next().unwrap();
        assert_eq!(s0.finalize(&confirms).err(), Some(Error::DkgInconsistentTranscript));
    }
}
