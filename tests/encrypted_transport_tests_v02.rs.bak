//! Integration tests for encrypted coordinator-to-party transport in threshold ML-DSA.
//!
//! Scenario covered in this module:
//! - 4 signing parties, threshold 3
//! - each party has its own ML-KEM key pair for transport encryption
//! - coordinator is provisioned with all 4 transport public keys
//! - coordinator sends Round-2 challenge messages encrypted per recipient

use std::collections::HashSet;

use chacha20poly1305::aead::{Aead, KeyInit, Payload};
use chacha20poly1305::{ChaCha20Poly1305, Nonce};
use kyber::{safe_encaps_derand, MlKemCiphertext, MlKemKeyPair, ML_KEM_768};
use rand::rngs::StdRng;
use rand::{RngCore, SeedableRng};
use sha3::{Digest, Sha3_256};

use threshold_ml_dsa::coordinator;
use threshold_ml_dsa::error::Error;
use threshold_ml_dsa::params::*;
use threshold_ml_dsa::poly::{MatrixA, PolyVecK, PolyVecL};
use threshold_ml_dsa::rss;
use threshold_ml_dsa::sign::Party;
use threshold_ml_dsa::verify;

const N_PARTIES: usize = 4;
const THRESHOLD: usize = 3;
const ROUND_CHALLENGE: u8 = 2;

#[derive(Clone)]
struct PlayerTransportIdentity {
    id: usize,
    keypair: MlKemKeyPair,
}

impl PlayerTransportIdentity {
    fn public_key(&self) -> &[u8] {
        self.keypair.public_key()
    }
}

#[derive(Clone)]
struct EncryptedEnvelope {
    receiver_id: usize,
    round: u8,
    counter: u64,
    session_id: [u8; 32],
    kem_ciphertext: Vec<u8>,
    nonce: [u8; 12],
    aad: Vec<u8>,
    ciphertext: Vec<u8>,
}

#[derive(Default)]
struct CoordinatorTransport {
    player_public_keys: Vec<Vec<u8>>,
    next_counter: u64,
}

impl CoordinatorTransport {
    fn new(players: &[PlayerTransportIdentity]) -> Self {
        let player_public_keys = players.iter().map(|p| p.public_key().to_vec()).collect();
        Self {
            player_public_keys,
            next_counter: 0,
        }
    }

    fn from_public_keys(player_public_keys: Vec<Vec<u8>>) -> Self {
        Self {
            player_public_keys,
            next_counter: 0,
        }
    }

    fn player_public_keys(&self) -> &[Vec<u8>] {
        &self.player_public_keys
    }

    fn encrypt_for_player<R: RngCore>(
        &mut self,
        receiver_id: usize,
        round: u8,
        session_id: [u8; 32],
        payload: &[u8],
        rng: &mut R,
    ) -> Result<EncryptedEnvelope, String> {
        if receiver_id >= self.player_public_keys.len() {
            return Err("receiver_id out of bounds".to_string());
        }

        let mut coins = [0u8; 32];
        rng.fill_bytes(&mut coins);
        let pk = &self.player_public_keys[receiver_id];
        let (kem_ct, shared_secret) =
            safe_encaps_derand(ML_KEM_768, pk, &coins).map_err(|e| format!("{e}"))?;

        let counter = self.next_counter;
        self.next_counter = self.next_counter.wrapping_add(1);

        let aad = build_aad(&session_id, round, receiver_id, counter);
        let nonce = build_nonce(receiver_id, counter);
        let key = derive_aead_key(
            shared_secret.as_bytes(),
            &session_id,
            round,
            receiver_id,
            counter,
        );

        let cipher = ChaCha20Poly1305::new_from_slice(&key).map_err(|e| format!("{e}"))?;
        let ciphertext = cipher
            .encrypt(
                Nonce::from_slice(&nonce),
                Payload {
                    msg: payload,
                    aad: &aad,
                },
            )
            .map_err(|_| "aead encrypt failed".to_string())?;

        Ok(EncryptedEnvelope {
            receiver_id,
            round,
            counter,
            session_id,
            kem_ciphertext: kem_ct.as_bytes().to_vec(),
            nonce,
            aad,
            ciphertext,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct ReplayToken {
    receiver_id: usize,
    round: u8,
    counter: u64,
    session_id: [u8; 32],
}

#[derive(Default)]
struct ReplayGuard {
    seen: HashSet<ReplayToken>,
}

impl ReplayGuard {
    fn accept_once(&mut self, env: &EncryptedEnvelope) -> bool {
        let token = ReplayToken {
            receiver_id: env.receiver_id,
            round: env.round,
            counter: env.counter,
            session_id: env.session_id,
        };
        self.seen.insert(token)
    }
}

fn decrypt_for_player(
    player: &PlayerTransportIdentity,
    envelope: &EncryptedEnvelope,
) -> Result<Vec<u8>, String> {
    if envelope.receiver_id != player.id {
        return Err("envelope receiver_id mismatch".to_string());
    }

    let ct = MlKemCiphertext::from_bytes(envelope.kem_ciphertext.clone());
    let shared_secret = player
        .keypair
        .decaps(&ct)
        .map_err(|e| format!("decaps failed: {e}"))?;

    let key = derive_aead_key(
        shared_secret.as_bytes(),
        &envelope.session_id,
        envelope.round,
        envelope.receiver_id,
        envelope.counter,
    );
    let cipher = ChaCha20Poly1305::new_from_slice(&key).map_err(|e| format!("{e}"))?;

    cipher
        .decrypt(
            Nonce::from_slice(&envelope.nonce),
            Payload {
                msg: &envelope.ciphertext,
                aad: &envelope.aad,
            },
        )
        .map_err(|_| "aead decrypt failed".to_string())
}

fn decrypt_for_player_with_replay_guard(
    player: &PlayerTransportIdentity,
    envelope: &EncryptedEnvelope,
    replay_guard: &mut ReplayGuard,
) -> Result<Vec<u8>, String> {
    if !replay_guard.accept_once(envelope) {
        return Err("replay detected".to_string());
    }
    decrypt_for_player(player, envelope)
}

fn build_aad(session_id: &[u8; 32], round: u8, receiver_id: usize, counter: u64) -> Vec<u8> {
    let mut aad = Vec::with_capacity(32 + 1 + 8 + 8 + 8);
    aad.extend_from_slice(b"TMLDSA-V1");
    aad.extend_from_slice(session_id);
    aad.push(round);
    aad.extend_from_slice(&(receiver_id as u64).to_le_bytes());
    aad.extend_from_slice(&counter.to_le_bytes());
    aad
}

fn build_nonce(receiver_id: usize, counter: u64) -> [u8; 12] {
    let mut nonce = [0u8; 12];
    nonce[..4].copy_from_slice(&(receiver_id as u32).to_le_bytes());
    nonce[4..].copy_from_slice(&counter.to_le_bytes());
    nonce
}

fn derive_aead_key(
    shared_secret: &[u8; 32],
    session_id: &[u8; 32],
    round: u8,
    receiver_id: usize,
    counter: u64,
) -> [u8; 32] {
    let mut h = Sha3_256::new();
    h.update(b"threshold-ml-dsa/kyber-chacha20poly1305/v1");
    h.update(shared_secret);
    h.update(session_id);
    h.update([round]);
    h.update((receiver_id as u64).to_le_bytes());
    h.update(counter.to_le_bytes());

    let digest = h.finalize();
    let mut key = [0u8; 32];
    key.copy_from_slice(&digest);
    key
}

fn key_material_from_seed(
    seed: [u8; 32],
) -> (Vec<u8>, [u8; TRBYTES], [u8; SEEDBYTES], PolyVecL, PolyVecK) {
    let (pk, sk) = verify::keygen(&seed);
    let mode = dilithium::params::ML_DSA_44;

    let mut rho_ref = [0u8; dilithium::params::SEEDBYTES];
    let mut tr_ref = [0u8; dilithium::params::TRBYTES];
    let mut key = [0u8; dilithium::params::SEEDBYTES];
    let mut t0 = dilithium::polyvec::PolyVecK::default();
    let mut s1_ref = dilithium::polyvec::PolyVecL::default();
    let mut s2_ref = dilithium::polyvec::PolyVecK::default();
    dilithium::packing::unpack_sk(
        mode,
        &mut rho_ref,
        &mut tr_ref,
        &mut key,
        &mut t0,
        &mut s1_ref,
        &mut s2_ref,
        &sk,
    );

    let mut s1 = PolyVecL::zero();
    let mut s2 = PolyVecK::zero();
    for i in 0..L {
        s1.polys[i].coeffs = s1_ref.vec[i].coeffs;
    }
    for i in 0..K {
        s2.polys[i].coeffs = s2_ref.vec[i].coeffs;
    }

    (pk, tr_ref, rho_ref, s1, s2)
}

fn compute_session_id(
    tr: &[u8; TRBYTES],
    msg: &[u8],
    active_party_ids: &[usize],
    attempt: usize,
) -> [u8; 32] {
    let mut h = Sha3_256::new();
    h.update(b"threshold-ml-dsa/encrypted-transport/session");
    h.update(tr);
    h.update(msg);
    h.update((attempt as u64).to_le_bytes());
    h.update((active_party_ids.len() as u64).to_le_bytes());
    for &id in active_party_ids {
        h.update((id as u64).to_le_bytes());
    }
    let digest = h.finalize();
    let mut out = [0u8; 32];
    out.copy_from_slice(&digest);
    out
}

fn generate_player_transport_identities(n: usize) -> Vec<PlayerTransportIdentity> {
    let mut out = Vec::with_capacity(n);
    for id in 0..n {
        let mut coins = [0u8; 64];

        let mut h0 = Sha3_256::new();
        h0.update(b"threshold-ml-dsa/player-transport-keygen/0");
        h0.update((id as u64).to_le_bytes());
        let d0 = h0.finalize();
        coins[..32].copy_from_slice(&d0);

        let mut h1 = Sha3_256::new();
        h1.update(b"threshold-ml-dsa/player-transport-keygen/1");
        h1.update((id as u64).to_le_bytes());
        let d1 = h1.finalize();
        coins[32..].copy_from_slice(&d1);

        out.push(PlayerTransportIdentity {
            id,
            keypair: MlKemKeyPair::generate_derand(ML_KEM_768, &coins),
        });
    }
    out
}

#[test]
fn test_encrypted_challenge_transport_n4_t3_fail_closed() {
    let mut rng = StdRng::seed_from_u64(0xE7C2_4A11_99);

    let (pk, tr, rho, s1, s2) = key_material_from_seed([7u8; 32]);
    let shares = rss::distribute_key(&s1, &s2, N_PARTIES, THRESHOLD, &mut rng).unwrap();
    let a_hat = MatrixA::expand(&rho);
    let mut parties: Vec<Party> = shares
        .iter()
        .map(|share| Party::new(share, a_hat.clone()))
        .collect();

    let players = generate_player_transport_identities(N_PARTIES);
    let mut transport = CoordinatorTransport::new(&players);
    assert_eq!(
        transport.player_public_keys().len(),
        N_PARTIES,
        "coordinator should hold one transport pk per player"
    );

    let active_indices = [0usize, 1usize, 2usize];
    let active_party_ids: Vec<usize> = active_indices.iter().map(|&idx| parties[idx].id).collect();
    let msg = b"encrypted challenge transport integration test";

    let mut saw_decrypted_challenge = false;
    let mut saw_safe_abort = false;
    let mut saw_success = false;

    for attempt in 0..128usize {
        let session_id = compute_session_id(&tr, msg, &active_party_ids, attempt);

        let mut hashes = Vec::with_capacity(THRESHOLD);
        for &idx in &active_indices {
            match parties[idx].precommit(&mut rng, &active_party_ids, &session_id) {
                Ok(hash) => hashes.push(hash),
                Err(_) => {
                    saw_safe_abort = true;
                    hashes.clear();
                    break;
                }
            }
        }
        if hashes.len() != THRESHOLD {
            continue;
        }

        let commitments = active_indices
            .iter()
            .filter_map(|&idx| parties[idx].reveal())
            .collect::<Vec<_>>();
        if commitments.len() != THRESHOLD {
            saw_safe_abort = true;
            continue;
        }

        let w = match coordinator::aggregate_commitments(&commitments, &hashes) {
            Ok(w) => w,
            Err(_) => {
                saw_safe_abort = true;
                continue;
            }
        };
        let c_tilde = coordinator::compute_challenge(&w, msg, &tr);

        let mut partials = Vec::with_capacity(THRESHOLD);
        let mut round_rejected = false;
        for &idx in &active_indices {
            let party_id = parties[idx].id;
            let envelope = transport
                .encrypt_for_player(party_id, ROUND_CHALLENGE, session_id, &c_tilde, &mut rng)
                .expect("challenge encryption should succeed");

            let decrypted_challenge = decrypt_for_player(&players[party_id], &envelope)
                .expect("challenge decryption should succeed");
            assert_eq!(
                decrypted_challenge, c_tilde,
                "decrypted challenge must match coordinator challenge"
            );
            saw_decrypted_challenge = true;

            match parties[idx].sign(&decrypted_challenge) {
                Ok(partial) => partials.push(partial),
                Err(Error::LocalRejectionAbort) => {
                    saw_safe_abort = true;
                    round_rejected = true;
                    break;
                }
                Err(e) => panic!("unexpected party signing error: {e}"),
            }
        }
        if round_rejected || partials.len() != THRESHOLD {
            continue;
        }

        match coordinator::aggregate_responses(&partials, &w, &c_tilde, THRESHOLD, &pk) {
            Ok(sig) => {
                let sig_bytes = sig.to_bytes();
                assert!(
                    verify::verify(&sig_bytes, msg, &pk),
                    "aggregated signature from encrypted flow must verify"
                );
                saw_success = true;
                break;
            }
            Err(Error::InvalidSignature) | Err(Error::InsufficientResponses) => {
                saw_safe_abort = true;
            }
            Err(e) => panic!("unexpected aggregation error: {e}"),
        }
    }

    assert!(
        saw_decrypted_challenge,
        "expected at least one successful encrypted challenge delivery"
    );
    assert!(
        saw_success || saw_safe_abort,
        "expected either a successful signature or a fail-closed abort"
    );
}

#[test]
fn test_encrypted_challenge_tamper_rejected() {
    let mut rng = StdRng::seed_from_u64(0xA1B2_C3D4);
    let players = generate_player_transport_identities(N_PARTIES);
    let mut transport = CoordinatorTransport::new(&players);

    let session_id = [0x55u8; 32];
    let challenge = [0xA5u8; CTILDEBYTES];
    let mut envelope = transport
        .encrypt_for_player(0, ROUND_CHALLENGE, session_id, &challenge, &mut rng)
        .expect("challenge encryption should succeed");

    envelope.ciphertext[0] ^= 0x80;

    let err = decrypt_for_player(&players[0], &envelope).unwrap_err();
    assert!(
        err.contains("decrypt failed"),
        "tampered ciphertext must fail authentication"
    );
}

#[test]
fn test_encrypted_challenge_wrong_public_key_mapping_fails() {
    let mut rng = StdRng::seed_from_u64(0x77_00_33_11);
    let players = generate_player_transport_identities(N_PARTIES);

    // Misconfigure coordinator mapping: swap pk[0] and pk[1].
    let mut mapped_pks: Vec<Vec<u8>> = players.iter().map(|p| p.public_key().to_vec()).collect();
    mapped_pks.swap(0, 1);
    let mut transport = CoordinatorTransport::from_public_keys(mapped_pks);

    let session_id = [0x33u8; 32];
    let challenge = [0x11u8; CTILDEBYTES];
    let envelope = transport
        .encrypt_for_player(0, ROUND_CHALLENGE, session_id, &challenge, &mut rng)
        .expect("challenge encryption should succeed");

    assert!(
        decrypt_for_player(&players[0], &envelope).is_err(),
        "decrypt must fail when coordinator uses wrong recipient public key"
    );
}

#[test]
fn test_encrypted_challenge_replay_rejected() {
    let mut rng = StdRng::seed_from_u64(0x42_42_42_42);
    let players = generate_player_transport_identities(N_PARTIES);
    let mut transport = CoordinatorTransport::new(&players);

    let session_id = [0x99u8; 32];
    let challenge = [0x22u8; CTILDEBYTES];
    let envelope = transport
        .encrypt_for_player(2, ROUND_CHALLENGE, session_id, &challenge, &mut rng)
        .expect("challenge encryption should succeed");

    let mut replay_guard = ReplayGuard::default();
    let first = decrypt_for_player_with_replay_guard(&players[2], &envelope, &mut replay_guard)
        .expect("first delivery should decrypt");
    assert_eq!(first, challenge.to_vec());

    let second = decrypt_for_player_with_replay_guard(&players[2], &envelope, &mut replay_guard);
    assert!(
        second.err().map(|e| e.contains("replay")).unwrap_or(false),
        "second delivery of same envelope must be rejected as replay"
    );
}
