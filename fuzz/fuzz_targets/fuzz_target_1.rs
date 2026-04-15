#![no_main]

use libfuzzer_sys::fuzz_target;
use threshold_ml_dsa::params::get_threshold_params;
use threshold_ml_dsa::partition;
use threshold_ml_dsa::rss;

fuzz_target!(|data: &[u8]| {
    if data.is_empty() {
        return;
    }

    // Build a deterministic seed from fuzz input.
    let mut seed = [0u8; 32];
    for (i, b) in seed.iter_mut().enumerate() {
        *b = data[i % data.len()];
    }

    // Map bytes to a valid (t, n) with 2 <= t <= n <= 6.
    let t = 2 + (data[0] % 5); // 2..=6
    let n = t + (data[data.len() / 2] % (7 - t)); // t..=6

    if let Some(params) = get_threshold_params(t, n) {
        let _ = rss::keygen_from_seed(&seed, &params);

        // Canonical active set and a malformed one (for error-path coverage).
        let active: Vec<u8> = (0..t).collect();
        let _ = partition::rss_recover(&active, n, t);

        if t >= 2 {
            let mut malformed = active.clone();
            malformed[0] = malformed[1]; // duplicate id
            let _ = partition::rss_recover(&malformed, n, t);
        }
    }
});
