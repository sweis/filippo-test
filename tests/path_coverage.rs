//! Diagnostic: count code-path coverage of the Wycheproof vectors against
//! the FIPS 204 verification algorithm. Always passes; output is on stderr.
//!
//! This duplicates verification logic so the production library stays
//! uninstrumented (which matters for mutation testing).

use mldsa_verify::{params::ParamSet, ML_DSA_44, ML_DSA_65, ML_DSA_87};
use serde::Deserialize;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128, Shake256,
};

const Q: i32 = 8_380_417;
const N: usize = 256;
const D: u32 = 13;

#[derive(Deserialize)]
struct TestFile {
    #[serde(rename = "testGroups")]
    test_groups: Vec<TestGroup>,
}
#[derive(Deserialize)]
struct TestGroup {
    #[serde(rename = "publicKey")]
    public_key: String,
    tests: Vec<TestCase>,
}
#[derive(Deserialize)]
struct TestCase {
    msg: String,
    #[serde(default)]
    ctx: Option<String>,
    sig: String,
    result: String,
}

#[derive(Default)]
struct Stats {
    z_bound_exceeded: u64,
    max_z_inf_norm: i32,
    hint_count_decreased: u64,
    decompose_edge: u64,
    decompose_edge_h1: u64,
    decompose_edge_r0_0_h1: u64,
    use_hint_r0_zero_h1: u64,
    use_hint_r0_pos_h1: u64,
    use_hint_r0_nonpos_h1: u64,
    rej_ntt_saw_q: u64,
}

fn ntt(p: &mut [i32; N]) {
    let mut k = 0;
    let mut len = 128;
    while len >= 1 {
        let mut start = 0;
        while start < N {
            k += 1;
            let zeta = ZETAS[k];
            for j in start..start + len {
                let t = (zeta as i64 * p[j + len] as i64 % Q as i64) as i32;
                let a = p[j];
                p[j + len] = (a - t).rem_euclid(Q);
                p[j] = (a + t).rem_euclid(Q);
            }
            start += 2 * len;
        }
        len /= 2;
    }
}

fn inv_ntt(p: &mut [i32; N]) {
    let n_inv: i64 = 8_347_681;
    let mut k = 256;
    let mut len = 1;
    while len < N {
        let mut start = 0;
        while start < N {
            k -= 1;
            let zeta = (Q - ZETAS[k]) as i64;
            for j in start..start + len {
                let t = p[j];
                p[j] = (t + p[j + len]).rem_euclid(Q);
                p[j + len] =
                    ((zeta * (t - p[j + len]) as i64) % Q as i64).rem_euclid(Q as i64) as i32;
            }
            start += 2 * len;
        }
        len *= 2;
    }
    for x in p.iter_mut() {
        *x = ((*x as i64 * n_inv) % Q as i64) as i32;
    }
}

fn expand_a(rho: &[u8; 32], k: usize, l: usize, stats: &mut Stats) -> Vec<Vec<[i32; N]>> {
    let mut rows = Vec::with_capacity(k);
    for r in 0..k {
        let mut row = Vec::with_capacity(l);
        for c in 0..l {
            let mut h = Shake128::default();
            h.update(rho);
            h.update(&[c as u8, r as u8]);
            let mut rdr = h.finalize_xof();
            let mut p = [0i32; N];
            let mut buf = [0u8; 3];
            let mut i = 0;
            while i < N {
                rdr.read(&mut buf);
                let v = (buf[0] as u32) | ((buf[1] as u32) << 8) | (((buf[2] & 0x7F) as u32) << 16);
                if v as i32 == Q {
                    stats.rej_ntt_saw_q += 1;
                }
                if (v as i32) < Q {
                    p[i] = v as i32;
                    i += 1;
                }
            }
            row.push(p);
        }
        rows.push(row);
    }
    rows
}

fn sample_in_ball(ct: &[u8], tau: usize) -> [i32; N] {
    let mut h = Shake256::default();
    h.update(ct);
    let mut rdr = h.finalize_xof();
    let mut sb = [0u8; 8];
    rdr.read(&mut sb);
    let mut signs = u64::from_le_bytes(sb);
    let mut c = [0i32; N];
    for i in (N - tau)..N {
        let mut jb = [0u8];
        loop {
            rdr.read(&mut jb);
            if (jb[0] as usize) <= i {
                break;
            }
        }
        let j = jb[0] as usize;
        c[i] = c[j];
        c[j] = if signs & 1 == 0 { 1 } else { Q - 1 };
        signs >>= 1;
    }
    c
}

fn decompose_edge(r: i32, g2: i32) -> (i32, i32, bool) {
    let rp = r.rem_euclid(Q);
    let tg2 = 2 * g2;
    let mut r0 = rp % tg2;
    if r0 > g2 {
        r0 -= tg2;
    }
    if rp - r0 == Q - 1 {
        (0, r0 - 1, true)
    } else {
        ((rp - r0) / tg2, r0, false)
    }
}

fn check_file(path: &str, ps: &ParamSet, stats: &mut Stats) {
    let data = std::fs::read_to_string(path).unwrap();
    let file: TestFile = serde_json::from_str(&data).unwrap();

    let k = ps.k;
    let l = ps.l;
    let ct_len = ps.c_tilde_bytes();
    let z_bits = 1 + (ps.gamma1 as u32).ilog2();
    let z_poly_bytes = (z_bits as usize) * N / 8;
    let omega = ps.omega;

    for tg in &file.test_groups {
        let pk = hex::decode(&tg.public_key).unwrap();
        if pk.len() != ps.pk_bytes() {
            continue;
        }
        for tc in &tg.tests {
            let sig = hex::decode(&tc.sig).unwrap();
            if sig.len() != ps.sig_bytes() {
                continue;
            }
            let ctx_bytes = tc.ctx.as_deref().map(|s| hex::decode(s).unwrap());
            let ctx: &[u8] = ctx_bytes.as_deref().unwrap_or(&[]);
            if ctx.len() > 255 {
                continue;
            }

            // --- z inf-norm ---
            let mut z = vec![[0i32; N]; l];
            let mut z_inf = 0i32;
            for i in 0..l {
                let off = ct_len + i * z_poly_bytes;
                let mask: u64 = (1 << z_bits) - 1;
                let mut acc: u64 = 0;
                let mut ab: u32 = 0;
                let mut bi = 0;
                for jj in 0..N {
                    while ab < z_bits {
                        acc |= (sig[off + bi] as u64) << ab;
                        bi += 1;
                        ab += 8;
                    }
                    let raw = (acc & mask) as i32;
                    acc >>= z_bits;
                    ab -= z_bits;
                    z[i][jj] = ps.gamma1 - raw;
                    if z[i][jj].abs() > z_inf {
                        z_inf = z[i][jj].abs();
                    }
                }
            }
            if z_inf > stats.max_z_inf_norm {
                stats.max_z_inf_norm = z_inf;
            }
            if z_inf >= ps.gamma1 - ps.beta {
                stats.z_bound_exceeded += 1;
            }

            // --- hints ---
            let hoff = ct_len + l * z_poly_bytes;
            let y = &sig[hoff..];
            let mut index = 0usize;
            let mut h = vec![[0i32; N]; k];
            let mut hint_ok = true;
            for i in 0..k {
                let end = y[omega + i] as usize;
                if end < index {
                    stats.hint_count_decreased += 1;
                    hint_ok = false;
                    break;
                }
                if end > omega {
                    hint_ok = false;
                    break;
                }
                let mut first = true;
                for j in index..end {
                    if !first && y[j] <= y[j - 1] {
                        hint_ok = false;
                        break;
                    }
                    first = false;
                    h[i][y[j] as usize] = 1;
                }
                if !hint_ok {
                    break;
                }
                index = end;
            }
            if hint_ok {
                for j in index..omega {
                    if y[j] != 0 {
                        hint_ok = false;
                        break;
                    }
                }
            }

            if tc.result != "valid" || !hint_ok {
                continue;
            }

            // --- Full verification path to reach UseHint ---
            let mut rho = [0u8; 32];
            rho.copy_from_slice(&pk[..32]);
            let mut t1 = vec![[0i32; N]; k];
            for i in 0..k {
                let off = 32 + i * 320;
                let mask: u64 = (1 << 10) - 1;
                let mut acc: u64 = 0;
                let mut ab: u32 = 0;
                let mut bi = 0;
                for jj in 0..N {
                    while ab < 10 {
                        acc |= (pk[off + bi] as u64) << ab;
                        bi += 1;
                        ab += 8;
                    }
                    t1[i][jj] = (acc & mask) as i32;
                    acc >>= 10;
                    ab -= 10;
                }
            }

            let a_hat = expand_a(&rho, k, l, stats);
            let c_tilde = &sig[..ct_len];
            let mut c = sample_in_ball(c_tilde, ps.tau);
            ntt(&mut c);

            let mut z_hat = z.clone();
            for p in z_hat.iter_mut() {
                for x in p.iter_mut() {
                    *x = x.rem_euclid(Q);
                }
                ntt(p);
            }
            let mut t1h = t1.clone();
            for p in t1h.iter_mut() {
                for x in p.iter_mut() {
                    *x = ((*x as i64 * (1i64 << D)) % Q as i64) as i32;
                }
                ntt(p);
            }

            for i in 0..k {
                let mut acc = [0i32; N];
                for j in 0..l {
                    for m in 0..N {
                        acc[m] = ((acc[m] as i64
                            + a_hat[i][j][m] as i64 * z_hat[j][m] as i64)
                            % Q as i64) as i32;
                    }
                }
                for m in 0..N {
                    let ct1 = (c[m] as i64 * t1h[i][m] as i64 % Q as i64) as i32;
                    acc[m] = (acc[m] - ct1).rem_euclid(Q);
                }
                inv_ntt(&mut acc);
                for m in 0..N {
                    let (_r1, r0, edge) = decompose_edge(acc[m], ps.gamma2);
                    if edge {
                        stats.decompose_edge += 1;
                        if h[i][m] == 1 {
                            stats.decompose_edge_h1 += 1;
                            // The pre-adjustment r0 was r0+1 (since the edge
                            // branch does r0-1). If the ORIGINAL r0 was 0,
                            // the returned value is -1, and we can detect it:
                            if r0 == -1 {
                                stats.decompose_edge_r0_0_h1 += 1;
                            }
                        }
                    }
                    if h[i][m] == 1 {
                        if r0 > 0 {
                            stats.use_hint_r0_pos_h1 += 1;
                        } else {
                            stats.use_hint_r0_nonpos_h1 += 1;
                        }
                        if r0 == 0 {
                            stats.use_hint_r0_zero_h1 += 1;
                        }
                    }
                }
            }
            // Silence unused-var lint for the message
            let _ = hex::decode(&tc.msg);
        }
    }
}

#[test]
fn show_path_coverage() {
    let mut stats = Stats::default();
    check_file(
        "testvectors/mldsa_44_verify_test.json",
        &ML_DSA_44,
        &mut stats,
    );
    check_file(
        "testvectors/mldsa_65_verify_test.json",
        &ML_DSA_65,
        &mut stats,
    );
    check_file(
        "testvectors/mldsa_87_verify_test.json",
        &ML_DSA_87,
        &mut stats,
    );
    eprintln!("\n=== Path coverage across all Wycheproof ML-DSA verify vectors ===");
    eprintln!(
        "Signatures with ||z||_inf >= gamma1 - beta: {}",
        stats.z_bound_exceeded
    );
    eprintln!(
        "Max ||z||_inf seen across all signatures:   {}",
        stats.max_z_inf_norm
    );
    eprintln!(
        "Hint cumulative count ever DECREASED:       {}",
        stats.hint_count_decreased
    );
    eprintln!(
        "RejNTTPoly saw value == q exactly:          {}",
        stats.rej_ntt_saw_q
    );
    eprintln!(
        "Decompose hit q-1 edge case:                {}",
        stats.decompose_edge
    );
    eprintln!(
        "  ...with a set hint bit at that position:  {}",
        stats.decompose_edge_h1
    );
    eprintln!(
        "  ...and input r == q-1 exactly (r0 was 0): {}",
        stats.decompose_edge_r0_0_h1
    );
    eprintln!(
        "UseHint with h=1 and r0 > 0:                {}",
        stats.use_hint_r0_pos_h1
    );
    eprintln!(
        "UseHint with h=1 and r0 <= 0:               {}",
        stats.use_hint_r0_nonpos_h1
    );
    eprintln!(
        "UseHint with h=1 and r0 == 0 exactly:       {}",
        stats.use_hint_r0_zero_h1
    );
}

const ZETAS: [i32; 256] = [
    0, 4808194, 3765607, 3761513, 5178923, 5496691, 5234739, 5178987, 7778734, 3542485, 2682288,
    2129892, 3764867, 7375178, 557458, 7159240, 5010068, 4317364, 2663378, 6705802, 4855975,
    7946292, 676590, 7044481, 5152541, 1714295, 2453983, 1460718, 7737789, 4795319, 2815639,
    2283733, 3602218, 3182878, 2740543, 4793971, 5269599, 2101410, 3704823, 1159875, 394148,
    928749, 1095468, 4874037, 2071829, 4361428, 3241972, 2156050, 3415069, 1759347, 7562881,
    4805951, 3756790, 6444618, 6663429, 4430364, 5483103, 3192354, 556856, 3870317, 2917338,
    1853806, 3345963, 1858416, 3073009, 1277625, 5744944, 3852015, 4183372, 5157610, 5258977,
    8106357, 2508980, 2028118, 1937570, 4564692, 2811291, 5396636, 7270901, 4158088, 1528066,
    482649, 1148858, 5418153, 7814814, 169688, 2462444, 5046034, 4213992, 4892034, 1987814,
    5183169, 1736313, 235407, 5130263, 3258457, 5801164, 1787943, 5989328, 6125690, 3482206,
    4197502, 7080401, 6018354, 7062739, 2461387, 3035980, 621164, 3901472, 7153756, 2925816,
    3374250, 1356448, 5604662, 2683270, 5601629, 4912752, 2312838, 7727142, 7921254, 348812,
    8052569, 1011223, 6026202, 4561790, 6458164, 6143691, 1744507, 1753, 6444997, 5720892,
    6924527, 2660408, 6600190, 8321269, 2772600, 1182243, 87208, 636927, 4415111, 4423672,
    6084020, 5095502, 4663471, 8352605, 822541, 1009365, 5926272, 6400920, 1596822, 4423473,
    4620952, 6695264, 4969849, 2678278, 4611469, 4829411, 635956, 8129971, 5925040, 4234153,
    6607829, 2192938, 6653329, 2387513, 4768667, 8111961, 5199961, 3747250, 2296099, 1239911,
    4541938, 3195676, 2642980, 1254190, 8368000, 2998219, 141835, 8291116, 2513018, 7025525,
    613238, 7070156, 6161950, 7921677, 6458423, 4040196, 4908348, 2039144, 6500539, 7561656,
    6201452, 6757063, 2105286, 6006015, 6346610, 586241, 7200804, 527981, 5637006, 6903432,
    1994046, 2491325, 6987258, 507927, 7192532, 7655613, 6545891, 5346675, 8041997, 2647994,
    3009748, 5767564, 4148469, 749577, 4357667, 3980599, 2569011, 6764887, 1723229, 1665318,
    2028038, 1163598, 5011144, 3994671, 8368538, 7009900, 3020393, 3363542, 214880, 545376,
    7609976, 3105558, 7277073, 508145, 7826699, 860144, 3430436, 140244, 6866265, 6195333,
    3123762, 2358373, 6187330, 5365997, 6663603, 2926054, 7987710, 8077412, 3531229, 4405932,
    4606686, 1900052, 7598542, 1054478, 7648983,
];
