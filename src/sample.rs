//! Rejection-sampling primitives from FIPS 204, Section 7.3.

use crate::params::{N, Q};
use crate::poly::{Poly, PolyMat};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128, Shake256,
};

/// Compute SHAKE256(input, out_len). Convenience wrapper for one-shot hashing.
pub fn shake256(inputs: &[&[u8]], out_len: usize) -> Vec<u8> {
    let mut h = Shake256::default();
    for part in inputs {
        h.update(part);
    }
    let mut out = vec![0u8; out_len];
    h.finalize_xof().read(&mut out);
    out
}

/// Generate one polynomial of A_hat by rejection sampling uniform values in
/// [0, q) from SHAKE128 output. See FIPS 204, Algorithm 30 (RejNTTPoly).
fn rej_ntt_poly(seed: &[u8; 32], col: u8, row: u8) -> Poly {
    let mut h = Shake128::default();
    h.update(seed);
    h.update(&[col, row]);
    let mut reader = h.finalize_xof();

    let mut p = Poly::zero();
    let mut buf = [0u8; 3];
    let mut i = 0usize;
    while i < N {
        reader.read(&mut buf);
        // Three bytes, mask the top bit to get a 23-bit value.
        let v = (buf[0] as u32) | ((buf[1] as u32) << 8) | (((buf[2] & 0x7F) as u32) << 16);
        if (v as i32) < Q {
            p.c[i] = v as i32;
            i += 1;
        }
    }
    p
}

/// Expand the public matrix A_hat (k rows, l cols) from seed rho. Each entry
/// is generated in NTT form directly (the rejection sampling produces values
/// already suited to NTT-domain use). See FIPS 204, Algorithm 32 (ExpandA).
pub fn expand_a(rho: &[u8; 32], k: usize, l: usize) -> PolyMat {
    let mut rows = Vec::with_capacity(k);
    for r in 0..k {
        let mut row = Vec::with_capacity(l);
        for c in 0..l {
            row.push(rej_ntt_poly(rho, c as u8, r as u8));
        }
        rows.push(row);
    }
    PolyMat { rows }
}

/// Sample the challenge polynomial c with exactly tau coefficients of ±1 and
/// the rest zero, from the commitment hash c_tilde. See FIPS 204, Algorithm 29
/// (SampleInBall).
pub fn sample_in_ball(c_tilde: &[u8], tau: usize) -> Poly {
    let mut h = Shake256::default();
    h.update(c_tilde);
    let mut reader = h.finalize_xof();

    // First 8 bytes give the sign bits (64 bits, enough for tau <= 64).
    let mut sign_bytes = [0u8; 8];
    reader.read(&mut sign_bytes);
    let mut signs: u64 = u64::from_le_bytes(sign_bytes);

    let mut c = Poly::zero();
    for i in (N - tau)..N {
        // Rejection-sample j in [0, i].
        let mut jbuf = [0u8];
        loop {
            reader.read(&mut jbuf);
            if (jbuf[0] as usize) <= i {
                break;
            }
        }
        let j = jbuf[0] as usize;
        c.c[i] = c.c[j];
        // The low bit of signs determines the sign: 0 -> +1, 1 -> -1.
        // Store in [0, Q) form: +1 as 1, -1 as Q-1.
        c.c[j] = if signs & 1 == 0 { 1 } else { Q - 1 };
        signs >>= 1;
    }
    c
}
