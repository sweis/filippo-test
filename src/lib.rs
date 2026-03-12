//! A clean-room ML-DSA (FIPS 204) signature verifier.
//!
//! Only verification is implemented; there is no key generation or signing.
//! The implementation follows the algorithms in FIPS 204 as directly as
//! possible so that each spec-level check maps to a single code site.

mod encode;
mod hint;
pub mod params;
mod poly;
mod sample;

use encode::{bit_unpack, hint_bit_unpack, simple_bit_pack, simple_bit_unpack};
use hint::use_hint;
use params::{ParamSet, D, N, Q};
use poly::{ntt, Poly, PolyVec};
use sample::{expand_a, sample_in_ball, shake256};

pub use params::{ML_DSA_44, ML_DSA_65, ML_DSA_87};

/// Reasons a verification may fail. Useful for testing to distinguish early
/// parse-time rejections from cryptographic mismatches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifyError {
    ContextTooLong,
    BadPublicKeyLength,
    BadSignatureLength,
    MalformedHint,
    ZOutOfRange,
    ChallengeMismatch,
}

/// FIPS 204, Algorithm 22 (pkDecode).
/// Returns (rho, t1) where t1 is a vector of k polynomials with
/// coefficients in [0, 2^10).
fn pk_decode(ps: &ParamSet, pk: &[u8]) -> Option<([u8; 32], PolyVec)> {
    if pk.len() != ps.pk_bytes() {
        return None;
    }
    let mut rho = [0u8; 32];
    rho.copy_from_slice(&pk[..32]);
    // Each t1[i] is packed with 10 bits per coefficient: bitlen(q-1) - d = 23 - 13.
    let t1_bits: u32 = 10;
    let t1_poly_bytes = (t1_bits as usize) * N / 8; // = 320
    let mut t1 = PolyVec::zero(ps.k);
    let mut off = 32;
    for i in 0..ps.k {
        t1.v[i] = simple_bit_unpack(&pk[off..off + t1_poly_bytes], t1_bits);
        off += t1_poly_bytes;
    }
    Some((rho, t1))
}

/// FIPS 204, Algorithm 26 (sigDecode).
/// Returns (c_tilde, z, h). The hint vector is None if the encoding is invalid.
fn sig_decode<'a>(
    ps: &ParamSet,
    sig: &'a [u8],
) -> Option<(&'a [u8], PolyVec, Option<Vec<Poly>>)> {
    if sig.len() != ps.sig_bytes() {
        return None;
    }
    let ct_len = ps.c_tilde_bytes();
    let c_tilde = &sig[..ct_len];

    let z_poly_bytes = ps.z_poly_bytes();
    let mut z = PolyVec::zero(ps.l);
    let mut off = ct_len;
    // z is encoded with BitPack(w, gamma1 - 1, gamma1): range [-(gamma1-1), gamma1].
    for i in 0..ps.l {
        z.v[i] = bit_unpack(&sig[off..off + z_poly_bytes], ps.gamma1 - 1, ps.gamma1);
        off += z_poly_bytes;
    }

    let h_bytes = &sig[off..off + ps.omega + ps.k];
    let h = hint_bit_unpack(h_bytes, ps.k, ps.omega);

    Some((c_tilde, z, h))
}

/// FIPS 204, Algorithm 28 (w1Encode). Pack each polynomial of w1 with
/// bitlen((q-1)/(2*gamma2) - 1) bits per coefficient.
fn w1_encode(ps: &ParamSet, w1: &PolyVec) -> Vec<u8> {
    let bits = ps.w1_bits_per_coeff();
    let mut out = Vec::with_capacity(ps.k * ps.w1_poly_bytes());
    for p in &w1.v {
        out.extend_from_slice(&simple_bit_pack(p, bits));
    }
    out
}

/// FIPS 204, Algorithm 8 (ML-DSA.Verify_internal).
/// Verifies a signature against the already-formatted message M'.
pub fn verify_internal(
    ps: &ParamSet,
    pk: &[u8],
    m_prime: &[u8],
    sig: &[u8],
) -> Result<(), VerifyError> {
    // 1. (rho, t1) <- pkDecode(pk)
    let (rho, t1) = pk_decode(ps, pk).ok_or(VerifyError::BadPublicKeyLength)?;

    // 2. (c_tilde, z, h) <- sigDecode(sig); if h = bot, return false
    let (c_tilde, z, h) = sig_decode(ps, sig).ok_or(VerifyError::BadSignatureLength)?;
    let h = h.ok_or(VerifyError::MalformedHint)?;

    // 3. if ||z||_inf >= gamma1 - beta, return false
    if z.inf_norm() >= ps.gamma1 - ps.beta {
        return Err(VerifyError::ZOutOfRange);
    }

    // 4. A_hat <- ExpandA(rho)
    let a_hat = expand_a(&rho, ps.k, ps.l);

    // 5. tr <- H(pk, 64)
    let tr = shake256(&[pk], 64);

    // 6. mu <- H(tr || M', 64)
    let mu = shake256(&[&tr, m_prime], 64);

    // 7. c <- SampleInBall(c_tilde)
    let c = sample_in_ball(c_tilde, ps.tau);

    // 8. w'_approx <- invNTT(A_hat * NTT(z) - NTT(c) * NTT(t1 * 2^d))
    // Note z is currently signed in (-gamma1, gamma1]; reduce to [0, q) first.
    let mut z_hat = z.clone();
    for p in z_hat.v.iter_mut() {
        p.reduce();
    }
    z_hat.ntt();

    let mut c_hat = c.clone();
    ntt(&mut c_hat);

    // t1 * 2^d
    let mut t1_shifted = PolyVec::zero(ps.k);
    for i in 0..ps.k {
        t1_shifted.v[i] = t1.v[i].shift_left(D);
    }
    t1_shifted.ntt();

    let az = a_hat.mul_vec(&z_hat);

    // c * (t1 * 2^d)
    let mut ct1 = PolyVec::zero(ps.k);
    for i in 0..ps.k {
        ct1.v[i] = c_hat.ntt_mul(&t1_shifted.v[i]);
    }

    let mut w_approx = az.sub(&ct1);
    w_approx.inv_ntt();

    // 9. w1' <- UseHint(h, w'_approx)
    let mut w1 = PolyVec::zero(ps.k);
    for i in 0..ps.k {
        for j in 0..N {
            w1.v[i].c[j] = use_hint(h[i].c[j], w_approx.v[i].c[j], ps.gamma2);
        }
    }

    // 10. c_tilde' <- H(mu || w1Encode(w1'), lambda/4)
    let w1_packed = w1_encode(ps, &w1);
    let c_tilde_prime = shake256(&[&mu, &w1_packed], ps.c_tilde_bytes());

    // 11. return [[c_tilde == c_tilde']]
    // (The ||z|| bound is already checked above; the hint count bound is
    // enforced by hint_bit_unpack refusing > omega entries.)
    if c_tilde == c_tilde_prime.as_slice() {
        Ok(())
    } else {
        Err(VerifyError::ChallengeMismatch)
    }
}

/// FIPS 204, Algorithm 3 (ML-DSA.Verify).
/// ctx is the optional context string; the empty slice is the default.
pub fn verify(
    ps: &ParamSet,
    pk: &[u8],
    msg: &[u8],
    sig: &[u8],
    ctx: &[u8],
) -> Result<(), VerifyError> {
    // 1. if |ctx| > 255, return false
    if ctx.len() > 255 {
        return Err(VerifyError::ContextTooLong);
    }
    // 2. M' <- BytesToBits(IntegerToBytes(0, 1) || IntegerToBytes(|ctx|, 1) || ctx || M)
    let mut m_prime = Vec::with_capacity(2 + ctx.len() + msg.len());
    m_prime.push(0);
    m_prime.push(ctx.len() as u8);
    m_prime.extend_from_slice(ctx);
    m_prime.extend_from_slice(msg);
    // 3. return ML-DSA.Verify_internal(pk, M', sig)
    verify_internal(ps, pk, &m_prime, sig)
}

// Reject the no-op t1 range check that FIPS 204 §7.2 mentions for pkDecode:
// SimpleBitUnpack with 10 bits can produce values up to 2^10 - 1 = 1023, but
// the valid range is [0, 2^(bitlen(q-1)-d) - 1] = [0, 1023], so every decoded
// value is in range by construction. We intentionally do not add a redundant
// check for it (it would be dead code and skew mutation results).
#[allow(dead_code)]
const _: () = {
    // Compile-time sanity: 2^10 - 1 is the max possible and is within range.
    assert!(((1 << 10) - 1) < (1 << ((32 - (Q as u32 - 1).leading_zeros()) - D)));
    // Actually the bound is 2^(bitlen(q-1) - d) = 2^(23 - 13) = 2^10 = 1024,
    // and max decoded is 1023 < 1024. So the assertion holds.
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn param_sizes() {
        // FIPS 204, Table 2: public key and signature sizes.
        assert_eq!(ML_DSA_44.pk_bytes(), 1312);
        assert_eq!(ML_DSA_44.sig_bytes(), 2420);
        assert_eq!(ML_DSA_65.pk_bytes(), 1952);
        assert_eq!(ML_DSA_65.sig_bytes(), 3309);
        assert_eq!(ML_DSA_87.pk_bytes(), 2592);
        assert_eq!(ML_DSA_87.sig_bytes(), 4627);
    }
}
