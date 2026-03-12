//! Rounding primitives from FIPS 204, Section 7.4.

use crate::params::Q;

/// Split r into a high part r1 and a low centered part r0, so that
/// r == r1 * (2*gamma2) + r0 mod q, with r0 in (-gamma2, gamma2].
/// See FIPS 204, Algorithm 36 (Decompose).
pub fn decompose(r: i32, gamma2: i32) -> (i32, i32) {
    let r_plus = r.rem_euclid(Q);
    // mod± (2*gamma2): the unique r0 in (-gamma2, gamma2] with r0 ≡ r_plus.
    let two_g2 = 2 * gamma2;
    let mut r0 = r_plus % two_g2;
    if r0 > gamma2 {
        r0 -= two_g2;
    }
    if r_plus - r0 == Q - 1 {
        // Edge case: when the high part would wrap to the top value, fold it
        // down to zero and shift r0 by one.
        (0, r0 - 1)
    } else {
        ((r_plus - r0) / two_g2, r0)
    }
}

/// Recover the high bits of the verifier's approximation using a 1-bit hint.
/// See FIPS 204, Algorithm 40 (UseHint).
pub fn use_hint(h: i32, r: i32, gamma2: i32) -> i32 {
    let m = (Q - 1) / (2 * gamma2);
    let (r1, r0) = decompose(r, gamma2);
    if h == 1 {
        if r0 > 0 {
            (r1 + 1).rem_euclid(m)
        } else {
            (r1 - 1).rem_euclid(m)
        }
    } else {
        r1
    }
}
