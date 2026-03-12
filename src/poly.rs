//! Polynomial arithmetic in Rq = Zq[X]/(X^256 + 1) and the associated NTT.
//! This follows FIPS 204, Section 7.5 and Algorithms 41/42.

use crate::params::{N, Q};

/// A polynomial in Rq with coefficients in [0, q).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly {
    pub c: [i32; N],
}

impl Poly {
    pub fn zero() -> Self {
        Poly { c: [0; N] }
    }

    /// Reduce all coefficients into the range [0, Q).
    pub fn reduce(&mut self) {
        for x in self.c.iter_mut() {
            *x = x.rem_euclid(Q);
        }
    }

    /// Pointwise multiply two polynomials already in NTT form.
    pub fn ntt_mul(&self, other: &Poly) -> Poly {
        let mut out = Poly::zero();
        for i in 0..N {
            out.c[i] = mul_mod(self.c[i], other.c[i]);
        }
        out
    }

    /// Add two polynomials, producing coefficients in [0, Q).
    pub fn add(&self, other: &Poly) -> Poly {
        let mut out = Poly::zero();
        for i in 0..N {
            out.c[i] = add_mod(self.c[i], other.c[i]);
        }
        out
    }

    /// Subtract other from self, producing coefficients in [0, Q).
    pub fn sub(&self, other: &Poly) -> Poly {
        let mut out = Poly::zero();
        for i in 0..N {
            out.c[i] = sub_mod(self.c[i], other.c[i]);
        }
        out
    }

    /// Multiply every coefficient by 2^d and reduce mod Q.
    pub fn shift_left(&self, d: u32) -> Poly {
        let mut out = Poly::zero();
        let factor = 1i64 << d;
        for i in 0..N {
            out.c[i] = ((self.c[i] as i64 * factor) % (Q as i64)) as i32;
        }
        out
    }

    /// Infinity norm: max of |c mod± q| over all coefficients.
    pub fn inf_norm(&self) -> i32 {
        let mut m = 0i32;
        for &x in &self.c {
            let v = centered(x).abs();
            if v > m {
                m = v;
            }
        }
        m
    }
}

/// A vector of polynomials.
#[derive(Clone, Debug)]
pub struct PolyVec {
    pub v: Vec<Poly>,
}

impl PolyVec {
    pub fn zero(len: usize) -> Self {
        PolyVec {
            v: (0..len).map(|_| Poly::zero()).collect(),
        }
    }

    pub fn ntt(&mut self) {
        for p in self.v.iter_mut() {
            ntt(p);
        }
    }

    pub fn inv_ntt(&mut self) {
        for p in self.v.iter_mut() {
            inv_ntt(p);
        }
    }

    pub fn sub(&self, other: &PolyVec) -> PolyVec {
        PolyVec {
            v: self.v.iter().zip(&other.v).map(|(a, b)| a.sub(b)).collect(),
        }
    }

    /// Infinity norm over all coefficients in all polynomials.
    pub fn inf_norm(&self) -> i32 {
        self.v.iter().map(Poly::inf_norm).max().unwrap_or(0)
    }
}

/// The public matrix A_hat with k rows and l columns, each entry already in NTT form.
pub struct PolyMat {
    pub rows: Vec<Vec<Poly>>,
}

impl PolyMat {
    /// Compute A_hat * z_hat (matrix-vector product in the NTT domain).
    /// All inputs must already be in NTT form; the output is also in NTT form.
    pub fn mul_vec(&self, z: &PolyVec) -> PolyVec {
        let k = self.rows.len();
        let l = z.v.len();
        let mut out = PolyVec::zero(k);
        for i in 0..k {
            let mut acc = Poly::zero();
            for j in 0..l {
                let prod = self.rows[i][j].ntt_mul(&z.v[j]);
                acc = acc.add(&prod);
            }
            out.v[i] = acc;
        }
        out
    }
}

/// Multiply two values in [0, Q) and reduce mod Q.
#[inline]
fn mul_mod(a: i32, b: i32) -> i32 {
    ((a as i64 * b as i64) % Q as i64) as i32
}

/// Add two values in [0, Q) and reduce mod Q.
#[inline]
fn add_mod(a: i32, b: i32) -> i32 {
    let s = a + b;
    if s >= Q {
        s - Q
    } else {
        s
    }
}

/// Subtract b from a, both in [0, Q), producing a value in [0, Q).
#[inline]
fn sub_mod(a: i32, b: i32) -> i32 {
    let d = a - b;
    if d < 0 {
        d + Q
    } else {
        d
    }
}

/// Return the centered representative of x mod q, i.e. a value in (-q/2, q/2].
#[inline]
fn centered(x: i32) -> i32 {
    let x = x.rem_euclid(Q);
    if x > Q / 2 {
        x - Q
    } else {
        x
    }
}

/// Precomputed powers of zeta for the NTT: ZETAS[k] = zeta^brv(k) mod q,
/// where zeta = 1753 is a primitive 512-th root of unity mod q and brv
/// is bit-reversal of an 8-bit integer. See FIPS 204, Appendix B.
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

/// The modular inverse of 256 mod q, used to normalise at the end of the inverse NTT.
const N_INV: i32 = 8_347_681;

/// Forward NTT (FIPS 204, Algorithm 41). Coefficients are mapped into the
/// NTT domain in bit-reversed order. Input coefficients must be in [0, Q).
pub fn ntt(p: &mut Poly) {
    let mut k = 0usize;
    let mut len = 128usize;
    while len >= 1 {
        let mut start = 0usize;
        while start < N {
            k += 1;
            let zeta = ZETAS[k];
            for j in start..start + len {
                let t = mul_mod(zeta, p.c[j + len]);
                p.c[j + len] = sub_mod(p.c[j], t);
                p.c[j] = add_mod(p.c[j], t);
            }
            start += 2 * len;
        }
        len /= 2;
    }
}

/// Inverse NTT (FIPS 204, Algorithm 42). Produces coefficients in [0, Q).
pub fn inv_ntt(p: &mut Poly) {
    let mut k = 256usize;
    let mut len = 1usize;
    while len < N {
        let mut start = 0usize;
        while start < N {
            k -= 1;
            // negate the zeta: in Z_q, -z = q - z (for z != 0).
            let zeta = Q - ZETAS[k];
            for j in start..start + len {
                let t = p.c[j];
                p.c[j] = add_mod(t, p.c[j + len]);
                p.c[j + len] = mul_mod(zeta, sub_mod(t, p.c[j + len]));
            }
            start += 2 * len;
        }
        len *= 2;
    }
    // Scale by 1/256.
    for x in p.c.iter_mut() {
        *x = mul_mod(*x, N_INV);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ntt_roundtrip() {
        let mut p = Poly::zero();
        for i in 0..N {
            p.c[i] = (i as i32 * 31 + 7) % Q;
        }
        let orig = p.clone();
        ntt(&mut p);
        inv_ntt(&mut p);
        assert_eq!(p, orig);
    }

    #[test]
    fn ntt_mul_matches_schoolbook() {
        // Multiply two small polynomials and compare against schoolbook
        // multiplication in Rq.
        let mut a = Poly::zero();
        let mut b = Poly::zero();
        for i in 0..N {
            a.c[i] = (i as i32 * 3 + 1) % Q;
            b.c[i] = (i as i32 * 7 + 5) % Q;
        }
        // Schoolbook: c = a * b mod (X^256 + 1)
        let mut school = [0i64; N];
        for i in 0..N {
            for j in 0..N {
                let prod = (a.c[i] as i64) * (b.c[j] as i64) % (Q as i64);
                if i + j < N {
                    school[i + j] = (school[i + j] + prod) % (Q as i64);
                } else {
                    school[i + j - N] = (school[i + j - N] - prod).rem_euclid(Q as i64);
                }
            }
        }
        // Via NTT
        ntt(&mut a);
        ntt(&mut b);
        let mut c = a.ntt_mul(&b);
        inv_ntt(&mut c);
        for i in 0..N {
            assert_eq!(c.c[i] as i64, school[i], "mismatch at coeff {}", i);
        }
    }
}
