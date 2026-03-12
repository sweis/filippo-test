//! Bit packing and unpacking routines from FIPS 204, Section 7.1.
//!
//! These translate polynomials to byte strings and back. Coefficients are
//! packed little-endian: the LSB of each coefficient is written first.

use crate::params::N;
use crate::poly::Poly;

/// Pack a polynomial whose coefficients are in [0, 2^bits). Returns
/// exactly bits * 256 / 8 bytes. See FIPS 204, Algorithm 16 (SimpleBitPack).
pub fn simple_bit_pack(p: &Poly, bits: u32) -> Vec<u8> {
    let mut out = Vec::with_capacity((bits as usize) * N / 8);
    let mut acc: u64 = 0;
    let mut acc_bits: u32 = 0;
    for &coeff in &p.c {
        acc |= (coeff as u64) << acc_bits;
        acc_bits += bits;
        while acc_bits >= 8 {
            out.push(acc as u8);
            acc >>= 8;
            acc_bits -= 8;
        }
    }
    debug_assert_eq!(acc_bits, 0);
    out
}

/// Inverse of simple_bit_pack. Reads coefficients in [0, 2^bits) from bytes.
/// See FIPS 204, Algorithm 18 (SimpleBitUnpack). Does not validate that the
/// resulting coefficients are within any particular bound.
pub fn simple_bit_unpack(bytes: &[u8], bits: u32) -> Poly {
    let mut p = Poly::zero();
    let mask: u64 = (1 << bits) - 1;
    let mut acc: u64 = 0;
    let mut acc_bits: u32 = 0;
    let mut byte_idx = 0usize;
    for i in 0..N {
        while acc_bits < bits {
            acc |= (bytes[byte_idx] as u64) << acc_bits;
            byte_idx += 1;
            acc_bits += 8;
        }
        p.c[i] = (acc & mask) as i32;
        acc >>= bits;
        acc_bits -= bits;
    }
    p
}

/// Inverse of BitPack (FIPS 204, Algorithm 19). The encoding maps each
/// coefficient w in [-a, b] to the integer (b - w), which is in [0, a+b],
/// and packs that with bitlen(a+b) bits. So to unpack, we read an unsigned
/// value and subtract it from b.
pub fn bit_unpack(bytes: &[u8], a: i32, b: i32) -> Poly {
    let bits = 32 - ((a + b) as u32).leading_zeros();
    let raw = simple_bit_unpack(bytes, bits);
    let mut out = Poly::zero();
    for i in 0..N {
        out.c[i] = b - raw.c[i];
    }
    out
}

/// Decode the hint vector h from its packed byte form (FIPS 204, Algorithm 21,
/// HintBitUnpack). Returns None if the encoding is invalid.
///
/// Encoding layout: omega bytes of hint positions, then k bytes of per-row
/// cumulative counts. Within each row the positions must be strictly increasing,
/// unused position slots must be zero, and the final count must not exceed omega.
pub fn hint_bit_unpack(y: &[u8], k: usize, omega: usize) -> Option<Vec<Poly>> {
    debug_assert_eq!(y.len(), omega + k);
    let mut h = vec![Poly::zero(); k];
    let mut index: usize = 0;
    for i in 0..k {
        // y[omega + i] is the cumulative count after row i.
        let end = y[omega + i] as usize;
        if end < index || end > omega {
            return None;
        }
        let mut first = true;
        for j in index..end {
            let pos = y[j] as usize;
            if !first && pos <= (y[j - 1] as usize) {
                // Positions must be strictly increasing within a row.
                return None;
            }
            first = false;
            h[i].c[pos] = 1;
        }
        index = end;
    }
    // The remaining position slots up to omega must be zero.
    for j in index..omega {
        if y[j] != 0 {
            return None;
        }
    }
    Some(h)
}
