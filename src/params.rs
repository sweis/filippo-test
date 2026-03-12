//! ML-DSA parameter sets from FIPS 204, Table 1.

/// The modulus q = 2^23 - 2^13 + 1 = 8380417.
pub const Q: i32 = 8_380_417;

/// Degree of the polynomial modulus X^n + 1.
pub const N: usize = 256;

/// Number of dropped bits in the public key compression.
pub const D: u32 = 13;

/// Parameter set for a security level of ML-DSA.
#[derive(Clone, Copy, Debug)]
pub struct ParamSet {
    /// Hamming weight of the challenge polynomial c.
    pub tau: usize,
    /// Collision strength in bits. c_tilde is lambda/4 bytes.
    pub lambda: usize,
    /// Coefficient range parameter for the signature vector z: z in [-gamma1+1, gamma1].
    pub gamma1: i32,
    /// Low-order rounding range: coefficients are reduced mod 2*gamma2.
    pub gamma2: i32,
    /// Number of rows of the public matrix A.
    pub k: usize,
    /// Number of columns of the public matrix A.
    pub l: usize,
    /// Infinity-norm bound on the secret key coefficients (not used in verification,
    /// but beta = tau * eta is the rejection bound).
    pub eta: i32,
    /// Rejection bound beta = tau * eta.
    pub beta: i32,
    /// Maximum number of 1 bits in the hint vector h.
    pub omega: usize,
}

/// ML-DSA-44 parameters (FIPS 204, Table 1).
pub const ML_DSA_44: ParamSet = ParamSet {
    tau: 39,
    lambda: 128,
    gamma1: 1 << 17,
    gamma2: (Q - 1) / 88, // = 95232
    k: 4,
    l: 4,
    eta: 2,
    beta: 39 * 2,
    omega: 80,
};

/// ML-DSA-65 parameters (FIPS 204, Table 1).
pub const ML_DSA_65: ParamSet = ParamSet {
    tau: 49,
    lambda: 192,
    gamma1: 1 << 19,
    gamma2: (Q - 1) / 32, // = 261888
    k: 6,
    l: 5,
    eta: 4,
    beta: 49 * 4,
    omega: 55,
};

/// ML-DSA-87 parameters (FIPS 204, Table 1).
pub const ML_DSA_87: ParamSet = ParamSet {
    tau: 60,
    lambda: 256,
    gamma1: 1 << 19,
    gamma2: (Q - 1) / 32, // = 261888
    k: 8,
    l: 7,
    eta: 2,
    beta: 60 * 2,
    omega: 75,
};

impl ParamSet {
    /// Length in bytes of the commitment hash c_tilde.
    pub fn c_tilde_bytes(&self) -> usize {
        self.lambda / 4
    }

    /// Number of bits used to encode a coefficient of z. Per FIPS 204, z is
    /// encoded with BitPack(z, gamma1 - 1, gamma1), so each coefficient uses
    /// bitlen(gamma1 - 1 + gamma1) = bitlen(2*gamma1 - 1) bits. Since gamma1
    /// is a power of two, this is 1 + log2(gamma1).
    pub fn z_bits_per_coeff(&self) -> u32 {
        1 + (self.gamma1 as u32).ilog2()
    }

    /// Bytes per polynomial in the z encoding.
    pub fn z_poly_bytes(&self) -> usize {
        (self.z_bits_per_coeff() as usize) * N / 8
    }

    /// Bits used to encode a coefficient of w1 = HighBits. This is
    /// bitlen((q-1)/(2*gamma2) - 1).
    pub fn w1_bits_per_coeff(&self) -> u32 {
        let m = ((Q - 1) / (2 * self.gamma2)) as u32;
        // m is the number of representatives; range is [0, m-1], so we need
        // bitlen(m-1) bits.
        32 - (m - 1).leading_zeros()
    }

    /// Bytes per polynomial in the w1 encoding.
    pub fn w1_poly_bytes(&self) -> usize {
        (self.w1_bits_per_coeff() as usize) * N / 8
    }

    /// Total public key size in bytes: 32 (rho) + k * 320 (t1 packed with 10 bits/coeff).
    pub fn pk_bytes(&self) -> usize {
        32 + self.k * (10 * N / 8)
    }

    /// Total signature size in bytes: c_tilde + l*z_poly + omega + k.
    pub fn sig_bytes(&self) -> usize {
        self.c_tilde_bytes() + self.l * self.z_poly_bytes() + self.omega + self.k
    }
}
