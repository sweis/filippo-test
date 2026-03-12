# FIPS 204 Compliance Cross-Check

The implementation was audited line-by-line against `NIST.FIPS.204.pdf`
(August 2024 final). Every algorithm on the verification path was compared.
No discrepancies were found.

## Summary

| Spec reference         | Implementation site            | Result |
|------------------------|--------------------------------|--------|
| Table 1 (parameters)   | `src/params.rs`                | ✓ exact |
| Table 2 (sizes)        | `params.rs::pk_bytes/sig_bytes`| ✓ exact (1312/2420, 1952/3309, 2592/4627) |
| §3.6.2 (length checks) | `lib.rs::pk_decode/sig_decode` | ✓ explicit early return on length mismatch |
| Algorithm 3            | `lib.rs::verify`               | ✓ domain-sep `0‖\|ctx\|‖ctx‖M`, ctx ≤ 255 |
| Algorithm 8            | `lib.rs::verify_internal`      | ✓ see note ① |
| Algorithm 14           | `sample.rs::rej_ntt_poly`      | ✓ `b0 \| b1<<8 \| (b2&0x7F)<<16`, `< q` |
| Algorithm 16/17        | `encode.rs::simple_bit_pack`   | ✓ little-endian, `bitlen b` bits/coeff |
| Algorithm 18/19        | `encode.rs::*_bit_unpack`      | ✓ inverse; spec confirms range always OK for our calls ② |
| Algorithm 21           | `encode.rs::hint_bit_unpack`   | ✓ all three ⊥ conditions (line 4, 9, 17) |
| Algorithm 23           | `lib.rs::pk_decode`            | ✓ 10 bits/t1-coeff; spec note: "always in correct range" |
| Algorithm 27           | `lib.rs::sig_decode`           | ✓ z via `BitUnpack(_, γ₁−1, γ₁)`; spec note: γ₁ power of 2 ⇒ range OK |
| Algorithm 28           | `lib.rs::w1_encode`            | ✓ `SimpleBitPack(_, (q−1)/(2γ₂)−1)` |
| Algorithm 29           | `sample.rs::sample_in_ball`    | ✓ 8 sign bytes LE, Fisher–Yates, `h[i+τ−256]` bit order |
| Algorithm 30           | `sample.rs::rej_ntt_poly`      | ✓ `G` = SHAKE128, 3-byte rejection |
| Algorithm 32           | `sample.rs::expand_a`          | ✓ seed `ρ‖s‖r` (col then row byte) |
| Algorithm 36           | `hint.rs::decompose`           | ✓ `mod± 2γ₂`, q−1 edge case with `r₀ ← r₀−1` |
| Algorithm 40           | `hint.rs::use_hint`            | ✓ `r₀ > 0` (strict), `mod m` |
| Algorithm 41           | `poly.rs::ntt`                 | ✓ `m←0`, `len←128→1`, butterfly matches |
| Algorithm 42           | `poly.rs::inv_ntt`             | ✓ `m←256`, `len←1→128`, `f=8347681` |
| Appendix B (zetas)     | `poly.rs::ZETAS`               | ✓ all 256 entries identical |
| H / G definitions      | `sample.rs`                    | ✓ H=SHAKE256, G=SHAKE128 |
| `mod±` definition      | `poly.rs::centered`, `hint.rs` | ✓ range `(−⌈α/2⌉, ⌊α/2⌋]` |
| `‖·‖∞` definition      | `poly.rs::inf_norm`            | ✓ `max \|centered(w_i)\|` |

## Notes

**① z-bound check ordering.** FIPS 204 Algorithm 8 evaluates the z-bound at
step 13 as the final conjunction `[[‖z‖∞ < γ₁−β]] ∧ [[c̃ = c̃′]]`. This
implementation checks the z-bound immediately after `sigDecode`, before any
hashing or NTT work. The two orderings are semantically equivalent (both
conditions must hold for acceptance). The early check is a conventional
short-circuit; it does not affect correctness or the mutation-test findings,
since no Wycheproof vector fails the bound either way.

**② Spec-confirmed range safety.** The PDF explicitly annotates that
`SimpleBitUnpack` in `pkDecode` (Alg 23 line 3) and `BitUnpack` in `sigDecode`
(Alg 27 line 3) *always* return in-range coefficients for the parameter
choices here. The implementation therefore adds no redundant range checks,
which would distort mutation results.

## Conclusion

The implementation is a faithful transcription of the FIPS 204 verification
path. No behavioural differences from the specification were identified. The
mutation-test findings in `REPORT.md` therefore reflect genuine gaps in the
Wycheproof corpus, not artefacts of this implementation.
