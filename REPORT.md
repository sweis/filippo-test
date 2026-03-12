# Wycheproof ML-DSA Verify Coverage: Mutation Test Report

## Method

A clean-room ML-DSA signature verifier was implemented in Rust directly from
the FIPS 204 specification, without reference to any existing ML-DSA
implementation or exposition. All 239 Wycheproof `mldsa_*_verify_test.json`
test cases pass (80 ML-DSA-44, 83 ML-DSA-65, 76 ML-DSA-87). `cargo-mutants`
was then run against this implementation with only the Wycheproof harness as
the oracle.

## Summary

| Outcome     | Count | Notes                                          |
|-------------|-------|------------------------------------------------|
| Caught      | 223   | Wycheproof kills the mutant.                   |
| **Missed**  | **28**| Mutant survives; see analysis below.           |
| Unviable    | 25    | Fails to compile (no `Default` on `Poly` etc). |
| Timeout     | 4     | Infinite rejection-sampling loops.             |
| **Total**   | 280   |                                                |

Of the 28 survivors, **20 indicate genuine gaps** in the Wycheproof corpus and
**8 are equivalent mutants** (no observable behaviour change).

---

## Genuine gaps

### G1. The ‖z‖∞ bound check is never exercised — 15 survivors

FIPS 204 Algorithm 8, step 3 mandates rejecting any signature where
`‖z‖∞ ≥ γ₁ − β`. Across every structurally-valid signature in the Wycheproof
set, the maximum observed `‖z‖∞` is **524 167** — exactly one below the
ML-DSA-87 bound of `2¹⁹ − 120 = 524 168`. **No signature triggers the check.**

| Surviving mutant                                        | Effect                                   |
|---------------------------------------------------------|------------------------------------------|
| `lib.rs:108` `γ₁ − β` → `γ₁ + β`                        | Bound doubled; never fires anyway.       |
| `poly.rs:63` `inf_norm` → `0`, `1`, or `−1` (×3)        | Norm computation deleted.                |
| `poly.rs:66` `>` → `==`, `<`, `>=` (×3)                 | Max tracker broken.                      |
| `poly.rs:107` `PolyVec::inf_norm` → `0`, `1`, `−1` (×3) | Norm computation deleted.                |
| `poly.rs:166` `centered` → `0`, `1`, `−1` (×3)          | Centered-remainder helper deleted.       |
| `poly.rs:167` `>` → `>=`                                | Off-by-one when x = ⌊q/2⌋.               |
| `poly.rs:168` `x − q` → `x / q`                         | Returns 0 instead of negative.           |

**Impact:** A verifier that **skips the z-bound check entirely** passes every
Wycheproof ML-DSA test. This is the most security-relevant gap — the z-bound
is what ties the signature to the short-vector hardness assumption.

**Suggested vector:** An invalid signature whose `z` decodes cleanly but has
one coefficient in `[γ₁ − β, γ₁]`. This requires crafting the raw bytes by
hand, since a correct signer will never emit such a `z`.

---

### G2. `HintBitUnpack` "decreasing cumulative count" branch is dead — 1 survivor

FIPS 204 Algorithm 21 validates that the per-polynomial cumulative hint count
`y[ω+i]` is non-decreasing. The Wycheproof `InvalidHintsEncoding` vectors
exercise three other malformations (non-monotonic positions, count > ω,
non-zero padding) but **never** a decreasing count.

| Surviving mutant          | Effect                                                                      |
|---------------------------|-----------------------------------------------------------------------------|
| `encode.rs:77` `‖` → `&&` | `end < index ‖ end > ω` → conjunction; each half can only fire with the other. |

The "too many hints" and "omega+1 hints" vectors both set `end > ω`, but the
mutated `&&` makes that half dead. The mutant survives because the
out-of-range position bytes happen to violate the *strict-ordering* check
first (positions[66] = 0 ≤ 252, positions[80] = 21 ≤ 254), so `MalformedHint`
is still returned via a different path.

**Suggested vector:** A hint encoding where `y[ω+i] < y[ω+i−1]` (count goes
*down*). All other validity conditions should hold so the cumulative-count
check is the *only* thing that fails.

---

### G3. `UseHint` never sees `r₀ = 0` with a set hint — 1 survivor

FIPS 204 Algorithm 40 branches on `r₀ > 0` versus `r₀ ≤ 0` when a hint bit is
set. Across 9 285 hint-set positions in valid signatures, 4 684 have `r₀ > 0`
and 4 601 have `r₀ < 0`. **None** have `r₀ = 0`.

| Surviving mutant                   | Effect                           |
|------------------------------------|----------------------------------|
| `hint.rs:31` `r₀ > 0` → `r₀ >= 0`  | Swaps branch when r₀ exactly 0.  |

**Suggested vector:** A valid signature where `w'_approx[i][j]` lands exactly
on a multiple of `2γ₂` at some position where `h[i][j] = 1`. Requires search
over nonces during signing.

---

### G4. `Decompose` edge case never hit with input exactly `q − 1` and hint set — 2 survivors

The `r₊ − r₀ = q − 1` edge case (FIPS 204 Algorithm 36, lines 5–7) fires 7 374
times overall, 184 times at a hint-set position. But the sub-case where the
*input* is exactly `q − 1` (so the pre-adjustment `r₀` is 0) never coincides
with a set hint.

| Surviving mutant                | Effect                                              |
|---------------------------------|-----------------------------------------------------|
| `hint.rs:19` `r₀ − 1` → `r₀ + 1` | Returns `+1` instead of `−1` when pre-adj `r₀ = 0`. |
| `hint.rs:19` `r₀ − 1` → `r₀ / 1` | Returns `0` instead of `−1` when pre-adj `r₀ = 0`.  |

Both mutations only change behaviour when the edge-case input is exactly
`q − 1` *and* the hint bit is set there — so `UseHint` would take a different
branch.

**Suggested vector:** A valid signature where `w'_approx[i][j] = q − 1` at a
position with `h[i][j] = 1`. Search-based as above.

---

### G5. `RejNTTPoly` never samples the exact value `q` — 1 survivor

The rejection sampler accepts 23-bit values `< q`. Across all matrix
expansions, the value `q = 8 380 417` never appears exactly in the SHAKE128
stream.

| Surviving mutant                  | Effect                                      |
|-----------------------------------|---------------------------------------------|
| `sample.rs:36` `v < Q` → `v <= Q` | Accepts `q` as valid (equivalent to 0 mod q). |

**Suggested vector:** A `ρ` seed (searched) such that for some `(row, col)`
pair, `SHAKE128(ρ ‖ col ‖ row)` outputs the three-byte pattern that decodes to
exactly `q`. Probability ≈ `1 / 2²³` per sample, so a handful of thousand
random `ρ` values should surface one.

---

## Equivalent mutants — 8 survivors, no action needed

These survive because the implementation has a second normalisation step that
corrects the mutation's damage, or because the mutation is a no-op for the
specific operands.

| Mutant                                                   | Why equivalent                                                                 |
|----------------------------------------------------------|--------------------------------------------------------------------------------|
| `poly.rs:19` `Poly::reduce` → `()` (no-op)               | Signed `z` coefficients propagate through NTT; `decompose`'s `rem_euclid` fixes them. |
| `poly.rs:56` `% Q` → `+ Q` in `shift_left`               | `t1` coeffs × 2¹³ are already < q; subsequent `mul_mod` applies `% Q` anyway. |
| `poly.rs:156` `< 0` → `> 0` / `<= 0` in `sub_mod` (×2)   | Out-of-range values get fixed by `decompose`'s `rem_euclid` or `mul_mod`'s `%`. |
| `hint.rs:28` `(Q−1)` → `(Q+1)` / `(Q/1)` (×2)            | Integer-divided by `2γ₂`, all three give the same quotient (44 or 16).        |
| `sample.rs:35` `\|` → `^` (×2)                           | OR vs XOR on disjoint bit-shifted values produces identical results.          |

These indicate defensive coding in the implementation, not test-vector gaps.

---

## Reproduction

```bash
cargo test                        # All 239 Wycheproof vectors pass.
cargo test --test path_coverage -- --nocapture   # Path-hit counters.
cargo mutants --jobs 4 --timeout 30              # Full mutation sweep (~2 min).
```

Raw mutation output is in `mutants.out/missed.txt` and friends.
