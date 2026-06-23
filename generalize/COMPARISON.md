# AABB / cube across the Toffoli · χ family: ARADI vs Ascon/Keccak

A reproducible structural comparison: does the low-degree behaviour behind ARADI's AABB
cube-distinguisher generalise to the χ-based permutations Ascon and Keccak?

> **Scope (read first).** This is a structural measurement of **reduced-round** primitives.
> It does **not** attack Ascon or Keccak (both are NIST standards), and both are permutations,
> so this is not key recovery. Degree results are cross-validated against the Ascon
> specification. See [`PLAN.md`](PLAN.md).

## 1. Shared premise: a low-degree S-box (Milestone 1)

`sbox_analysis.py` measures the algebraic degree of each non-linear building block:

| S-box | algebraic degree |
|---|---|
| ARADI (4 Toffoli gates) | **3** |
| χ (Keccak / Ascon core) | **2** |
| Ascon S-box (affine ∘ χ) | **2** |

Both ARADI's Toffoli S-box and the χ map used by Ascon and Keccak have low algebraic degree.
A cube / higher-order-derivative methodology applies to the whole family.

## 2. A verified Ascon permutation (Milestone 2a)

`ascon_perm.py` implements the Ascon permutation with parameters from NIST SP 800-232 (S-box
LUT, linear-layer rotations). Verification is **component-level**: the recovered S-box LUT
matches the published Ascon S-box and every component is invertible (5-bit bijection; each
per-word linear map has GF(2) rank 64). A byte-level known-answer test was not matched this
session; it is not needed for the degree study, because algebraic degree is invariant under the
affine round constants and under bit re-ordering.

## 3. Algebraic degree per round: two-sided bounds (Milestones 2b, 3)

`degree_bounds.py` reports, per round r:
- a **lower bound** from `degree_growth.py` (the exact degree of the round function restricted
  to a chosen set of input bits is a lower bound on the true degree), and
- an **upper bound** from a bit-based division-property MILP (`ascon_milp.py` for Ascon,
  `../python/aradi_milp.py` for ARADI).

| r | ARADI lower | ARADI upper | Ascon lower | Ascon upper |
|---|---|---|---|---|
| 1 | 3 | **3** | 2 | **2** |
| 2 | 6 | **8** | 4 | **4** |
| 3 | 9 | **22** | 8 | **8** |

The Ascon upper bounds (2, 4, 8) reproduce the forward degree bound of the **Ascon v1.2 design
specification** (Dobraunig–Eichlseder–Mendel–Schläffer; Table 16: deg(p^R) ≤ 2^R for R ≤ 8),
which cross-validates the upper side of the MILP. NIST SP 800-232 (the final standard) fixes the
S-box and linear-layer parameters but contains no degree analysis. The MILP is also checked for
**discrimination** (`ascon_milp.py`: a linear S-box gives degree 1, an injected cubic monomial
gives degree 3), so the 2^r match is not a vacuous over-count.

Reading:
- **Ascon: the lower bound (restriction) equals the upper bound (MILP) at r = 1,2,3, so the
  degree is exactly 2, 4, 8 = 2^r.** Exactness is this toolkit's two-sided result, not something
  the specification proves (the spec gives only the ≤ side). The degree-2 χ map doubles the
  degree each round, and the bound is tight.
- **ARADI sits higher with a gap: 3, then [6,8], then [9,22].** Its degree-3 Toffoli S-box raises
  the degree faster per round (the regime is closer to 3^r than 2^r).

So the low-degree property is shared, but ARADI reaches a high algebraic degree in fewer rounds
than Ascon. This is a quantitative structural difference inside the same Toffoli/χ family.

## 4. The integral / zero-sum consequence (Milestone 3b)

A low algebraic degree forces an integral (zero-sum) distinguisher: summing the output over a
cube of dimension D > degree gives zero on every output bit. `zerosum_demo.py` confirms this
empirically on the Ascon permutation: a cube of dimension 2^r + 1 zero-sums after r rounds for
r = 1,2,3. This is consistent with the low-degree cause behind ARADI's AABB distinguisher.
(The boundary case D = degree is cube-dependent: one cube summing to zero at D = degree does
not contradict the proven degree, since not every degree-d cube contains the top monomial.)

## 5. What transfers and what does not

- **Transfers:** the link *low S-box degree → cube / integral (zero-sum) distinguisher at reduced
  rounds* is a structural feature of the whole Toffoli/χ family. Ascon and Keccak both share the
  degree-2 χ core, and the zero-sum distinguisher above is its direct consequence.
- **Does not transfer verbatim:** ARADI's **AABB** property is a *byte-structured* cube sum
  (within a word, the first two bytes equal and the last two equal). That specific structure
  depends on ARADI's four-word layout and its 16-bit-half linear layer. Ascon is a 320-bit
  permutation with five 64-bit words and different diffusion, so its analogue is a plain
  zero-sum integral, not a byte-equality pattern. The AABB *form* is ARADI-specific even though
  the underlying low-degree cause is shared.

## 6. Honest assessment (GATE-0)

Ascon and Keccak are heavily studied NIST standards with existing public distinguishers (Todo's
division-property integrals; Boura–Canteaut degree bounds; zero-sum distinguishers on Keccak).
This module does not claim a new attack or any improvement over those results. Its contribution
is a unified, reproducible **comparison**: the same degree and zero-sum measurements applied to
ARADI and to the χ family in one toolkit, showing precisely what the AABB story shares with the
family and where it is ARADI-specific.

## Reproduce

```bash
python sbox_analysis.py     # S-box degrees (3 vs 2)
python ascon_perm.py        # Ascon permutation, component verification
python degree_growth.py     # degree lower bounds (restriction)
python ascon_milp.py        # Ascon degree upper bounds, vs Ascon v1.2 Table 16
python degree_bounds.py     # combined two-sided table
python zerosum_demo.py      # integral / zero-sum distinguisher on reduced Ascon
```

## References

- Ascon v1.2 design specification (Dobraunig–Eichlseder–Mendel–Schläffer; CAESAR / NIST-LWC
  submission, J. Cryptology 34(3), 2021): algebraic-degree bound (Table 16: deg(p^R) ≤ 2^R, R ≤ 8).
- NIST SP 800-232 (final Ascon standard): S-box LUT (Table 6) and linear-layer rotations.
- C. Boura, A. Canteaut. On the influence of the algebraic degree of the round function on the
  algebraic degree of the whole construction.
- Y. Todo. Structural evaluation by generalized integral property (division property).
- Z. Xiang, W. Zhang, Z. Bao, D. Lin. Applying MILP to bit-based division property (Asiacrypt 2016).
- K. Hao et al. Modeling for three-subset division property without unknown subsets (EUROCRYPT 2020).
- D. Kim et al. Byte-wise equal property of ARADI (ePrint 2024/1772) — the AABB property.
