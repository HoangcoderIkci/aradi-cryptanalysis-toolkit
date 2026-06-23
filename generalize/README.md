# Module ②: AABB / cube across the Toffoli · χ family

Does the low-degree behaviour behind ARADI's AABB cube-distinguisher generalise to the χ-based
permutations Ascon and Keccak? **Full comparison: [`COMPARISON.md`](COMPARISON.md).** Design and
honest scope: [`PLAN.md`](PLAN.md).

> **Scope:** reduced-round *structural distinguishers / degree measurements* only. This does
> **not** attack Ascon or Keccak (both NIST standards); they are permutations, so this is not
> key recovery. Cross-validated against the Ascon specification.

## Results

**Milestone 1 — S-box degree** (`sbox_analysis.py`): ARADI (Toffoli) **3**, χ (Ascon/Keccak) **2**.

**Milestone 2a — Ascon permutation** (`ascon_perm.py`): S-box LUT matches the published Ascon
S-box, all components invertible. Component-level verification (no byte-level KAT this session;
the degree study is invariant under round constants and bit order).

**Milestones 2b + 3 — algebraic degree, two-sided bounds** (`degree_growth.py`,
`ascon_milp.py`, `degree_bounds.py`):

| r | ARADI lower | ARADI upper | Ascon lower | Ascon upper |
|---|---|---|---|---|
| 1 | 3 | **3** | 2 | **2** |
| 2 | 6 | **8** | 4 | **4** |
| 3 | 9 | **22** | 8 | **8** |

Upper bounds via bit-based division-property MILP; the Ascon column reproduces the Ascon v1.2
design-spec degree bound (Table 16: deg ≤ 2^r for r ≤ 8). **Ascon: lower = upper = 2^r → degree
exact (2,4,8)** (exactness is the toolkit's two-sided result, not proven by the spec). **ARADI:
higher, with a gap → its degree-3 S-box raises the degree faster per round.**

**Milestone 3b — integral / zero-sum** (`zerosum_demo.py`): a cube of dimension 2^r + 1 zero-sums
after r rounds of Ascon (r = 1,2,3), the family-wide consequence of low degree.

## What transfers

The link *low S-box degree → cube / zero-sum distinguisher* is shared across the Toffoli/χ family.
ARADI's **AABB byte-equality** form, however, depends on its four-word layout and 16-bit-half
linear layer, so it is ARADI-specific; Ascon's analogue is a plain zero-sum integral. Details and
references in [`COMPARISON.md`](COMPARISON.md).

## Files

```
sbox_analysis.py   S-box algebraic degrees (M1)
ascon_perm.py      verified Ascon permutation (M2a)
degree_growth.py   degree lower bounds via restriction (M2b)
ascon_milp.py      Ascon degree upper bounds, division-property MILP (M3) — vs Ascon v1.2 Table 16
degree_bounds.py   combined two-sided table (M2b + M3)
zerosum_demo.py    integral / zero-sum distinguisher on reduced Ascon (M3b)
PLAN.md  COMPARISON.md
```
