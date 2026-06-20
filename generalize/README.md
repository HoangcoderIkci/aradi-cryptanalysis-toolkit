# Module ②: AABB / cube across the Toffoli · χ family

Does the AABB-style structured cube-sum that ARADI exhibits generalise to χ-based
permutations (Ascon, Keccak)? Full design, honest scope and milestones: **[`PLAN.md`](PLAN.md)**.

> **Scope:** reduced-round *structural distinguishers / degree measurements* only. This does
> **not** attack Ascon or Keccak (both NIST standards), and they are permutations — so this
> is not key recovery. See [`PLAN.md`](PLAN.md) §"Honest scope".

## Milestone 1 — algebraic-degree baseline of the S-boxes ✅

`python sbox_analysis.py` (bijection + degree checks pass):

| S-box | bits | algebraic degree |
|---|---|---|
| ARADI (4 Toffoli gates) | 4 | **3** |
| χ (Keccak / Ascon core) | 5 | **2** |
| Ascon S-box (affine ∘ χ) | 5 | **2** |

→ The whole family is **low algebraic degree** — the premise that makes a cube methodology
applicable across it.

## Milestone 2a — faithful Ascon permutation, verified ✅

`python ascon_perm.py`. Parameters (rotations, round constants) from the official spec +
NIST SP 800-232, cross-checked against OpenTitan. Verified by component invertibility: the
recovered S-box LUT `[0x4,0xb,0x1f,0x14,…]` **matches the published Ascon S-box**, and each
per-word linear map has GF(2) rank 64. (A byte-level KAT was not needed — algebraic degree is
invariant under the affine round constants and bit re-ordering.)

## Milestone 2b — algebraic-degree growth (lower bounds) ✅

`python degree_growth.py` (~5 s). We measure the exact degree of the output restricted to a
chosen set of input bits (all bits of a few state columns), a **lower bound** on the true
degree. Sanity-checked: round 0 → degree 1, round 1 → exactly the S-box degree.

| rounds | ARADI (Toffoli, S-box deg 3) | Ascon (χ, S-box deg 2) |
|---|---|---|
| 0 | 1 | 1 |
| 1 | **3** | **2** |
| 2 | 6 | 4 |
| 3 | 9 | 8 |

→ At round 1 the degree equals the S-box degree (3 vs 2). Beyond that, ARADI's higher-degree
Toffoli S-box drives the lower bound up faster (1,3,6,9) than Ascon's χ (1,2,4,8). These are
**lower bounds**, not a fixed per-round multiplier.

## Next

- **Tight upper bounds** on degree via the division property (MILP), reusing
  [`../python/aradi_milp.py`](../python/aradi_milp.py) — turns the lower bounds above into
  two-sided bounds (Milestone 3).
- **Structured cube-sum search** for an AABB analogue at reduced rounds, then the comparative
  write-up (Milestones 3–4).
