# Module ②: AABB / cube across the Toffoli · χ family

Does the AABB-style structured cube-sum that ARADI exhibits generalise to χ-based
permutations (Ascon, Keccak)? Full design, honest scope and milestones: **[`PLAN.md`](PLAN.md)**.

> **Scope:** reduced-round *structural distinguishers / degree measurements* only. This does
> **not** attack Ascon or Keccak (both NIST standards), and they are permutations — so this
> is not key recovery. See [`PLAN.md`](PLAN.md) §"Honest scope".

## Milestone 1 — algebraic-degree baseline ✅

`python sbox_analysis.py` (verified, bijection + degree checks pass):

| S-box | bits | algebraic degree |
|---|---|---|
| ARADI (4 Toffoli gates) | 4 | **3** |
| χ (Keccak / Ascon core) | 5 | **2** |
| Ascon S-box (affine ∘ χ) | 5 | **2** |

→ The whole family is **low algebraic degree**, which is the premise that makes a cube
methodology applicable across it. Next: degree growth over reduced rounds (Milestone 2).
