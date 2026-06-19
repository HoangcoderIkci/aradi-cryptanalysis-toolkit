# Module ②: generalising AABB/cube to the Toffoli / χ family — plan

## Hypothesis

ARADI's AABB cube-distinguisher exists because its S-box (4 Toffoli gates) has **low
algebraic degree** (3). Ascon and Keccak use the **χ (chi)** non-linear map, which is also
low-degree (2) and AND-based — structurally close to a Toffoli gate. **Question:** does an
analogous *low-degree, structured cube-sum* distinguisher appear in χ-based permutations at
**reduced rounds**, i.e. is the AABB phenomenon a feature of the **Toffoli/χ family** rather
than of ARADI alone?

## GATE-0 — honest assessment (read before investing weeks)

Ascon and Keccak are **heavily studied** and already have public distinguishers:
- Ascon: cube / cube-like attacks (Dobraunig–Eichlseder–Mendel–Schläffer and follow-ups).
- Keccak: zero-sum distinguishers (Aumasson–Meier; Boura–Canteaut–De Cannière).

So the realistic contribution of ② is **NOT** a new attack and **NOT** beating those results.
The defensible contribution is a **unified comparative lens + reproducible tooling**:
> "Does the *specific* AABB-style structured cube-sum, and the cube methodology used on ARADI,
> transfer to χ-based permutations at reduced rounds — measured the same way, in one toolkit?"

That is a **methodology / reproduction** result. It is worth doing for a portfolio (breadth
across a cipher family + reproducibility), but do **not** expect novelty over existing
Ascon/Keccak cryptanalysis. Frame every claim accordingly.

## Honest scope boundary

- Reduced-round **structural distinguishers / degree measurements** only.
- **Never** "break Ascon" or "break Keccak" — both are NIST standards (Ascon = Lightweight
  Crypto 2023; Keccak = SHA-3). Ascon/Keccak are **permutations** (no key schedule), so this
  is **not** key recovery — it measures degree growth and structured cube-sums of the
  permutation.

## Targets

- The **χ map** (5-bit), shared by Ascon and Keccak.
- **Ascon permutation** (320-bit state) at reduced rounds.
- A **small Keccak variant** (Keccak-f[200] or the toy Keccak-f[25]) for tractable
  exhaustive/large-cube experiments.

## Milestones

1. **S-box degree baseline:** algebraic degree of ARADI's S-box vs χ vs the Ascon S-box →
   confirm the "low-degree common to the family" premise. → [`sbox_analysis.py`](sbox_analysis.py).
2a. **Faithful Ascon permutation, verified:** implement Ascon p (S-box + linear + constants
   from the cited spec) and verify it. Verification here = every component is invertible
   (S-box layer is a 5-bit bijection whose recovered LUT matches the published Ascon S-box;
   each per-word linear map has GF(2) rank 64). A byte-level bare-permutation KAT was not
   publicly available this session; it is not needed for 2b because algebraic degree is
   invariant under the affine round constants and under bit re-ordering.
   → [`ascon_perm.py`](ascon_perm.py).
2b. **Degree growth over reduced rounds:** measure algebraic-degree growth per round —
   rigorously (cube/division-property design, reusing [`../python/aradi_milp.py`](../python/aradi_milp.py)
   for upper bounds), labelling lower vs upper bounds honestly. NOT yet done.
3. **Structured cube-sum search:** look for AABB-style structured / zero cube-sums at reduced
   rounds — empirically first, then via the MILP model.
4. **Comparative write-up:** tabulate "after how many rounds does structure vanish" for
   ARADI vs Ascon vs Keccak; relate back to ARADI's AABB.

## Status

- [x] Milestone 1 — S-box degree baseline.
- [x] Milestone 2a — faithful Ascon permutation, components verified.
- [ ] Milestone 2b — degree growth over reduced rounds (rigorous bounds).
- [ ] Milestone 3 — cube-sum search.
- [ ] Milestone 4 — comparative write-up.
