# Generalising AABB to the Toffoli / χ cipher family  (roadmap — module ②)

**Question.** The AABB cube-distinguisher in ARADI comes from its Toffoli-gate S-box
(algebraic degree 3, built from 4 Toffoli gates). Is an AABB-style low-degree /
byte-structure property a **structural feature of the whole Toffoli / χ family**, rather
than a quirk of ARADI?

**Plan.**
1. Pick 2–3 other Toffoli-gate / χ-based primitives at a *reduced-round* scope, e.g.
   Simon (AND-based), Xoodoo and Ascon/Keccak's χ layer.
2. Reuse the MILP propagation model in [`../python/aradi_milp.py`](../python/aradi_milp.py)
   to search for analogous low-degree cube distinguishers in each.
3. Tabulate: does a byte/word-structured cube sum appear, after how many rounds, with what
   random-case probability?
4. Write up as a small comparative study: *what the Toffoli/χ structure buys an attacker.*

**Honesty boundary.** This is a **structural comparison**, not an attack on any of those
ciphers. Report distinguishers on reduced rounds only; never claim to break a primitive.

**Status:** not started. Owner: —.
