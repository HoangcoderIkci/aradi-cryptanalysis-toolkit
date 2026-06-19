# -*- coding: utf-8 -*-
"""Milestone 2a of module 2: a faithful Ascon permutation, with verification.

Ascon permutation p (320-bit state = five 64-bit words x0..x4). One round =
  pC  add a 1-byte round constant to x2,
  pS  the 5-bit S-box applied bit-sliced across the 64 columns,
  pL  per-word linear diffusion (xor of two right-rotated copies).

Parameters are taken from authoritative sources (official Ascon spec +
NIST SP 800-232, cross-checked against the OpenTitan Ascon documentation):
  - linear rotations: x0(19,28) x1(61,39) x2(1,6) x3(10,17) x4(7,41)
  - p12 round constants: f0 e1 d2 c3 b4 a5 96 87 78 69 5a 4b
    (a reduced p^a uses the LAST a constants; e.g. p6 uses 96..4b)
  - S-box: the canonical bit-sliced (affine o chi o affine) form, already
    verified in sbox_analysis.py (5-bit bijection, algebraic degree 2).

VERIFICATION LEVEL (honest): we verify every component is invertible — the
S-box layer (a 5-bit bijection) and each per-word linear map (GF(2) rank 64) —
so the whole permutation is a bijection, and the parameters match the cited
spec. We do NOT match a byte-level known-answer test here (no authoritative
bare-permutation KAT was available this session). For the degree-growth study
in 2b this is sufficient: round constants are affine and bit ordering is a
relabelling, so neither affects algebraic degree — which depends only on the
(verified) S-box and the (cited) linear layer.
"""
from __future__ import annotations

MASK64 = (1 << 64) - 1
ROT = [(19, 28), (61, 39), (1, 6), (10, 17), (7, 41)]
RC = [0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87, 0x78, 0x69, 0x5a, 0x4b]


def rotr(x: int, n: int) -> int:
    return ((x >> n) | (x << (64 - n))) & MASK64


def _not(x: int) -> int:
    return x ^ MASK64


def sbox_layer(s: list[int]) -> list[int]:
    x0, x1, x2, x3, x4 = s
    x0 ^= x4; x4 ^= x3; x2 ^= x1
    t0 = _not(x0) & x1
    t1 = _not(x1) & x2
    t2 = _not(x2) & x3
    t3 = _not(x3) & x4
    t4 = _not(x4) & x0
    x0 ^= t1; x1 ^= t2; x2 ^= t3; x3 ^= t4; x4 ^= t0
    x1 ^= x0; x0 ^= x4; x3 ^= x2; x2 = _not(x2)
    return [x0, x1, x2, x3, x4]


def linear_layer(s: list[int]) -> list[int]:
    return [s[i] ^ rotr(s[i], ROT[i][0]) ^ rotr(s[i], ROT[i][1]) for i in range(5)]


def permutation(state: list[int], rounds: int = 12) -> list[int]:
    """Apply `rounds` rounds of the Ascon permutation (uses the last `rounds` constants)."""
    s = list(state)
    for c in RC[12 - rounds:]:
        s[2] ^= c
        s = sbox_layer(s)
        s = linear_layer(s)
    return [w & MASK64 for w in s]


# ---------------- verification ----------------

def _sbox_lut() -> list[int]:
    """Recover the 5-bit S-box LUT from the bit-sliced layer (put each column bit in lane 0)."""
    lut = []
    for v in range(32):
        s = [((v >> (4 - i)) & 1) for i in range(5)]  # x0..x4 each 0/1 in lane 0
        o = sbox_layer(s)
        lut.append(sum((o[i] & 1) << (4 - i) for i in range(5)))
    return lut


def _gf2_rank64(cols: list[int]) -> int:
    """Rank over GF(2) of 64 column vectors (each a 64-bit int)."""
    basis = []
    for c in cols:
        for b in basis:
            c = min(c, c ^ b)
        if c:
            basis.append(c)
            basis.sort(reverse=True)
    return len(basis)


def verify() -> None:
    # 1. S-box layer is a 5-bit bijection.
    lut = _sbox_lut()
    assert sorted(lut) == list(range(32)), "S-box layer is not a bijection"
    # 2. Each per-word linear map x -> x ^ rotr(x,r1) ^ rotr(x,r2) is invertible (rank 64).
    for i, (r1, r2) in enumerate(ROT):
        cols = [(1 << k) ^ rotr(1 << k, r1) ^ rotr(1 << k, r2) for k in range(64)]
        assert _gf2_rank64(cols) == 64, f"linear map for word {i} is singular"
    # 3. Therefore the full permutation is a bijection (each component invertible).
    print("Ascon permutation VERIFIED (components invertible):")
    print(f"  S-box layer : 5-bit bijection  (LUT = {[hex(v) for v in lut]})")
    print(f"  linear maps : all rank 64 (invertible) for rotations {ROT}")
    print(f"  constants   : p12 = {[hex(c) for c in RC]}  (reduced p^a uses last a)")
    print("  NOTE: parameters from cited spec; byte-level KAT match deferred to 2b "
          "(does not affect the affine-invariant degree study).")


if __name__ == "__main__":
    verify()
