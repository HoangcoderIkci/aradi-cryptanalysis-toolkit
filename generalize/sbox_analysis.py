# -*- coding: utf-8 -*-
"""Milestone 1 of module 2: algebraic-degree baseline across the Toffoli / chi family.

Establishes the premise behind generalising ARADI's AABB cube-distinguisher: the
ARADI S-box (4 Toffoli gates) and the chi map used by Ascon / Keccak are all
*low algebraic degree*, which is what makes a cube methodology applicable.

For each S-box we compute, via the Mobius (ANF) transform, the algebraic degree
of every output coordinate and take the maximum. We also check bijectivity as a
sanity test that the implementation is a valid permutation.

Honest scope: this is a structural comparison of building blocks, NOT an attack
on Ascon or Keccak. See PLAN.md.
"""
from __future__ import annotations

from math import comb  # noqa: F401  (kept for readers; not strictly needed)


def anf_coeffs(truth_table: list[int]) -> list[int]:
    """Mobius transform: truth table (length 2^n, GF(2) values) -> ANF coeffs."""
    a = list(truth_table)
    size = len(a)
    step = 1
    while step < size:
        for base in range(0, size, step * 2):
            for k in range(base, base + step):
                a[k + step] ^= a[k]
        step *= 2
    return a


def boolean_degree(truth_table: list[int]) -> int:
    """Algebraic degree of one Boolean function given by its truth table."""
    coeffs = anf_coeffs(truth_table)
    deg = 0
    for monomial, c in enumerate(coeffs):
        if c:
            deg = max(deg, bin(monomial).count("1"))
    return deg


def sbox_degree(lut: list[int], n: int) -> int:
    """Algebraic degree of an n-bit S-box = max degree over its output bits."""
    best = 0
    for bit in range(n):
        tt = [(lut[v] >> bit) & 1 for v in range(1 << n)]
        best = max(best, boolean_degree(tt))
    return best


def is_bijection(lut: list[int]) -> bool:
    return sorted(lut) == list(range(len(lut)))


# ---------------- S-box definitions ----------------

def aradi_sbox_lut() -> list[int]:
    """ARADI 4-bit S-box (4 Toffoli gates), mirrors python/aradi_ref.py:sbox_layer."""
    lut = []
    for v in range(16):
        w, x, y, z = (v >> 3) & 1, (v >> 2) & 1, (v >> 1) & 1, v & 1
        x ^= w & y
        z ^= x & y
        y ^= w & z
        w ^= x & z
        lut.append((w << 3) | (x << 2) | (y << 1) | z)
    return lut


def chi_lut() -> list[int]:
    """Keccak/Ascon chi map on a 5-bit row: y_i = x_i XOR ((NOT x_{i+1}) AND x_{i+2})."""
    lut = []
    for v in range(32):
        x = [(v >> (4 - i)) & 1 for i in range(5)]
        y = [x[i] ^ ((x[(i + 1) % 5] ^ 1) & x[(i + 2) % 5]) for i in range(5)]
        out = 0
        for bit in y:
            out = (out << 1) | bit
        lut.append(out)
    return lut


def ascon_sbox_lut() -> list[int]:
    """Ascon 5-bit S-box = affine o chi o affine (canonical bit-sliced form).

    Affine wrappers do not change algebraic degree, so this stays degree 2 like chi.
    """
    nt = lambda b: b ^ 1
    lut = []
    for v in range(32):
        x0, x1, x2, x3, x4 = (v >> 4) & 1, (v >> 3) & 1, (v >> 2) & 1, (v >> 1) & 1, v & 1
        x0 ^= x4; x4 ^= x3; x2 ^= x1
        t0 = nt(x0) & x1; t1 = nt(x1) & x2; t2 = nt(x2) & x3
        t3 = nt(x3) & x4; t4 = nt(x4) & x0
        x0 ^= t1; x1 ^= t2; x2 ^= t3; x3 ^= t4; x4 ^= t0
        x1 ^= x0; x0 ^= x4; x3 ^= x2; x2 = nt(x2)
        lut.append((x0 << 4) | (x1 << 3) | (x2 << 2) | (x3 << 1) | x4)
    return lut


def main() -> None:
    rows = [
        ("ARADI S-box (4 Toffoli gates)", aradi_sbox_lut(), 4, 3),
        ("chi map (Keccak / Ascon core)", chi_lut(), 5, 2),
        ("Ascon S-box (affine o chi)",    ascon_sbox_lut(), 5, 2),
    ]
    print(f"{'S-box':<32} {'bits':>4} {'bijective':>10} {'alg. degree':>12} {'expected':>9}")
    print("-" * 72)
    all_ok = True
    for name, lut, n, expected in rows:
        bij = is_bijection(lut)
        deg = sbox_degree(lut, n)
        ok = bij and deg == expected
        all_ok = all_ok and ok
        print(f"{name:<32} {n:>4} {str(bij):>10} {deg:>12} {expected:>9}  {'OK' if ok else 'MISMATCH'}")
    print("-" * 72)
    print("Premise: ARADI (Toffoli, deg 3) and the Ascon/Keccak chi family (deg 2) are")
    print("all LOW algebraic degree -> a cube methodology applies across the family.")
    print("Module 2 then asks whether the *structured* AABB-style cube-sum also transfers")
    print("(see PLAN.md). This file does NOT attack Ascon or Keccak.")
    assert all_ok, "a sanity check failed (bijection or degree mismatch)"
    print("\nALL CHECKS PASS")


if __name__ == "__main__":
    main()
