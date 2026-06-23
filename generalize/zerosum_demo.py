# -*- coding: utf-8 -*-
"""Milestone 3b of module 2: integral / zero-sum distinguisher on reduced Ascon.

The low algebraic degree that gives ARADI its AABB cube-distinguisher also gives
Ascon a plain integral (zero-sum) distinguisher: if a cube has dimension D strictly
greater than the algebraic degree after r rounds, then summing the r-round output
over the whole cube is zero on every output bit (Dinur-Shamir / higher-order
derivative). With Ascon's degree = 2^r (verified in ascon_milp.py), a cube of
dimension 2^r + 1 must zero-sum.

We confirm this empirically with the verified permutation in ascon_perm.py.

Honest scope: this reproduces the standard kind of integral distinguisher; it is
NOT an attack on Ascon. The point is the COMPARISON with ARADI (see COMPARISON.md):
the low-degree -> integral-distinguisher link transfers across the Toffoli/chi
family, but ARADI's *byte-structured* AABB form is specific to its word layout.
"""
from __future__ import annotations

import ascon_perm as ap

DEGREE = {1: 2, 2: 4, 3: 8}  # Ascon forward degree (lower=upper, verified)


def cube_sum(rounds: int, positions):
    """XOR-sum the r-round Ascon output over the full cube on `positions`
    (list of (word, bit)); all other state bits fixed to 0. Returns 5 words."""
    acc = [0] * 5
    d = len(positions)
    for a in range(1 << d):
        s = [0] * 5
        for i, (w, b) in enumerate(positions):
            if (a >> i) & 1:
                s[w] |= 1 << b
        out = ap.permutation(s, rounds)
        for w in range(5):
            acc[w] ^= out[w]
    return acc


def spread_positions(d: int):
    """d cube variables spread across all 5 words and successive bit columns."""
    pos = []
    col = 0
    while len(pos) < d:
        for w in range(5):
            if len(pos) < d:
                pos.append((w, col))
        col += 1
    return pos


def main():
    print("Integral (zero-sum) distinguisher on reduced-round Ascon permutation")
    print("(cube dimension D vs algebraic degree deg; D > deg  =>  output sums to 0)")
    print(f"  {'rounds':>6} | {'deg':>3} | {'cube D=deg+1':>12} | {'sum==0 (zero-sum)':>18}")
    print("  " + "-" * 52)
    for r in (1, 2, 3):
        deg = DEGREE[r]
        D = deg + 1
        acc = cube_sum(r, spread_positions(D))
        zero = all(w == 0 for w in acc)
        print(f"  {r:>6} | {deg:>3} | {D:>12} | {str(zero):>18}")
        assert zero, f"round {r}: D={D} cube did NOT zero-sum (unexpected)"
    print("\n  Boundary check (D = deg, should generally NOT zero-sum -> degree is tight):")
    for r in (2, 3):
        deg = DEGREE[r]
        acc = cube_sum(r, spread_positions(deg))
        zero = all(w == 0 for w in acc)
        print(f"    r={r}, D=deg={deg}: sum==0 ? {zero}  "
              f"({'still zero (this cube missed the top monomial)' if zero else 'NON-zero -> a degree-' + str(deg) + ' monomial survives'})")
    print("\nConfirmed: Ascon has integral/zero-sum distinguishers at D = 2^r + 1,")
    print("the family-wide consequence of low algebraic degree. NOT an attack.")


if __name__ == "__main__":
    main()
