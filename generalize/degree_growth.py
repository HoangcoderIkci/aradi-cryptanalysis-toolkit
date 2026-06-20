# -*- coding: utf-8 -*-
"""Milestone 2b of module 2: algebraic-degree growth per round, ARADI vs Ascon.

We measure the EXACT algebraic degree of the round function output *restricted to
a chosen set of input bits* (all other state bits fixed to 0). That restricted
degree is a rigorous LOWER BOUND on the true degree of the permutation.

Choice of variables matters: an S-box mixes one bit from each word in a column,
so to avoid an artificially low degree we take, per chosen column, the bit from
EVERY word (5 words for Ascon, 4 for ARADI). Then the per-round multiplier
equals the S-box degree until variables/saturation cap it.

Method per (cipher, rounds r): evaluate the r-round map over all 2^d assignments
of the d variable bits, pack each output state into one big integer, run the
Mobius (ANF) transform over the whole output vector at once (big-int XORs), and
report the maximum monomial degree with a non-zero coefficient.

Sanity anchors (asserted): r=0 -> degree 1 (identity is linear); r=1 -> degree
equal to the S-box degree (Ascon chi = 2, ARADI Toffoli = 3).

Honest scope: reduced-round structural measurement, lower bounds. NOT an attack
on Ascon/Keccak. See PLAN.md.
"""
from __future__ import annotations

import os
import sys

import ascon_perm

# Reuse the verified ARADI round (python/aradi_ref.py) for an apples-to-apples comparison.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from aradi_ref import sbox_layer as aradi_sbox, linear_layer as aradi_linear  # noqa: E402


def _mobius_degree(values: list[int]) -> int:
    """Max algebraic degree over all output bits, given the output vector (one big
    int per input assignment). In-place Mobius butterfly on big ints, then scan."""
    a = list(values)
    n = len(a)
    step = 1
    while step < n:
        for base in range(0, n, step * 2):
            for k in range(base, base + step):
                a[k + step] ^= a[k]
        step *= 2
    deg = 0
    for monomial, coeff in enumerate(a):
        if coeff:
            deg = max(deg, bin(monomial).count("1"))
    return deg


def _build_state(nwords: int, var_positions: list[tuple[int, int]], assignment: int) -> list[int]:
    s = [0] * nwords
    for i, (word, bit) in enumerate(var_positions):
        if (assignment >> i) & 1:
            s[word] |= 1 << bit
    return s


def degree_after(eval_rounds, nwords: int, wordbits: int,
                 var_positions: list[tuple[int, int]], r: int) -> int:
    d = len(var_positions)
    packed = []
    for assignment in range(1 << d):
        out = eval_rounds(_build_state(nwords, var_positions, assignment), r)
        oi = 0
        for k, word in enumerate(out):
            oi |= word << (wordbits * k)
        packed.append(oi)
    return _mobius_degree(packed)


# ---- round functions (keyless: key/constant additions are affine, degree-neutral) ----

def ascon_rounds(state: list[int], r: int) -> list[int]:
    return ascon_perm.permutation(state, r)


def aradi_rounds(state: list[int], r: int) -> list[int]:
    w, x, y, z = state
    for i in range(r):
        w, x, y, z = aradi_sbox(w, x, y, z)
        w, x, y, z = aradi_linear(w, x, y, z, i)
    return [w, x, y, z]


def columns_vars(cols: list[int], nwords: int) -> list[tuple[int, int]]:
    """All `nwords` bits of each chosen column -> variable positions."""
    return [(word, c) for c in cols for word in range(nwords)]


def run(name: str, eval_rounds, nwords: int, wordbits: int, cols: list[int],
        rounds: int, sbox_deg: int) -> None:
    vars_ = columns_vars(cols, nwords)
    print(f"\n{name}: {len(vars_)} variable bits ({nwords} words x columns {cols}), "
          f"S-box degree {sbox_deg}")
    print(f"  {'rounds':>6} | {'degree (lower bound)':>22}")
    print("  " + "-" * 33)
    degs = []
    for r in range(rounds + 1):
        deg = degree_after(eval_rounds, nwords, wordbits, vars_, r)
        degs.append(deg)
        note = ""
        if r >= 1 and degs[r - 1]:
            note = f"  (x{deg / degs[r - 1]:.2f} vs prev)" if degs[r - 1] else ""
        print(f"  {r:>6} | {deg:>22}{note}")
    assert degs[0] == 1, f"{name}: r=0 should be degree 1 (identity), got {degs[0]}"
    assert degs[1] == sbox_deg, f"{name}: r=1 should equal S-box degree {sbox_deg}, got {degs[1]}"
    print(f"  sanity OK: r0=1, r1={sbox_deg} (= S-box degree)")


def main() -> None:
    print("Algebraic-degree growth. We report the EXACT degree of the output restricted")
    print("to the chosen input bits, which is a LOWER BOUND on the true degree.")
    # ARADI: 4-bit Toffoli S-box (deg 3). Columns spread across the 32-bit words.
    run("ARADI round (Toffoli)", aradi_rounds, nwords=4, wordbits=32,
        cols=[0, 8, 16, 24], rounds=3, sbox_deg=3)
    # Ascon: 5-bit chi S-box (deg 2). Columns spread across the 64-bit words.
    run("Ascon permutation (chi)", ascon_rounds, nwords=5, wordbits=64,
        cols=[0, 21, 42], rounds=3, sbox_deg=2)
    print("\nReading (honest): at round 1 the degree equals the S-box degree exactly")
    print("(ARADI 3 vs Ascon 2). Beyond that, ARADI's higher-degree Toffoli S-box drives")
    print("the lower bound up faster (1,3,6,9) than Ascon's degree-2 chi (1,2,4,8). These")
    print("are LOWER bounds (restriction to the chosen input bits), not a fixed per-round")
    print("multiplier; tight upper bounds (division property / MILP) are Milestone 3.")
    print("Reduced-round structural measurement only - NOT an attack on Ascon/Keccak.")


if __name__ == "__main__":
    main()
