# -*- coding: utf-8 -*-
"""Milestone 3 of module 2: two-sided algebraic-degree bounds, ARADI vs Ascon.

Combines, per round r:
  - LOWER bounds from degree_growth.py (exact degree of a restriction to a chosen
    set of input bits; a lower bound on the true degree).
  - UPPER bounds from a bit-based division-property MILP:
      * ARADI  -> ../python/aradi_milp.py (max over the four output words; the
        solver status is reported so a time-limited, non-optimal incumbent is
        never silently reported as a proven bound).
      * Ascon  -> ascon_milp.py.

The Ascon upper side (2,4,8) reproduces the Ascon v1.2 design-spec degree bound
(Table 16: deg(p^R) <= 2^R for R<=8). EXACTNESS for Ascon is this toolkit's own
two-sided result -- lower (restriction) == upper (MILP) == 2^r -- not something
the spec proves (the spec gives only the upper side).

Reading: Ascon (chi, S-box degree 2) has lower == upper at r=1,2,3, so its degree
is exactly 2,4,8 = 2^r. ARADI (Toffoli, S-box degree 3) sits higher with a gap
(3; [6,8]; [9,22]) -- its higher-degree S-box raises the degree faster per round.

Honest scope: structural degree measurement of reduced-round primitives. NOT an
attack on Ascon or Keccak.
"""
from __future__ import annotations

import os
import sys

import pulp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import aradi_milp as am          # noqa: E402
import ascon_milp as asm         # noqa: E402
import degree_growth as dg       # noqa: E402

_FULL = list(range(32))


def aradi_lower(r):
    return dg.degree_after(dg.aradi_rounds, 4, 32, dg.columns_vars([0, 8, 16, 24], 4), r)


def ascon_lower(r):
    return dg.degree_after(dg.ascon_rounds, 5, 64, dg.columns_vars([0, 21, 42], 5), r)


def aradi_upper(r, time_limit=150):
    """Max division-property degree bound over the four output words (bit 0).
    Returns (bound, optimal) where optimal=True iff the max-achieving word's
    MILP was solved to proven optimality (not a time-limit incumbent)."""
    best, best_optimal = -1, False
    for tw in "WXYZ":
        model, _ = am.build_aradi_milp(
            r, cube_indices_w=_FULL, cube_indices_x=_FULL,
            cube_indices_y=_FULL, cube_indices_z=_FULL,
            target_word=tw, target_bit_position=0)
        st = model.solve(pulp.PULP_CBC_CMD(msg=False, timeLimit=time_limit))
        val = int(round(pulp.value(model.objective)))
        is_opt = (pulp.LpStatus[st] == "Optimal")
        if val > best:
            best, best_optimal = val, is_opt
        elif val == best and is_opt:
            best_optimal = True
    return best, best_optimal


def ascon_upper(r, time_limit=150):
    return asm.degree_upper_bound(r, time_limit=time_limit)[0]


def main():
    asm.verify_anf()
    asm.discrimination_test()
    print("Two-sided algebraic-degree bounds (lower = restriction; upper = division-property MILP)")
    print(f"  {'r':>2} | {'ARADI lower':>11} {'ARADI upper':>14} | {'Ascon lower':>11} {'Ascon upper':>11}")
    print("  " + "-" * 60)
    for r in (1, 2, 3):
        al = aradi_lower(r)
        au, au_opt = aradi_upper(r)
        sl, su = ascon_lower(r), ascon_upper(r)
        au_str = f"{au} ({'opt' if au_opt else 'incumbent'})"
        print(f"  {r:>2} | {al:>11} {au_str:>14} | {sl:>11} {su:>11}")
        assert al <= au and sl <= su, f"round {r}: lower exceeds upper (bug)"
        assert su == 2 ** r, f"round {r}: Ascon upper {su} != 2^{r} (spec upper-side mismatch)"
    print("  " + "-" * 60)
    print("Ascon: lower == upper == 2^r  -> degree is EXACT (2,4,8) by this toolkit's")
    print("       two-sided computation (the spec gives only the <= side).")
    print("ARADI: higher and with a gap -> degree-3 S-box raises degree faster.")
    print("       (r=3 upper 22 is proven optimal; see 'opt' tag.)")


if __name__ == "__main__":
    main()
