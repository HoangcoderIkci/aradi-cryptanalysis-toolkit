"""
6-round MILP attempt for ARADI under PuLP + CBC.

Two cubes are tried:
  * WonlyDim28        -- I_W = {0,...,27}, others empty (cube dim 28)
  * AllWordsDim112    -- I_W = I_X = I_Y = I_Z = {0,...,27} (cube dim 112)

For each cube we attempt the degree bound of X^6_0 and Z^6_0 with a hard
per-call wall-clock budget.  When the budget is exceeded the script falls
back to the CBC incumbent value (best feasible solution found so far) if
one is available, otherwise reports "no result".

NOTE.  Kim et al. (Table 1 of IACR ePrint 2024/1772) get the 6-round
distinguisher with data complexity 2^77 using Gurobi 11 and the modified
"Algorithm 3" of their paper, which performs polynomial expansion in the
penultimate round.  Reproducing this exact bound with vanilla Algorithm 1
+ free CBC is computationally very hard; the goal here is to document the
attempt and report whatever upper bound the solver returns within the
available compute budget.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

from aradi_milp import compute_degree_upper_bound

CUBES = [
    ("WonlyDim28", dict(cube_indices_w=set(range(28)))),
    ("AllWordsDim112", dict(
        cube_indices_w=set(range(28)),
        cube_indices_x=set(range(28)),
        cube_indices_y=set(range(28)),
        cube_indices_z=set(range(28)),
    )),
]

# A *short* per-call wall budget so the script terminates reliably.  6-round
# CBC may produce an incumbent before this, or it may not converge at all.
TIME_LIMIT = 120

SAMPLE_BITS = [('X', 0), ('Z', 0)]

RESULTS_FILE = Path(__file__).with_name("results_6round.txt")


def main():
    lines = []
    def out(s=""):
        print(s); sys.stdout.flush(); lines.append(s)

    out("=" * 70)
    out("3SBDP MILP -- 6-round ARADI degree upper bounds (vanilla Algorithm 1)")
    out("=" * 70)
    out(f"Per-call CBC wall budget: {TIME_LIMIT} s")
    out(f"Sample target bits      : {SAMPLE_BITS}")
    out("")
    out("NOTE: Kim et al. use Gurobi 11 + Algorithm 3 (polynomial-expansion)")
    out("      to obtain the dim-77 distinguisher reported in their Table 1.")
    out("      The vanilla Algorithm 1 + free CBC reproduction below may")
    out("      either time out or give a looser bound; results are reported")
    out("      verbatim from the solver.")
    out("")

    total_start = time.time()
    for cube_name, cube_args in CUBES:
        cube_dim = sum(len(s) for s in cube_args.values())
        out(f"=== Cube {cube_name} -- total dimension {cube_dim} ===")
        for word, bit in SAMPLE_BITS:
            t0 = time.time()
            try:
                ub, t, nv = compute_degree_upper_bound(
                    num_rounds=6,
                    **cube_args,
                    target_bit_position=bit,
                    target_word=word,
                    time_limit=TIME_LIMIT,
                )
            except Exception as exc:
                wall = time.time() - t0
                out(f"  6-round {word}^6_{bit:>2}  cube_dim={cube_dim:3d}  "
                    f"ERROR ({type(exc).__name__}: {exc})  wall={wall:.1f}s")
                continue
            wall = time.time() - t0
            verdict = "BALANCED" if 0 <= ub < cube_dim else "NOT proved"
            out(f"  6-round {word}^6_{bit:>2}  cube_dim={cube_dim:3d}  "
                f"ub={ub}  t={t:.1f}s  wall={wall:.1f}s  nvars={nv}  {verdict}")
        out("")
        # If total time is getting absurd, stop early.
        if time.time() - total_start > 1200:
            out("(skipping remaining cubes -- 20 min wall budget reached)")
            break

    out("=" * 70)
    out("End of 6-round attempt.")
    out("=" * 70)
    RESULTS_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nSaved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
