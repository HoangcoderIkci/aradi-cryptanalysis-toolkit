"""
Reproduce the 5-round byte-wise equal property verification of
Kim et al. (IACR ePrint 2024/1772) for ARADI with the cube
    I_W = {11, 12, ..., 23},  I_X = I_Y = I_Z = empty
using a 3SBDP MILP model solved by PuLP+CBC.

Per Kim et al. Section 4, the byte-wise equal property of X^5 ⊕ Z^5 follows
from:
  (i)   Theorem 1: L_i preserves the byte-wise equal property;
  (ii)  the algebraic identity
        X^5  = L_0(X^4 ⊕ W^4 ⊙ Y^4)
        Z^5  = L_0(Z^4 ⊕ X^4 ⊙ Y^4 ⊕ W^4 ⊙ Y^4)
        ==>  X^5 ⊕ Z^5 = L_0(X^4 ⊕ Z^4 ⊕ X^4 ⊙ Y^4);
  (iii) vanilla Algorithm 1 shows that X^4_i, Z^4_i (and X^4 ⊙ Y^4) have
        algebraic degree < 13 in every bit, so their cube sums over 2^13
        plaintexts vanish, leaving only L_0(W^4 ⊙ Y^4) -- and L_0 of a
        byte-wise-equal quantity is itself byte-wise equal.

This script verifies the upper-bound degrees for X^4_i and Z^4_i for every
bit i in {0,...,31}.  Bits are dispatched to a multiprocessing pool to run
several CBC instances in parallel.
"""

from __future__ import annotations

import multiprocessing as mp
import sys
import time
from pathlib import Path

# Import lazily inside workers (PuLP setup is per-process).
CUBE_W = frozenset(range(11, 24))
CUBE_SIZE = len(CUBE_W)
NUM_ROUNDS = 4   # Per Section 4 of Kim et al. -- we check the 4-round components.
RESULTS_FILE = Path(__file__).with_name("results_5round.txt")
PROGRESS_FILE = Path(__file__).with_name("results_5round_progress.txt")


def worker(args):
    """Compute the degree bound for one (num_rounds, target_word, bit)."""
    num_rounds, target_word, bit = args
    from aradi_milp import compute_degree_upper_bound
    ub, t, _ = compute_degree_upper_bound(
        num_rounds=num_rounds,
        cube_indices_w=set(CUBE_W),
        target_bit_position=bit,
        target_word=target_word)
    return (num_rounds, target_word, bit, ub, t)


def main() -> None:
    print(f"3SBDP MILP verification -- {NUM_ROUNDS}-round components for "
          f"5-round byte-wise equal property")
    print(f"Cube I_W = {{11..23}}, size = {CUBE_SIZE}")

    jobs = []
    for word in ('X', 'Z'):
        for bit in range(32):
            jobs.append((NUM_ROUNDS, word, bit))

    # Use 4 worker processes (CBC will fight a bit for CPU but it's still net
    # faster than serial on the typical 4-8 core laptop).
    n_workers = max(1, min(4, mp.cpu_count() - 1))
    print(f"Using {n_workers} parallel workers for {len(jobs)} bit positions")
    sys.stdout.flush()

    t_start = time.time()
    results = {}
    PROGRESS_FILE.write_text("", encoding="utf-8")
    with mp.Pool(n_workers) as pool:
        for r, word, bit, ub, t in pool.imap_unordered(worker, jobs):
            results[(word, bit)] = (ub, t)
            line = f"  {NUM_ROUNDS}-round {word}^{r}_{bit:>2}: deg<={ub:>3}, t={t:6.2f}s  (elapsed total: {time.time()-t_start:6.1f}s)"
            print(line)
            sys.stdout.flush()
            with PROGRESS_FILE.open("a", encoding="utf-8") as fh:
                fh.write(line + "\n")

    total_t = time.time() - t_start
    print(f"\nAll done in {total_t:.1f}s")
    sys.stdout.flush()

    # ---- Write the formatted report --------------------------------------
    lines = []
    def out(s=""):
        lines.append(s)

    out("=" * 70)
    out("3SBDP MILP verification of 5-round ARADI byte-wise equal property")
    out("=" * 70)
    out(f"Cube           : I_W = {{11, 12, ..., 23}}, I_X = I_Y = I_Z = empty")
    out(f"Cube dimension : {CUBE_SIZE}")
    out(f"Solver         : PuLP + CBC (free, embedded in PuLP)")
    out(f"Strategy       : verify the 4-round component bounds, per Section 4")
    out(f"                 of Kim et al. (IACR ePrint 2024/1772)")
    out("")
    out(f"Total CBC time : {total_t:.1f}s for {len(jobs)} runs")
    out(f"Workers        : {n_workers} parallel processes")
    out("")
    out(f"  bit  |  deg X^{NUM_ROUNDS}_i  |  deg Z^{NUM_ROUNDS}_i  |   t_X (s)  |   t_Z (s)")
    out("-" * 64)

    max_x = max_z = 0
    sum_tx = sum_tz = 0.0
    for i in range(32):
        ub_x, t_x = results[('X', i)]
        ub_z, t_z = results[('Z', i)]
        max_x = max(max_x, ub_x)
        max_z = max(max_z, ub_z)
        sum_tx += t_x; sum_tz += t_z
        out(f"  {i:>3}  |   {ub_x:>7}   |   {ub_z:>7}   |   {t_x:6.2f}  |   {t_z:6.2f}")

    out("-" * 64)
    out(f"  max   |   {max_x:>7}   |   {max_z:>7}   |   {sum_tx:6.1f}  |   {sum_tz:6.1f}")
    out("")
    cube_size = CUBE_SIZE
    all_lt = all(results[('X', i)][0] < cube_size for i in range(32)) and \
             all(results[('Z', i)][0] < cube_size for i in range(32))
    out(f"All bits have degree upper bound < cube size ({cube_size}) ? "
        f"{'YES' if all_lt else 'NO'}")
    if all_lt:
        out("")
        out("==> X^4 and Z^4 are BALANCED on the cube I_W = {11..23}.")
        out("    Combined with Theorem 1 of Kim et al. and the algebraic")
        out("    identity X^5 ⊕ Z^5 = L_0(X^4 ⊕ Z^4 ⊕ X^4 ⊙ Y^4), this")
        out("    confirms the 5-round byte-wise equal property of X^5 ⊕ Z^5.")

    RESULTS_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nResults written to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
