"""
Scaled re-run of the Phase-2 multi-cube experiment on 6-round modified ARADI.

Identical methodology to multicube_attack.py (same functions, same fixed seed
0xC0FFEE), but with NUM_TRIALS configurable (default 100) so the professor's
"20 is too few" concern can be addressed with 100 random master keys.

Because the RNG is seeded once with the SAME seed as the original 20-key run,
the first 20 master keys reproduce the original experiment exactly; trials
21..100 are fresh — i.e. this is a strict superset, fully reproducible.

Recovers the 16-bit base half l of the pair (u,l) — NOT the 256-bit master key
(the projection K_base = XOR K_i is non-invertible).

Usage:  python run_multicube_100.py [num_trials]
Output: multicube_results_100.txt  (does NOT clobber multicube_results.txt)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import multicube_attack as mca

OUT_PATH = Path(__file__).with_name("multicube_results_100.txt")

NUM_ROUNDS = 6
SEED = 0xC0FFEE          # SAME seed as the original 20-key run -> reproducible superset
K_VALUES = (1, 2, 3, 5)


def main():
    num_trials = int(sys.argv[1]) if len(sys.argv) > 1 else 100

    lines = []

    def out(s=""):
        print(s)
        sys.stdout.flush()
        lines.append(s)

    out("=" * 76)
    out(" Multi-cube optimisation of phase-1 attack on 6-round modified ARADI")
    out(" (scaled re-run: 100 random master keys)")
    out("=" * 76)
    out(f"Cube starts (first K used per trial): {mca.DEFAULT_CUBE_STARTS}")
    out(f"K values: {', '.join(str(k) for k in K_VALUES)}")
    out(f"Trials per K: {num_trials}")
    out(f"Fixed RNG seed: 0x{SEED:X}  (first 20 trials reproduce the original run)")
    out("Conditions: 100 random master keys, Phase 2 on Python+NumPy, "
        "single thread, Intel x86_64.")
    out("Phase 1 recovers the 16-bit base half l of the pair (u,l); "
        "the 256-bit master key is NOT recovered (K_base projection is non-invertible).")
    out("")
    out("Two benchmark modes are run:")
    out("  MODE A -- Full AABB check (per-candidate FP rate 2^{-32}).")
    out("            Headline recovery claim: l_true isolated, 0 false positives.")
    out("  MODE B -- Half-AABB check (per-candidate FP rate 2^{-16}).")
    out("            Validates the multi-cube 2^{-16K} false-positive law.")
    out("")

    t_wall = time.time()

    results_A = mca.run_bench_mode(
        "Full AABB (4 byte-pairs equal per word)",
        mca.is_byte_wise_equal, K_VALUES, num_trials, NUM_ROUNDS, SEED, out)
    results_B = mca.run_bench_mode(
        "Half-AABB (only upper half of each word)",
        mca.is_half_byte_wise_equal, K_VALUES, num_trials, NUM_ROUNDS, SEED, out)

    out("")
    out("=" * 76)
    out(" COMBINED SUMMARY")
    out("=" * 76)
    out(f"{'K':>3} | {'A: avg |s|':>10} {'A: max |s|':>10} {'A: time s':>10} | "
        f"{'B: avg |s|':>10} {'B: max |s|':>10} {'B: time s':>10}")
    out("-" * 76)
    for K in K_VALUES:
        a_avg = sum(results_A[K]['survivors']) / len(results_A[K]['survivors'])
        a_max = max(results_A[K]['survivors'])
        a_t = sum(results_A[K]['times']) / len(results_A[K]['times'])
        b_avg = sum(results_B[K]['survivors']) / len(results_B[K]['survivors'])
        b_max = max(results_B[K]['survivors'])
        b_t = sum(results_B[K]['times']) / len(results_B[K]['times'])
        out(f"{K:>3} | {a_avg:>10.3f} {a_max:>10d} {a_t:>10.2f} | "
            f"{b_avg:>10.3f} {b_max:>10d} {b_t:>10.2f}")
    out("-" * 76)
    out("")

    # ---- Headline verification block (Mode A, the real attack) ----
    a1 = results_A[1]
    succ = sum(a1['correct'])
    n = len(a1['correct'])
    # In Mode A the only survivor should be the true l (0 false positives).
    total_fp = sum(s - c for s, c in zip(a1['survivors'], a1['correct']))
    avg_t = sum(a1['times']) / len(a1['times'])
    out("HEADLINE (Mode A, K=1 -- full AABB, the deployed attack):")
    out(f"  Trials (random master keys) : {n}")
    out(f"  Successful (u,l) recoveries  : {succ}/{n}")
    out(f"  False positives (total)      : {total_fp}")
    out(f"  Average Phase 2 time         : {avg_t:.2f} s (Python+NumPy, single thread)")
    out(f"  Working memory per cube       : 2^13 x 4 words x 4 bytes = 128 KB "
        "(independent of trial count)")
    verdict = "PASS" if (succ == n and total_fp == 0) else "FAIL"
    out(f"  VERDICT                      : {verdict} "
        f"({'all recovered, 0 FP' if verdict=='PASS' else 'see failures above'})")
    out("")
    out(f"Total wall-clock: {time.time()-t_wall:.1f}s")

    OUT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nSaved to {OUT_PATH}")
    print(f"HEADLINE_VERDICT={verdict} SUCC={succ}/{n} FP={total_fp} AVG_T={avg_t:.2f}")


if __name__ == "__main__":
    main()
