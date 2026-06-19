"""
Classification of 13-bit cubes for the 5-round byte-wise equal property of
ARADI.

Approach (original contribution of the thesis, extending Kim et al. 2024/1772
which considered only the single cube I_W = {11,...,23}):

  Strategy 1 -- contiguous cubes in word W:
      For every starting index k in {0,1,...,19} we test the 13-bit cube
      I_W = {k, k+1, ..., k+12}, I_X = I_Y = I_Z = empty.

  Strategy 2 -- contiguous cubes in the other three words:
      Same 20 starting positions, but the cube lives in word X, Y, or Z
      (with the other three words empty).

  Strategy 3 -- random non-contiguous 13-bit subsets of word W:
      A small sample of 20 uniformly-random 13-element subsets of {0..31}
      placed in word W.

Verification pipeline per cube:
  1) Experimental cube-sum check.  For 3 fresh random 256-bit master keys we
     sum the 5-round ARADI ciphertexts over the 2^13 cube and test whether the
     resulting X-word and Z-word both satisfy the byte-wise equal property
     ("AABB", i.e. high half = rotl16(high half, 8) and similarly for low
     half).  Cubes for which X and Z are byte-wise equal across all three keys
     are flagged as "byte-wise equal".

  2) MILP confirmation.  For a *representative* output bit (X^4_0) and the
     antipodal bit (Z^4_0) we run the 3SBDP MILP of aradi_milp.py and record
     the algebraic-degree upper bound.  A degree upper bound strictly less
     than 13 = |cube| at *every* bit would imply the cube-sum vanishes; we
     check the two sample bits and report the bound.  CBC is given a 60-s
     wall budget per call.

The full output table is written to milp/classification_results.txt.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import random
import secrets
import sys
import time
from pathlib import Path

from aradi_ref import aradi_encrypt, key_schedule, MASK32

CUBE_BITS = 13
NUM_ROUNDS = 5            # cipher rounds used in the experimental check
MILP_ROUNDS = 4           # 4-round component bound feeds the 5-round property
NUM_KEYS = 3              # random keys for experimental verification
MILP_TIME_LIMIT = 60      # seconds per CBC call
NUM_RANDOM = 20           # Strategy 3 cubes

OUT_PATH = Path(__file__).with_name("classification_results.txt")


# ---------------------------------------------------------------------------
# Experimental check helpers
# ---------------------------------------------------------------------------

def has_aabb(word: int) -> bool:
    """Byte-wise equal property: rotl16-by-8 fixes both halves."""
    u = (word >> 16) & 0xFFFF
    l = word & 0xFFFF
    return ((u << 8) | (u >> 8)) & 0xFFFF == u and \
           ((l << 8) | (l >> 8)) & 0xFFFF == l


def cube_sum(cube_indices: list, target_word_letter: str, key, num_rounds: int):
    """
    Sum (XOR) the ciphertext four-tuple over the 2^|cube| inputs where each
    bit of `cube_indices` runs over {0,1} and the cube lives in word
    `target_word_letter` (the other three words are forced to zero).
    Returns the tuple (sum_W, sum_X, sum_Y, sum_Z) of XOR sums.
    """
    rk = key_schedule(key)
    sum_w = sum_x = sum_y = sum_z = 0
    cube_count = 1 << len(cube_indices)
    for c in range(cube_count):
        w_in = x_in = y_in = z_in = 0
        bits = 0
        for k, bit_pos in enumerate(cube_indices):
            if (c >> k) & 1:
                bits |= (1 << bit_pos)
        if target_word_letter == 'W':
            w_in = bits
        elif target_word_letter == 'X':
            x_in = bits
        elif target_word_letter == 'Y':
            y_in = bits
        else:
            z_in = bits
        w_out, x_out, y_out, z_out = aradi_encrypt(
            w_in, x_in, y_in, z_in, rk, num_rounds=num_rounds)
        sum_w ^= w_out
        sum_x ^= x_out
        sum_y ^= y_out
        sum_z ^= z_out
    return sum_w, sum_x, sum_y, sum_z


def experimental_check(cube_indices, target_word_letter, num_keys=NUM_KEYS,
                       num_rounds=NUM_ROUNDS, rng=None):
    """
    For `num_keys` fresh random master keys, compute the cube sum (over the
    cube placed in `target_word_letter`) and test whether both X- and Z-
    output words are byte-wise equal.  Also record the sums of the first key
    for diagnostic display.
    Returns (all_aabb_xz, x_sum, z_sum, all_sums) where all_sums is a list of
    (w,x,y,z) tuples per key.
    """
    if rng is None:
        rng = random.Random()
    all_aabb_xz = True
    all_sums = []
    for _ in range(num_keys):
        key = [rng.getrandbits(32) for _ in range(8)]
        sw, sx, sy, sz = cube_sum(cube_indices, target_word_letter, key, num_rounds)
        all_sums.append((sw, sx, sy, sz))
        if not (has_aabb(sx) and has_aabb(sz)):
            all_aabb_xz = False
    return all_aabb_xz, all_sums


# ---------------------------------------------------------------------------
# MILP confirmation
# ---------------------------------------------------------------------------

def milp_worker(args):
    """One CBC run for (cube_indices, cube_word_letter, target_word, target_bit)."""
    cube_indices, cube_word, target_word, target_bit, rounds, time_limit = args
    from aradi_milp import compute_degree_upper_bound
    kwargs = dict(num_rounds=rounds, target_bit_position=target_bit,
                  target_word=target_word, time_limit=time_limit)
    if cube_word == 'W':
        kwargs['cube_indices_w'] = set(cube_indices)
    elif cube_word == 'X':
        kwargs['cube_indices_x'] = set(cube_indices)
    elif cube_word == 'Y':
        kwargs['cube_indices_y'] = set(cube_indices)
    else:
        kwargs['cube_indices_z'] = set(cube_indices)
    ub, t, _ = compute_degree_upper_bound(**kwargs)
    return (cube_word, tuple(cube_indices), target_word, target_bit, ub, t)


# ---------------------------------------------------------------------------
# Cube enumerators
# ---------------------------------------------------------------------------

def contiguous_cubes(word_size: int = 32, cube_size: int = CUBE_BITS):
    """Return all contiguous cubes {k,k+1,...,k+12} for k in 0..(32-13)."""
    return [list(range(k, k + cube_size)) for k in range(word_size - cube_size + 1)]


def random_cubes(n: int, seed: int = 0xA8AD1):
    """Return n random 13-element subsets of {0..31}, deterministic via seed."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        sample = rng.sample(range(32), CUBE_BITS)
        sample.sort()
        out.append(sample)
    return out


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    lines = []
    def out(s=""):
        print(s); sys.stdout.flush()
        lines.append(s)

    out("=" * 76)
    out(" Classification of 13-bit cubes -- 5-round ARADI byte-wise equal property")
    out("=" * 76)
    out(f"Cube size           : {CUBE_BITS}")
    out(f"Cipher rounds       : {NUM_ROUNDS} (experimental), {MILP_ROUNDS} (MILP)")
    out(f"Random keys per cube: {NUM_KEYS}")
    out(f"MILP wall budget    : {MILP_TIME_LIMIT}s per CBC call")
    out("")

    # ---- Build candidate list ---------------------------------------------
    candidates = []
    for cube in contiguous_cubes():
        for word in ('W', 'X', 'Y', 'Z'):
            candidates.append(('contig', word, cube))
    for cube in random_cubes(NUM_RANDOM):
        candidates.append(('random_W', 'W', cube))

    out(f"Total candidate cubes: {len(candidates)}  "
        f"(80 contiguous + {NUM_RANDOM} random)")
    out("")

    # ---- Stage 1: experimental check for every cube ----------------------
    out("Stage 1: experimental cube-sum check (this is the cheap stage).")
    out("-" * 76)
    rng = random.Random(0xC0FFEE)
    exper_results = {}  # key = (cube_word, tuple(cube)) -> dict
    t_stage1 = time.time()
    for idx, (kind, word, cube) in enumerate(candidates):
        t0 = time.time()
        all_ok, sums = experimental_check(cube, word, num_keys=NUM_KEYS,
                                          num_rounds=NUM_ROUNDS, rng=rng)
        elapsed = time.time() - t0
        sw0, sx0, sy0, sz0 = sums[0]
        x_aabb = has_aabb(sx0)
        z_aabb = has_aabb(sz0)
        exper_results[(word, tuple(cube))] = {
            'kind': kind,
            'cube_word': word,
            'cube': cube,
            'all_aabb_xz': all_ok,
            'sums_key1': sums[0],
            'sums_all': sums,
            'x_aabb_key1': x_aabb,
            'z_aabb_key1': z_aabb,
            'time': elapsed,
        }
        label = f"{kind:>8} {word} {{{cube[0]:>2}..{cube[-1]:>2}}}" \
            if (cube == list(range(cube[0], cube[0] + CUBE_BITS))) \
            else f"{kind:>8} {word} {cube}"
        out(f"  [{idx+1:>3}/{len(candidates)}] {label:<35}  byte-wise={all_ok!s:<5}  "
            f"X1={sx0:08x} Z1={sz0:08x}  t={elapsed:.2f}s")
    out("-" * 76)
    out(f"Stage 1 elapsed: {time.time() - t_stage1:.1f}s")
    out("")

    # ---- Stage 2: MILP degree-bound confirmation -------------------------
    # To stay within the compute budget we focus the MILP runs:
    #   * every cube (any strategy) that passed the experimental check
    #     -> sample bits X^4_0 and Z^4_0 to *confirm* the 4-round bound
    #        governing the 5-round property
    #   * 2 representative "failed" contiguous cubes from word W
    #        (e.g. boundary positions k=0 and k=19) to record a deg_ub >= 13
    #        as a contrast
    milp_jobs = []
    seen = set()
    for entry in exper_results.values():
        key = (entry['cube_word'], tuple(entry['cube']))
        if key in seen:
            continue
        if not entry['all_aabb_xz']:
            continue  # confirmed-failure cubes get no MILP run
        seen.add(key)
        for tw in ('X', 'Z'):
            milp_jobs.append((
                list(entry['cube']),
                entry['cube_word'],
                tw,
                0,
                MILP_ROUNDS,
                MILP_TIME_LIMIT,
            ))

    # add a couple of "failed" cubes as MILP contrast (if any)
    contrasts_added = 0
    for entry in exper_results.values():
        if entry['kind'] != 'contig' or entry['all_aabb_xz']:
            continue
        if entry['cube_word'] != 'W':
            continue
        key = (entry['cube_word'], tuple(entry['cube']))
        if key in seen:
            continue
        seen.add(key)
        for tw in ('X', 'Z'):
            milp_jobs.append((list(entry['cube']), entry['cube_word'],
                              tw, 0, MILP_ROUNDS, MILP_TIME_LIMIT))
        contrasts_added += 1
        if contrasts_added >= 2:
            break

    out(f"Stage 2: MILP degree-bound confirmation (2 bits/cube, {MILP_TIME_LIMIT}s budget)")
    out(f"  Total MILP runs: {len(milp_jobs)}")
    out("-" * 76)
    milp_results = {}  # (cube_word, tuple(cube), target_word, bit) -> (ub, t)

    n_workers = max(1, min(4, mp.cpu_count() - 1))
    out(f"  Workers: {n_workers} parallel CBC processes")
    t_stage2 = time.time()
    with mp.Pool(n_workers) as pool:
        for cube_word, cube_tuple, tw, tb, ub, t in pool.imap_unordered(milp_worker, milp_jobs):
            milp_results[(cube_word, cube_tuple, tw, tb)] = (ub, t)
            label = f"{cube_word}={{{cube_tuple[0]:>2}..{cube_tuple[-1]:>2}}}" \
                if (list(cube_tuple) == list(range(cube_tuple[0], cube_tuple[0] + len(cube_tuple)))) \
                else f"{cube_word}={cube_tuple}"
            out(f"    {label:<25}  {tw}^4_{tb}  ub={ub:>3}  t={t:6.2f}s  "
                f"(elapsed: {time.time() - t_stage2:.1f}s)")
    out("-" * 76)
    out(f"Stage 2 elapsed: {time.time() - t_stage2:.1f}s")
    out("")

    # ---- Final classification table --------------------------------------
    out("=" * 76)
    out(" CLASSIFICATION TABLE")
    out("=" * 76)
    out(f"{'Strategy':<10} {'Word':<5} {'Cube':<24} {'BW-equal':<9} "
        f"{'ub(X4_0)':<9} {'ub(Z4_0)':<9}")
    out("-" * 76)

    summary_count = {'contig_W_bw': 0, 'contig_X_bw': 0, 'contig_Y_bw': 0, 'contig_Z_bw': 0,
                     'random_W_bw': 0,
                     'contig_W_total': 0, 'contig_X_total': 0,
                     'contig_Y_total': 0, 'contig_Z_total': 0,
                     'random_W_total': 0}

    rows_for_report = []
    for entry in exper_results.values():
        cube_word = entry['cube_word']
        cube = entry['cube']
        cube_tuple = tuple(cube)
        bw = entry['all_aabb_xz']

        ub_x = milp_results.get((cube_word, cube_tuple, 'X', 0), (None, None))[0]
        ub_z = milp_results.get((cube_word, cube_tuple, 'Z', 0), (None, None))[0]
        ub_x_str = "n/a" if ub_x is None else str(ub_x)
        ub_z_str = "n/a" if ub_z is None else str(ub_z)

        # cube display
        if list(cube) == list(range(cube[0], cube[0] + len(cube))):
            cube_str = f"{{{cube[0]},...,{cube[-1]}}}"
        else:
            cube_str = "{" + ",".join(str(b) for b in cube) + "}"

        strategy = entry['kind']
        out(f"{strategy:<10} {cube_word:<5} {cube_str:<24} "
            f"{('YES' if bw else 'no'):<9} {ub_x_str:<9} {ub_z_str:<9}")
        rows_for_report.append(dict(entry=entry, cube_str=cube_str,
                                     ub_x=ub_x, ub_z=ub_z))

        ctype = entry['kind']
        if ctype == 'contig':
            summary_count[f'contig_{cube_word}_total'] += 1
            if bw:
                summary_count[f'contig_{cube_word}_bw'] += 1
        else:
            summary_count['random_W_total'] += 1
            if bw:
                summary_count['random_W_bw'] += 1

    out("=" * 76)
    out(" SUMMARY")
    out("=" * 76)
    for word in 'WXYZ':
        bw = summary_count[f'contig_{word}_bw']
        tot = summary_count[f'contig_{word}_total']
        out(f"  contiguous in {word}: {bw}/{tot} cubes have byte-wise equal X & Z")
    out(f"  random in W       : {summary_count['random_W_bw']}/"
        f"{summary_count['random_W_total']} cubes have byte-wise equal X & Z")
    out("")

    # ---- Persist ----------------------------------------------------------
    OUT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nSaved to {OUT_PATH}")


if __name__ == "__main__":
    main()
