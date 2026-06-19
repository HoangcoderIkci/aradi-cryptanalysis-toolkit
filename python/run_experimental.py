"""
Experimental verification of the 5-round byte-wise equal property of ARADI.

For 5 random 256-bit master keys, sum (over the 2^13 cube I_W = {11,...,23})
the 5-round ARADI ciphertext word by word, and print the resulting four
32-bit words.  Per Kim et al. (and the original observation of Bellini et
al.), the X-word and the Z-word of the cube sum will be byte-wise equal:
both halves of each byte pattern AABB (e.g. 9d 9d c5 c5) -- though Kim
formalises this as: every 16-bit half equals its rotl16-by-8 (i.e.
S^8_16(u)=u and S^8_16(l)=l), which is exactly the AABB pattern.
"""

from __future__ import annotations

import os
import secrets
import time
from pathlib import Path

from aradi_ref import aradi_encrypt, key_schedule, MASK32

NUM_KEYS = 5
NUM_ROUNDS = 5
CUBE_W = list(range(11, 24))           # I_W = {11, 12, ..., 23}, |I| = 13
CUBE_SIZE_BITS = len(CUBE_W)
RESULTS_FILE = Path(__file__).with_name("experimental_5round.txt")


def has_aabb(word: int) -> bool:
    """Test whether a 32-bit word has the byte-wise equal property AABB."""
    # u = high 16 bits, l = low 16 bits.  Property: S^8_16(u)=u and S^8_16(l)=l
    u = (word >> 16) & 0xFFFF
    l = word & 0xFFFF
    def rotl8(x):
        return ((x << 8) | (x >> 8)) & 0xFFFF
    return rotl8(u) == u and rotl8(l) == l


def main() -> None:
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 70)
    out("Experimental verification of 5-round ARADI byte-wise equal property")
    out("=" * 70)
    out(f"Cube I_W = {{11, 12, ..., 23}}, |I| = {CUBE_SIZE_BITS}  -- cube size = 2^{CUBE_SIZE_BITS} = {1 << CUBE_SIZE_BITS}")
    out(f"Rounds: {NUM_ROUNDS}")
    out(f"Random master keys tested: {NUM_KEYS}")
    out("")
    out(f"{'#':>2} | {'cube sum (W, X, Y, Z)':^55} | X_AABB Z_AABB")
    out("-" * 90)

    all_aabb_xz = True
    for run in range(NUM_KEYS):
        # Random 256-bit master key (8 x 32-bit).
        key = [int.from_bytes(secrets.token_bytes(4), 'big') for _ in range(8)]
        rk = key_schedule(key)

        sum_w = sum_x = sum_y = sum_z = 0
        cube_count = 1 << CUBE_SIZE_BITS

        t0 = time.time()
        for c in range(cube_count):
            # Place the 13-bit cube value into bits 11..23 of W.
            w_in = 0
            for k, bit_pos in enumerate(CUBE_W):
                w_in |= ((c >> k) & 1) << bit_pos
            # X = Y = Z = 0 (non-cube indices fixed to zero).
            w_out, x_out, y_out, z_out = aradi_encrypt(
                w_in, 0, 0, 0, rk, num_rounds=NUM_ROUNDS)
            sum_w ^= w_out
            sum_x ^= x_out
            sum_y ^= y_out
            sum_z ^= z_out

        elapsed = time.time() - t0

        x_aabb = has_aabb(sum_x)
        z_aabb = has_aabb(sum_z)
        all_aabb_xz = all_aabb_xz and x_aabb and z_aabb

        out(f"{run+1:>2} | {sum_w:08x} {sum_x:08x} {sum_y:08x} {sum_z:08x}  ({elapsed:5.1f}s) | "
            f"  {'YES' if x_aabb else 'no '}    {'YES' if z_aabb else 'no '}")
        out(f"   |   X == Z ?  {'YES' if sum_x == sum_z else 'no '}")

    out("")
    out("=" * 70)
    out(f"X (cube sum) is byte-wise equal (AABB) for every test? "
        f"{'YES' if all_aabb_xz else 'NO'}")
    out(f"Z (cube sum) is byte-wise equal (AABB) for every test? "
        f"{'YES' if all_aabb_xz else 'NO'}")
    out("=" * 70)

    RESULTS_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nSaved to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
