"""
Reference Python implementation of ARADI block cipher.
Used for experimental verification of the byte-wise equal property.

ARADI specification (NSA, ePrint 2024/1240):
  - 128-bit block, 256-bit key, 16 rounds.
  - State: four 32-bit words (W, X, Y, Z); they map to x[0..3] in the C
    reference (i.e. W = x[0], X = x[1], Y = x[2], Z = x[3]).
  - Round = key_addition o S-box o Linear; final post-whitening key XOR
    after round 16.

This implementation passes the official NSA test vector
(eprint 2024/1240):

    key bytes = 0x00, 0x01, ..., 0x1f     (LE-decoded: K[0]=0x03020100, ...)
    plaintext = 0
    ciphertext words = (0x3f09abf4, 0x00e3bd74, 0x03260def, 0xb7c53912)

Public API (kept stable for multicube_attack.py and other experiments):
    LINEAR_PARAMS          -- list of 4 tuples (a, b, c), to be used with rotl16
    rotl16(x, n)           -- 16-bit left rotation
    L_word(word, a, b, c)  -- one-word linear transform (rotl16 convention)
    sbox_layer(w, x, y, z) -- 4-Toffoli S-box layer
    linear_layer(w, x, y, z, round_idx)
    key_schedule(K)        -- K = list of 8 uint32, returns 17 round keys
    aradi_encrypt(W, X, Y, Z, rk, num_rounds=16) -- returns 4-tuple ciphertext

NOTE on linear-layer convention. The NSA C reference uses ROR16 (right
rotation) with parameters a={5,6,7,8}, b={8,7,12,7}, c={2,5,2,9}. Since
ROR16(x, n) == rotl16(x, 16 - n), the mathematically identical L_i can be
expressed with rotl16 and parameters (16 - a, 16 - b, 16 - c). That is the
convention used here (and in multicube_attack.py), giving

    LINEAR_PARAMS = [(11, 8, 14), (10, 9, 11), (9, 4, 14), (8, 9, 7)]

The bug that previously broke the test vector was NOT in the linear layer
but in the key schedule (M-transform structure and pair ordering).
"""

MASK32 = 0xFFFFFFFF
MASK16 = 0xFFFF

# Linear-layer rotation parameters per (i mod 4): (a_i, b_i, c_i).
# These are the rotl16 equivalents of the NSA C reference's ROR16 parameters
# (a={5,6,7,8}, b={8,7,12,7}, c={2,5,2,9}). See module docstring.
LINEAR_PARAMS = [(11, 8, 14), (10, 9, 11), (9, 4, 14), (8, 9, 7)]


def rotl16(x: int, n: int) -> int:
    """16-bit left rotation."""
    n &= 15
    return ((x << n) | (x >> (16 - n))) & MASK16


def rotl32(x: int, n: int) -> int:
    """32-bit left rotation."""
    n &= 31
    return ((x << n) | (x >> (32 - n))) & MASK32


def rotr32(x: int, n: int) -> int:
    """32-bit right rotation."""
    n &= 31
    return ((x >> n) | (x << (32 - n))) & MASK32


def L_word(word: int, a: int, b: int, c: int) -> int:
    """
    Linear layer L_i applied to one 32-bit word, using rotl16 convention.

    Split 32-bit word into u (high 16 bits) and l (low 16 bits) and compute
        u' = u XOR rotl16(u, a) XOR rotl16(l, c)
        l' = l XOR rotl16(l, a) XOR rotl16(u, b)
    The result is (u' << 16) | l'.

    With parameters LINEAR_PARAMS this is equivalent to the NSA C reference
    L_i which uses ROR16 with (16-a, 16-b, 16-c).
    """
    u = (word >> 16) & MASK16
    l = word & MASK16
    u2 = (u ^ rotl16(u, a) ^ rotl16(l, c)) & MASK16
    l2 = (l ^ rotl16(l, a) ^ rotl16(u, b)) & MASK16
    return ((u2 << 16) | l2) & MASK32


def sbox_layer(w: int, x: int, y: int, z: int):
    """
    S-box layer (4 Toffoli gates), bit-sliced across 32 lanes.
    Order matches the NSA C reference (x[0..3] = (W, X, Y, Z)):
        x[1] ^= x[0] & x[2]    --  X ^= W & Y
        x[3] ^= x[1] & x[2]    --  Z ^= X & Y
        x[2] ^= x[0] & x[3]    --  Y ^= W & Z
        x[0] ^= x[1] & x[3]    --  W ^= X & Z
    """
    x ^= (w & y)
    z ^= (x & y)
    y ^= (w & z)
    w ^= (x & z)
    return w & MASK32, x & MASK32, y & MASK32, z & MASK32


def linear_layer(w: int, x: int, y: int, z: int, round_idx: int):
    """Apply L_i (round_idx mod 4 selects parameters)."""
    a, b, c = LINEAR_PARAMS[round_idx & 3]
    return (L_word(w, a, b, c), L_word(x, a, b, c),
            L_word(y, a, b, c), L_word(z, a, b, c))


def _m_transform(K: list, x_idx: int, r1: int, r2: int) -> None:
    """
    In-place M-transform on the pair (K[x_idx], K[x_idx+1]).

    Mirrors the NSA C reference exactly:
        t = K[x+1]
        K[x+1] = ROR32(K[x], r1) ^ ROR32(K[x+1], r2) ^ K[x+1]
        K[x]   = ROR32(K[x], r1) ^ t
    """
    a = K[x_idx]
    b = K[x_idx + 1]
    K[x_idx + 1] = (rotr32(a, r1) ^ rotr32(b, r2) ^ b) & MASK32
    K[x_idx]     = (rotr32(a, r1) ^ b) & MASK32


def key_schedule(K) -> list:
    """
    ARADI key schedule.

    K is a list of eight 32-bit words K[0..7] in the same order as the C
    reference: when the key bytes are 0x00,0x01,...,0x1f, the little-endian
    decoding gives K[0]=0x03020100, K[1]=0x07060504, ..., K[7]=0x1f1e1d1c.

    Returns a list of 17 round keys rk[0..16], each a 4-tuple (w, x, y, z).
    rk[i][0..3] corresponds to k[i][0..3] = (K[4j], K[4j+1], K[4j+2], K[4j+3])
    with j = i mod 2 at the start of iteration i.
    """
    state = [w & MASK32 for w in K]
    rk = []

    for i in range(17):
        j = (i & 1) << 2  # 0 or 4
        rk.append((state[j + 0], state[j + 1], state[j + 2], state[j + 3]))
        if i == 16:
            break
        # M-transforms on the four pairs (K[0],K[1]), (K[2],K[3]),
        # (K[4],K[5]), (K[6],K[7]).
        # k = 0, 2 (even): r1=31, r2=29  (ROR by 31,29 = ROL by 1,3)
        # k = 1, 3 (odd):  r1=23, r2= 4  (ROR by 23, 4 = ROL by 9,28)
        for k in range(4):
            x = k + k
            r1 = 23 if (k & 1) else 31
            r2 = 4  if (k & 1) else 29
            _m_transform(state, x, r1, r2)
        # Counter XOR
        state[7] = (state[7] ^ i) & MASK32
        # Permutation: jj = 0 if i even, jj = 2 if i odd
        jj = (i & 1) << 1
        state[1], state[2 + jj] = state[2 + jj], state[1]
        state[5 - jj], state[6] = state[6], state[5 - jj]
    return rk


def aradi_encrypt(W: int, X: int, Y: int, Z: int,
                  rk: list, num_rounds: int = 16):
    """
    Encrypt (W, X, Y, Z) with the given list of round keys.

    (W, X, Y, Z) correspond to x[0], x[1], x[2], x[3] in the C reference.
    For the standard cipher num_rounds=16 and len(rk) must be 17; the last
    entry rk[16] is the post-whitening key.
    """
    w, x, y, z = W & MASK32, X & MASK32, Y & MASK32, Z & MASK32
    for i in range(num_rounds):
        rk_w, rk_x, rk_y, rk_z = rk[i]
        w ^= rk_w
        x ^= rk_x
        y ^= rk_y
        z ^= rk_z
        w, x, y, z = sbox_layer(w, x, y, z)
        w, x, y, z = linear_layer(w, x, y, z, i)
    # Final post-whitening key XOR.
    rk_w, rk_x, rk_y, rk_z = rk[num_rounds]
    return ((w ^ rk_w) & MASK32, (x ^ rk_x) & MASK32,
            (y ^ rk_y) & MASK32, (z ^ rk_z) & MASK32)


if __name__ == "__main__":
    # Official NSA test vector (eprint 2024/1240)
    key = [0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c,
           0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c]
    rk = key_schedule(key)
    ct = aradi_encrypt(0, 0, 0, 0, rk, 16)
    expected = (0x3f09abf4, 0x00e3bd74, 0x03260def, 0xb7c53912)
    print("ARADI test vector (key 0x00..0x1f, plaintext 0):")
    print(f"  computed: {ct[0]:08x} {ct[1]:08x} {ct[2]:08x} {ct[3]:08x}")
    print(f"  expected: {expected[0]:08x} {expected[1]:08x} {expected[2]:08x} {expected[3]:08x}")
    print(f"  MATCH = {ct == expected}")
    assert ct == expected, "ARADI implementation FAILED the NSA test vector"
