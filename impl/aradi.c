/*
 * aradi.c -- Reference C implementation of the ARADI block cipher.
 *
 * Mirrors the Python reference in milp/aradi_ref.py, which has been verified
 * against the official NSA test vector (ePrint 2024/1240). See aradi.h for
 * the public API.
 *
 * Notation:
 *   - State word indices 0..3 = (W, X, Y, Z).
 *   - The linear layer uses left rotations with parameters
 *         (11,8,14), (10,9,11), (9,4,14), (8,9,7),
 *     which are the rotl16 equivalents of the NSA C reference's ROR16
 *     parameters a={5,6,7,8}, b={8,7,12,7}, c={2,5,2,9}.
 */

#include "aradi.h"

/* ------------------------------------------------------------------------- */
/*  Linear layer L_i on a single 32-bit word.                                */
/* ------------------------------------------------------------------------- */

static const uint8_t LINEAR_A[4] = {11, 10,  9,  8};
static const uint8_t LINEAR_B[4] = { 8,  9,  4,  9};
static const uint8_t LINEAR_C[4] = {14, 11, 14,  7};

static uint32_t L_word(uint32_t word, unsigned a, unsigned b, unsigned c)
{
    uint16_t u = (uint16_t)(word >> 16);
    uint16_t l = (uint16_t)(word & 0xFFFFu);
    uint16_t u2 = (uint16_t)(u ^ ROTL16(u, a) ^ ROTL16(l, c));
    uint16_t l2 = (uint16_t)(l ^ ROTL16(l, a) ^ ROTL16(u, b));
    return ((uint32_t)u2 << 16) | (uint32_t)l2;
}

static void linear_layer(uint32_t *w, uint32_t *x, uint32_t *y, uint32_t *z,
                         int round_idx)
{
    unsigned i = (unsigned)(round_idx & 3);
    unsigned a = LINEAR_A[i], b = LINEAR_B[i], c = LINEAR_C[i];
    *w = L_word(*w, a, b, c);
    *x = L_word(*x, a, b, c);
    *y = L_word(*y, a, b, c);
    *z = L_word(*z, a, b, c);
}

/* ------------------------------------------------------------------------- */
/*  S-box layer (4 Toffoli gates, bit-sliced across 32 lanes).               */
/*                                                                           */
/*    X ^= W & Y                                                             */
/*    Z ^= X & Y                                                             */
/*    Y ^= W & Z                                                             */
/*    W ^= X & Z                                                             */
/*                                                                           */
/*  Each gate is its own inverse, so the inverse layer applies them in       */
/*  reverse order.                                                           */
/* ------------------------------------------------------------------------- */

static void sbox_layer(uint32_t *w, uint32_t *x, uint32_t *y, uint32_t *z)
{
    *x ^= (*w) & (*y);
    *z ^= (*x) & (*y);
    *y ^= (*w) & (*z);
    *w ^= (*x) & (*z);
}

static void sbox_layer_inv(uint32_t *w, uint32_t *x, uint32_t *y, uint32_t *z)
{
    *w ^= (*x) & (*z);
    *y ^= (*w) & (*z);
    *z ^= (*x) & (*y);
    *x ^= (*w) & (*y);
}

/* ------------------------------------------------------------------------- */
/*  Key schedule.                                                            */
/*                                                                           */
/*  In-place M-transform on the pair (K[x], K[x+1]):                         */
/*     t = K[x+1]                                                            */
/*     K[x+1] = ROR32(K[x], r1) ^ ROR32(K[x+1], r2) ^ K[x+1]                 */
/*     K[x]   = ROR32(K[x], r1) ^ t                                          */
/*                                                                           */
/*  Even pairs (k = 0, 2): r1 = 31, r2 = 29.                                 */
/*  Odd  pairs (k = 1, 3): r1 = 23, r2 =  4.                                 */
/* ------------------------------------------------------------------------- */

static void m_transform(uint32_t state[8], unsigned x_idx,
                        unsigned r1, unsigned r2)
{
    uint32_t a = state[x_idx];
    uint32_t b = state[x_idx + 1];
    state[x_idx + 1] = ROTR32(a, r1) ^ ROTR32(b, r2) ^ b;
    state[x_idx]     = ROTR32(a, r1) ^ b;
}

void key_schedule(const uint32_t K[8], uint32_t rk[ARADI_NUM_RKS][4])
{
    uint32_t state[8];
    int i;
    for (i = 0; i < 8; ++i) state[i] = K[i];

    for (i = 0; i < ARADI_NUM_RKS; ++i) {
        unsigned j = ((unsigned)i & 1u) << 2;   /* 0 or 4 */
        rk[i][0] = state[j + 0];
        rk[i][1] = state[j + 1];
        rk[i][2] = state[j + 2];
        rk[i][3] = state[j + 3];
        if (i == ARADI_ROUNDS) break;

        /* M-transforms on the four pairs. */
        for (unsigned k = 0; k < 4; ++k) {
            unsigned x = k + k;
            unsigned r1 = (k & 1u) ? 23u : 31u;
            unsigned r2 = (k & 1u) ?  4u : 29u;
            m_transform(state, x, r1, r2);
        }

        /* Counter XOR into K[7]. */
        state[7] ^= (uint32_t)i;

        /* Permutation: jj = 0 if i even, jj = 2 if i odd. */
        unsigned jj = ((unsigned)i & 1u) << 1;
        uint32_t t;
        t = state[1];        state[1] = state[2 + jj]; state[2 + jj] = t;
        t = state[5 - jj];   state[5 - jj] = state[6]; state[6] = t;
    }
}

/* ------------------------------------------------------------------------- */
/*  Encrypt / decrypt.                                                       */
/* ------------------------------------------------------------------------- */

void aradi_encrypt(const uint32_t pt[4],
                   const uint32_t rk[ARADI_NUM_RKS][4],
                   uint32_t ct[4],
                   int num_rounds)
{
    uint32_t w = pt[0], x = pt[1], y = pt[2], z = pt[3];
    int i;
    for (i = 0; i < num_rounds; ++i) {
        w ^= rk[i][0]; x ^= rk[i][1]; y ^= rk[i][2]; z ^= rk[i][3];
        sbox_layer(&w, &x, &y, &z);
        linear_layer(&w, &x, &y, &z, i);
    }
    /* Post-whitening key XOR (rk[num_rounds]). */
    ct[0] = w ^ rk[num_rounds][0];
    ct[1] = x ^ rk[num_rounds][1];
    ct[2] = y ^ rk[num_rounds][2];
    ct[3] = z ^ rk[num_rounds][3];
}

/*
 * The linear layer L_i is an involution over GF(2)^32 (rotation parameters
 * chosen so that L_i o L_i = identity in the rotl16 convention used here),
 * so we re-apply linear_layer with the same round index on the way back.
 * The S-box layer is inverted by running the four Toffoli gates in reverse.
 */
void aradi_decrypt(const uint32_t ct[4],
                   const uint32_t rk[ARADI_NUM_RKS][4],
                   uint32_t pt[4],
                   int num_rounds)
{
    uint32_t w = ct[0] ^ rk[num_rounds][0];
    uint32_t x = ct[1] ^ rk[num_rounds][1];
    uint32_t y = ct[2] ^ rk[num_rounds][2];
    uint32_t z = ct[3] ^ rk[num_rounds][3];

    int i;
    for (i = num_rounds - 1; i >= 0; --i) {
        linear_layer(&w, &x, &y, &z, i);    /* L_i is an involution */
        sbox_layer_inv(&w, &x, &y, &z);
        w ^= rk[i][0]; x ^= rk[i][1]; y ^= rk[i][2]; z ^= rk[i][3];
    }
    pt[0] = w; pt[1] = x; pt[2] = y; pt[3] = z;
}
