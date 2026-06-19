/*
 * aradi.h -- Reference C implementation of the ARADI block cipher.
 *
 * ARADI is a 128-bit block cipher with a 256-bit key and 16 rounds, designed
 * by the NSA (Greene, Motley, Weeks, ePrint 2024/1240). The state is four
 * 32-bit words (W, X, Y, Z), corresponding to x[0..3] in the NSA C reference.
 *
 * This implementation accompanies a 2026 B.Sc. thesis on cube cryptanalysis
 * of ARADI, and passes the official NSA test vector:
 *
 *     key bytes  = 0x00, 0x01, ..., 0x1f
 *                  (LE-decoded: K[0]=0x03020100, ..., K[7]=0x1f1e1d1c)
 *     plaintext  = (0, 0, 0, 0)
 *     ciphertext = (0x3f09abf4, 0x00e3bd74, 0x03260def, 0xb7c53912)
 *
 * Public API:
 *     key_schedule(K, rk)      -- expand 8-word key into 17 round keys
 *     aradi_encrypt(pt, rk, ct, nr) -- encrypt one 128-bit block
 *     aradi_decrypt(ct, rk, pt, nr) -- decrypt one 128-bit block
 */

#ifndef ARADI_H
#define ARADI_H

#include <stdint.h>

#define ARADI_ROUNDS     16
#define ARADI_NUM_RKS    (ARADI_ROUNDS + 1)   /* 17: includes post-whitening */

/* 16-bit and 32-bit left rotations (n must be in 0..15 / 0..31). */
#define ROTL16(x, n) ((uint16_t)(((uint16_t)(x) << ((n) & 15)) | \
                                 ((uint16_t)(x) >> ((16 - ((n) & 15)) & 15))))
#define ROTL32(x, n) ((uint32_t)(((uint32_t)(x) << ((n) & 31)) | \
                                 ((uint32_t)(x) >> ((32 - ((n) & 31)) & 31))))
#define ROTR32(x, n) ((uint32_t)(((uint32_t)(x) >> ((n) & 31)) | \
                                 ((uint32_t)(x) << ((32 - ((n) & 31)) & 31))))

/*
 * Expand the 256-bit master key K[0..7] into 17 round keys rk[0..16][0..3].
 * rk[i][j] is k[i][j] in the NSA notation; rk[16] is the post-whitening key.
 */
void key_schedule(const uint32_t K[8], uint32_t rk[ARADI_NUM_RKS][4]);

/*
 * Encrypt a 128-bit plaintext block pt = (W, X, Y, Z) under the round keys
 * rk, producing ciphertext ct. num_rounds is normally ARADI_ROUNDS (=16);
 * smaller values are useful for round-reduced cryptanalysis.
 *
 * pt and ct may alias.
 */
void aradi_encrypt(const uint32_t pt[4],
                   const uint32_t rk[ARADI_NUM_RKS][4],
                   uint32_t ct[4],
                   int num_rounds);

/*
 * Decrypt ct -> pt. Uses the inverse S-box layer (each Toffoli is its own
 * inverse, applied in reverse order) and the involution property of the
 * linear layer.
 */
void aradi_decrypt(const uint32_t ct[4],
                   const uint32_t rk[ARADI_NUM_RKS][4],
                   uint32_t pt[4],
                   int num_rounds);

#endif /* ARADI_H */
