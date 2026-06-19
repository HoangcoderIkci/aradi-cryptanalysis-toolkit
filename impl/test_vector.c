/*
 * test_vector.c -- Verify the C ARADI implementation against the official
 * NSA test vector (ePrint 2024/1240).
 *
 *     key bytes  = 0x00, 0x01, ..., 0x1f       (little-endian into K[0..7])
 *     plaintext  = (0, 0, 0, 0)
 *     ciphertext = (0x3f09abf4, 0x00e3bd74, 0x03260def, 0xb7c53912)
 *
 * Exit code: 0 on PASS, 1 on FAIL.
 */

#include <stdint.h>
#include <stdio.h>

#include "aradi.h"

int main(void)
{
    /* Master key, little-endian decoding of bytes 0x00..0x1f. */
    const uint32_t key[8] = {
        0x03020100u, 0x07060504u, 0x0b0a0908u, 0x0f0e0d0cu,
        0x13121110u, 0x17161514u, 0x1b1a1918u, 0x1f1e1d1cu
    };
    const uint32_t pt[4] = {0u, 0u, 0u, 0u};
    const uint32_t expected[4] = {
        0x3f09abf4u, 0x00e3bd74u, 0x03260defu, 0xb7c53912u
    };

    uint32_t rk[ARADI_NUM_RKS][4];
    uint32_t ct[4];
    key_schedule(key, rk);
    aradi_encrypt(pt, rk, ct, ARADI_ROUNDS);

    printf("ARADI NSA test vector (key 0x00..0x1f, plaintext = 0):\n");
    printf("  computed: %08x %08x %08x %08x\n",
           ct[0], ct[1], ct[2], ct[3]);
    printf("  expected: %08x %08x %08x %08x\n",
           expected[0], expected[1], expected[2], expected[3]);

    int ok = 1;
    for (int i = 0; i < 4; ++i) {
        if (ct[i] != expected[i]) ok = 0;
    }

    /* Round-trip sanity check: decrypt(encrypt(pt)) must equal pt. */
    uint32_t recovered[4];
    aradi_decrypt(ct, rk, recovered, ARADI_ROUNDS);
    int rt_ok = (recovered[0] == pt[0]) && (recovered[1] == pt[1]) &&
                (recovered[2] == pt[2]) && (recovered[3] == pt[3]);

    if (ok && rt_ok) {
        printf("PASS  (encryption matches and decrypt round-trip ok)\n");
        return 0;
    }
    if (!ok)    printf("FAIL: ciphertext does not match NSA test vector\n");
    if (!rt_ok) printf("FAIL: decrypt(encrypt(pt)) != pt\n");
    return 1;
}
