/* CS2 serial baseline: k-mer counting.
 * usage: ./serial L k skew      (L = sequence length, default 100000000)
 */
#include <stdio.h>
#include <string.h>
#include "common.h"

int main(int argc, char **argv) {
    size_t L = argc > 1 ? strtoull(argv[1], 0, 10) : 100000000ULL;
    int k = argc > 2 ? atoi(argv[2]) : 8;
    int skew = argc > 3 ? atoi(argv[3]) : 0;
    uint8_t *s = gen_sequence(L, skew, 12345);
    size_t bins = (size_t)1 << (2 * k);
    uint32_t *h = (uint32_t *)calloc(bins, sizeof(uint32_t));
    uint32_t mask = (uint32_t)(bins - 1);

    double t0 = now();
    uint32_t code = 0;
    for (int i = 0; i < k - 1; i++) code = (code << 2) | s[i];
    for (size_t i = k - 1; i < L; i++) {
        code = ((code << 2) | s[i]) & mask;
        h[code]++;
    }
    double t1 = now();

    uint64_t total, cs = hist_checksum(h, bins, &total);
    double secs = t1 - t0;
    printf("RESULT cs2 variant=serial L=%zu k=%d skew=%d p=1 time=%.6f total=%llu expected=%zu checksum=%llu Mkmers_s=%.2f\n",
           L, k, skew, secs, (unsigned long long)total, L - k + 1, (unsigned long long)cs,
           (L - k + 1) / secs / 1e6);
    free(h); free(s);
    return 0;
}
