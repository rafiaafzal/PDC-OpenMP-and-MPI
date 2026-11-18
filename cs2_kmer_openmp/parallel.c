/* CS2 parallel k-mer counting with OpenMP.
 * usage: ./parallel L k skew MODE threads
 *   MODE:
 *     A  data race on purpose (plain shared ++)          -> WRONG result
 *     B  #pragma omp atomic on the shared histogram
 *     C  per-thread private histograms (cache-line aligned), merged at the end
 *     F  per-thread histograms packed WITHOUT padding (false sharing demo, use k=1 or 2)
 *     P  same as F but each private histogram padded to a multiple of 64 bytes
 *
 * Chunking: the k-mer START positions 0..L-k are split into contiguous ranges, one
 * per thread. A thread owning start positions [s,e) reads bases s..e+k-2, so a k-mer
 * that straddles two chunks of *bases* is still counted exactly once, by the thread
 * that owns its first base. Each thread primes its rolling code with k-1 bases first.
 */
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "common.h"

int main(int argc, char **argv) {
    size_t L = argc > 1 ? strtoull(argv[1], 0, 10) : 100000000ULL;
    int k = argc > 2 ? atoi(argv[2]) : 8;
    int skew = argc > 3 ? atoi(argv[3]) : 0;
    char mode = argc > 4 ? argv[4][0] : (char)67;
    int nt = argc > 5 ? atoi(argv[5]) : omp_get_max_threads();
    omp_set_num_threads(nt);

    uint8_t *s = gen_sequence(L, skew, 12345);
    size_t bins = (size_t)1 << (2 * k);
    uint32_t mask = (uint32_t)(bins - 1);
    size_t nk = L - k + 1;                       /* number of k-mers */

    uint32_t *shared = NULL, *priv = NULL;
    size_t stride = 0;
    if (mode == 65 || mode == 66) {              /* A or B */
        shared = (uint32_t *)calloc(bins, sizeof(uint32_t));
    } else {
        size_t padded = (bins + 15) & ~(size_t)15;   /* 16 x 4B = 64B cache line */
        stride = (mode == 70) ? bins : padded;       /* F: packed, else padded */
        priv = (uint32_t *)aligned_alloc(64, ((stride * nt * sizeof(uint32_t)) + 63) & ~(size_t)63);
        memset(priv, 0, stride * nt * sizeof(uint32_t));
    }

    int actual = 1;
    double t0 = now();
    #pragma omp parallel
    {
        int t = omp_get_thread_num();
        #pragma omp single
        actual = omp_get_num_threads();
        int P = omp_get_num_threads();
        size_t b = nk / P, r = nk % P;
        size_t st = t * b + (t < (int)r ? t : r);
        size_t en = st + b + (t < (int)r ? 1 : 0);
        uint32_t code = 0;
        for (int i = 0; i < k - 1; i++) code = (code << 2) | s[st + i];
        if (mode == 65) {
            for (size_t i = st; i < en; i++) {
                code = ((code << 2) | s[i + k - 1]) & mask;
                shared[code]++;                                  /* RACE */
            }
        } else if (mode == 66) {
            for (size_t i = st; i < en; i++) {
                code = ((code << 2) | s[i + k - 1]) & mask;
                #pragma omp atomic
                shared[code]++;
            }
        } else {
            uint32_t *h = priv + (size_t)t * stride;
            for (size_t i = st; i < en; i++) {
                code = ((code << 2) | s[i + k - 1]) & mask;
                h[code]++;
            }
        }
    }
    uint32_t *hist = shared;
    if (priv) {                                   /* merge private histograms */
        hist = (uint32_t *)calloc(bins, sizeof(uint32_t));
        #pragma omp parallel for schedule(static)
        for (long long b = 0; b < (long long)bins; b++) {
            uint32_t sum = 0;
            for (int t = 0; t < actual; t++) sum += priv[(size_t)t * stride + b];
            hist[b] = sum;
        }
    }
    double t1 = now();

    uint64_t total, cs = hist_checksum(hist, bins, &total);
    double secs = t1 - t0;
    double extra_mb = priv ? (double)stride * actual * 4 / 1048576.0 : 0.0;
    printf("RESULT cs2 variant=%c L=%zu k=%d skew=%d p=%d time=%.6f total=%llu expected=%zu checksum=%llu Mkmers_s=%.2f extraMB=%.2f %s\n",
           mode, L, k, skew, actual, secs, (unsigned long long)total, nk, (unsigned long long)cs,
           nk / secs / 1e6, extra_mb, total == nk ? "COUNT_OK" : "COUNT_WRONG");
    return 0;
}
