/* Shared helpers for CS2: timer, RNG, synthetic DNA generator. */
#ifndef COMMON_H
#define COMMON_H
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

static inline double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline uint64_t splitmix64(uint64_t *s) {
    uint64_t z = (*s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* Bases stored one per byte as 0..3 (A,C,G,T). skew=0: uniform.
 * skew=1: P(A)=0.40, P(T)=0.30, P(C)=0.15, P(G)=0.15 -> non-uniform histogram. */
static inline uint8_t *gen_sequence(size_t L, int skew, uint64_t seed) {
    uint8_t *s = (uint8_t *)malloc(L);
    uint64_t st = seed;
    for (size_t i = 0; i < L; i++) {
        uint64_t r = splitmix64(&st);
        if (!skew) s[i] = (uint8_t)(r & 3);
        else {
            unsigned u = (unsigned)((r >> 11) % 100);
            s[i] = u < 40 ? 0 : (u < 70 ? 3 : (u < 85 ? 1 : 2));
        }
    }
    return s;
}

/* order-dependent checksum of a histogram, used to compare serial vs parallel */
static inline uint64_t hist_checksum(const uint32_t *h, size_t bins, uint64_t *total) {
    uint64_t sum = 0, tot = 0;
    for (size_t i = 0; i < bins; i++) { tot += h[i]; sum += (uint64_t)h[i] * (i * 2654435761ULL + 1); }
    *total = tot;
    return sum;
}
#endif
