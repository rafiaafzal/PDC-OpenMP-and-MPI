/* Key generation shared by serial.cpp and parallel.cpp.
 * dist 0 = uniform 32-bit; dist 1 = skewed (Zipf-like: u^4 concentrates keys near 0);
 * dist 2 = many duplicates (50% of keys are the single value 123456, rest uniform). */
#ifndef KEYGEN_H
#define KEYGEN_H
#include <cstdint>
#include <cmath>
#include <vector>
#include <chrono>
static inline uint64_t sm64(uint64_t &x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
/* generate `count` keys for virtual rank `rank` (seed depends only on rank) */
static inline void gen_keys(std::vector<uint32_t> &v, uint64_t count, int rank, int dist) {
    uint64_t st = 0xC0FFEE00ULL + 7919ULL * (uint64_t)rank;
    v.resize(count);
    for (uint64_t i = 0; i < count; i++) {
        uint64_t r = sm64(st);
        double u = (r >> 11) * (1.0 / 9007199254740992.0);
        if (dist == 0) v[i] = (uint32_t)(r >> 32);
        else if (dist == 1) { double w = u * u * u * u; v[i] = (uint32_t)(w * 4294967295.0); }
        else v[i] = (u < 0.5) ? 123456u : (uint32_t)(sm64(st) >> 32);
    }
}
static inline uint64_t my_count(uint64_t N, int p, int rank) { return N / p + ((uint64_t)rank < N % p ? 1 : 0); }
static inline double wall() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}
#endif
