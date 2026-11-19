// CS5 serial baseline: std::sort of the whole data set.
// usage: ./serial N dist P
//   P = number of "virtual ranks": the data is generated exactly as a P-rank run of
//   parallel would generate it (same per-rank seeds), so checksums can be compared.
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "keygen.h"
int main(int argc, char **argv) {
    uint64_t N = argc > 1 ? strtoull(argv[1], 0, 10) : 100000000ULL;
    int dist = argc > 2 ? atoi(argv[2]) : 0;
    int P = argc > 3 ? atoi(argv[3]) : 1;
    std::vector<uint32_t> all; all.reserve(N);
    std::vector<uint32_t> part;
    for (int r = 0; r < P; r++) {
        gen_keys(part, my_count(N, P, r), r, dist);
        all.insert(all.end(), part.begin(), part.end());
    }
    std::vector<uint32_t>().swap(part);
    uint64_t sum = 0; for (uint32_t x : all) sum += x;
    double t0 = wall();
    std::sort(all.begin(), all.end());
    double t1 = wall();
    printf("RESULT cs5 variant=serial N=%llu dist=%d p=1 time=%.6f sorted=%d checksum=%llu\n",
           (unsigned long long)N, dist, t1 - t0, (int)std::is_sorted(all.begin(), all.end()),
           (unsigned long long)sum);
    return 0;
}
