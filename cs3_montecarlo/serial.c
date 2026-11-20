/* CS3 serial Monte Carlo European call. usage: ./serial M seed */
#include <stdio.h>
#include <stdlib.h>
#include "rng.h"
int main(int argc, char **argv) {
    long long M = argc > 1 ? atoll(argv[1]) : 100000000LL;
    uint64_t seed = argc > 2 ? strtoull(argv[2], 0, 10) : 42;
    xo_t g = xo_stream(seed, 0);
    double sum = 0, sumsq = 0;
    double t0 = now();
    mc_batch(&g, M, &sum, &sumsq);
    double t1 = now();
    mc_report("serial", M, 1, sum, sumsq, t1 - t0, "");
    return 0;
}
