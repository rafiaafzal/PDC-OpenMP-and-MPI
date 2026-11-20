/* CS3 OpenMP Monte Carlo. usage: ./parallel M threads seed MODE
 *   MODE 0 = correct: one independent xoshiro stream per thread (jump-ahead),
 *            partial sums combined in thread-index order -> deterministic for fixed (seed,p)
 *   MODE 1 = BAD: ONE shared generator protected by omp critical (serialisation demo;
 *            use a small M such as 1000000, it is very slow)
 *   MODE 2 = BAD: every thread uses the SAME seed (identical streams -> correlated samples,
 *            the reported standard error is far too small)
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "rng.h"

#define MAXT 1024
int main(int argc, char **argv) {
    long long M = argc > 1 ? atoll(argv[1]) : 100000000LL;
    int nt = argc > 2 ? atoi(argv[2]) : omp_get_max_threads();
    uint64_t seed = argc > 3 ? strtoull(argv[3], 0, 10) : 42;
    int mode = argc > 4 ? atoi(argv[4]) : 0;
    if (nt > MAXT) nt = MAXT;
    omp_set_num_threads(nt);
    static double ps[MAXT], ps2[MAXT];
    int actual = 1;
    xo_t shared_g = xo_stream(seed, 0);

    double t0 = now();
    #pragma omp parallel
    {
        int t = omp_get_thread_num(), P = omp_get_num_threads();
        #pragma omp single
        actual = P;
        long long m = M / P + (t < M % P ? 1 : 0);
        double s = 0, s2 = 0;
        if (mode == 0) {
            xo_t g = xo_stream(seed, t);
            mc_batch(&g, m, &s, &s2);
        } else if (mode == 2) {
            xo_t g = xo_stream(seed, 0);          /* same stream in every thread */
            mc_batch(&g, m, &s, &s2);
        } else {
            for (long long i = 0; i < m; i += 2) {
                #pragma omp critical(rng)
                {
                    mc_batch(&shared_g, 2, &s, &s2);
                }
            }
        }
        ps[t] = s; ps2[t] = s2;
    }
    double sum = 0, sumsq = 0;
    for (int t = 0; t < actual; t++) { sum += ps[t]; sumsq += ps2[t]; }   /* fixed order */
    double t1 = now();
    mc_report(mode == 0 ? "omp_independent" : mode == 1 ? "omp_shared_rng" : "omp_same_seed",
              M, actual, sum, sumsq, t1 - t0, "");
    return 0;
}
