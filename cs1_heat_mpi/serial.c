/* CS1 serial baseline: 2-D heat diffusion, 5-point stencil.
 * usage: ./serial N T K
 *   N = grid is N x N (boundary included), T = time steps,
 *   K = compute max-change every K steps (kept for parity with the MPI code).
 * Boundary: top row fixed at 100, all other boundaries fixed at 0.
 * Interior initial value 0. alpha = 0.2 (stable for alpha <= 0.25).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(int argc, char **argv) {
    int N = argc > 1 ? atoi(argv[1]) : 1024;
    int T = argc > 2 ? atoi(argv[2]) : 500;
    int K = argc > 3 ? atoi(argv[3]) : 50;
    const double alpha = 0.2;
    double *cur = calloc((size_t)N * N, sizeof(double));
    double *nxt = calloc((size_t)N * N, sizeof(double));
    if (!cur || !nxt) { fprintf(stderr, "alloc failed\n"); return 1; }
    for (int j = 0; j < N; j++) cur[j] = nxt[j] = 100.0; /* hot top edge */

    double maxdiff = 0.0;
    double t0 = now();
    for (int t = 1; t <= T; t++) {
        int check = (t % K == 0) || (t == T);
        double md = 0.0;
        for (int i = 1; i < N - 1; i++) {
            const double *c = cur + (size_t)i * N;
            const double *up = c - N, *dn = c + N;
            double *o = nxt + (size_t)i * N;
            for (int j = 1; j < N - 1; j++) {
                double v = c[j] + alpha * (up[j] + dn[j] + c[j - 1] + c[j + 1] - 4.0 * c[j]);
                o[j] = v;
                if (check) { double d = fabs(v - c[j]); if (d > md) md = d; }
            }
        }
        if (check) maxdiff = md;
        double *tmp = cur; cur = nxt; nxt = tmp;
    }
    double t1 = now();
    double sum = 0.0;
    for (size_t i = 0; i < (size_t)N * N; i++) sum += cur[i];
    printf("RESULT cs1 variant=serial N=%d T=%d p=1 K=%d time=%.6f comm=0 checksum=%.9e maxdiff=%.6e\n",
           N, T, K, t1 - t0, sum, maxdiff);
    free(cur); free(nxt);
    return 0;
}
