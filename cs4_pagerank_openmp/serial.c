/* CS4 serial PageRank baseline.
 * usage: ./serial [-s scale] [-e edgefactor] [-f snap_edge_list.txt] [-i maxiter] [-d degrees.csv]
 * Default graph: R-MAT scale 20 (1,048,576 vertices), edge factor 8 (~8.4M edges). */
#include "pr_common.h"
int main(int argc, char **argv) {
    int scale = atoi(argval(argc, argv, "-s", "20"));
    int ef = atoi(argval(argc, argv, "-e", "8"));
    const char *file = argval(argc, argv, "-f", NULL);
    int maxit = atoi(argval(argc, argv, "-i", "100"));
    const char *deg = argval(argc, argv, "-d", NULL);
    graph_t g;
    if (file) { if (!graph_load(&g, file)) { fprintf(stderr, "cannot open %s\n", file); return 1; } }
    else graph_rmat(&g, scale, ef, 2024);
    if (deg) dump_degrees(&g, deg);
    double *r = (double *)malloc(g.n * sizeof(double));
    double t0 = now();
    int it = pagerank_serial(&g, 0.85, 1e-8, maxit, r);
    double t1 = now();
    double sum = 0; for (uint32_t v = 0; v < g.n; v++) sum += r[v];
    printf("RESULT cs4 variant=serial n=%u m=%llu p=1 iters=%d time=%.6f time_per_iter=%.6f medges_s=%.2f rank_sum=%.9f\n",
           g.n, (unsigned long long)g.m, it, t1 - t0, (t1 - t0) / it, g.m * (double)it / (t1 - t0) / 1e6, sum);
    return 0;
}
