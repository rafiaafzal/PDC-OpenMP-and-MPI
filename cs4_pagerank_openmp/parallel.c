/* CS4 OpenMP PageRank (pull and push).
 * usage: ./parallel [-s scale] [-e ef] [-f snap.txt] [-m pull|push] [-S static|dynamic|guided]
 *                   [-c chunk] [-t threads] [-v 1] [-i maxiter] [-d degrees.csv]
 *   -m pull : each thread owns destination vertices, reads in-neighbours (no write conflicts)
 *   -m push : each thread pushes d*r_u/outdeg(u) to out-neighbours with #pragma omp atomic
 *   -S/-c   : loop schedule of the main vertex loop (schedule(runtime)); chunk 0 = default
 *   -v 1    : also run the serial version in-process, compare ranks (L1) and report speedup
 * Reports the load-imbalance ratio R = slowest-thread time / mean-thread time (main loop only).
 */
#include <omp.h>
#include "pr_common.h"

int main(int argc, char **argv) {
    int scale = atoi(argval(argc, argv, "-s", "20"));
    int ef = atoi(argval(argc, argv, "-e", "8"));
    const char *file = argval(argc, argv, "-f", NULL);
    const char *meth = argval(argc, argv, "-m", "pull");
    const char *sch = argval(argc, argv, "-S", "static");
    int chunk = atoi(argval(argc, argv, "-c", "0"));
    int nt = atoi(argval(argc, argv, "-t", "0"));
    int verify = atoi(argval(argc, argv, "-v", "1"));
    int maxit = atoi(argval(argc, argv, "-i", "100"));
    const char *deg = argval(argc, argv, "-d", NULL);
    int push = !strcmp(meth, "push");
    if (nt > 0) omp_set_num_threads(nt);
    omp_sched_t kind = !strcmp(sch, "dynamic") ? omp_sched_dynamic :
                       !strcmp(sch, "guided") ? omp_sched_guided : omp_sched_static;
    omp_set_schedule(kind, chunk);

    graph_t g;
    if (file) { if (!graph_load(&g, file)) { fprintf(stderr, "cannot open %s\n", file); return 1; } }
    else graph_rmat(&g, scale, ef, 2024);
    if (deg) dump_degrees(&g, deg);
    const uint32_t n = g.n;
    const double d = 0.85, eps = 1e-8;
    int P = omp_get_max_threads();

    /* analytic imbalance of schedule(static): edges (work) per thread for contiguous blocks */
    double static_imb;
    {
        double mx = 0, tot = 0;
        for (int t = 0; t < P; t++) {
            uint32_t b = n / P, r = n % P;
            uint32_t s = t * b + (t < (int)r ? t : r), e = s + b + (t < (int)r ? 1 : 0);
            double w = (double)(g.in_ptr[e] - g.in_ptr[s]);
            if (w > mx) mx = w; tot += w;
        }
        static_imb = mx / (tot / P);
    }

    double *r = (double *)malloc(n * sizeof(double));
    double *nr = (double *)malloc(n * sizeof(double));
    double *contrib = (double *)malloc(n * sizeof(double));
    double *tt = (double *)calloc(P, sizeof(double));
    #pragma omp parallel for schedule(static)
    for (long v = 0; v < (long)n; v++) r[v] = 1.0 / n;

    int it = 0;
    double t0 = omp_get_wtime();
    while (it < maxit) {
        double dang = 0;
        #pragma omp parallel for schedule(static) reduction(+:dang)
        for (long u = 0; u < (long)n; u++) {
            uint64_t od = g.out_ptr[u + 1] - g.out_ptr[u];
            if (od) contrib[u] = r[u] / (double)od; else { contrib[u] = 0; dang += r[u]; }
        }
        double base = (1.0 - d) / n + d * dang / n, diff = 0;
        if (!push) {
            #pragma omp parallel
            {
                int t = omp_get_thread_num();
                double ld = 0, s0 = omp_get_wtime();
                #pragma omp for schedule(runtime) nowait
                for (long v = 0; v < (long)n; v++) {
                    double s = 0;
                    for (uint64_t e = g.in_ptr[v]; e < g.in_ptr[v + 1]; e++) s += contrib[g.in_col[e]];
                    double x = base + d * s;
                    ld += fabs(x - r[v]);
                    nr[v] = x;
                }
                tt[t] += omp_get_wtime() - s0;
                #pragma omp atomic
                diff += ld;
            }
        } else {
            #pragma omp parallel for schedule(static)
            for (long v = 0; v < (long)n; v++) nr[v] = base;
            #pragma omp parallel
            {
                int t = omp_get_thread_num();
                double s0 = omp_get_wtime();
                #pragma omp for schedule(runtime) nowait
                for (long u = 0; u < (long)n; u++) {
                    if (contrib[u] == 0.0) continue;
                    double c = d * contrib[u];
                    for (uint64_t e = g.out_ptr[u]; e < g.out_ptr[u + 1]; e++) {
                        #pragma omp atomic
                        nr[g.out_col[e]] += c;
                    }
                }
                tt[t] += omp_get_wtime() - s0;
            }
            #pragma omp parallel for schedule(static) reduction(+:diff)
            for (long v = 0; v < (long)n; v++) diff += fabs(nr[v] - r[v]);
        }
        double *tmp = r; r = nr; nr = tmp;
        it++;
        if (diff < eps) break;
    }
    double t1 = omp_get_wtime();

    double mx = 0, mean = 0;
    for (int t = 0; t < P; t++) { if (tt[t] > mx) mx = tt[t]; mean += tt[t]; }
    mean /= P;
    double secs = t1 - t0;
    double bpe = push ? 20.0 : 12.0;     /* bytes moved per edge (see README) */
    double l1 = -1, tser = -1;
    if (verify) {
        double *rs = (double *)malloc(n * sizeof(double));
        double s0 = now();
        pagerank_serial(&g, d, eps, maxit, rs);
        tser = now() - s0;
        l1 = 0;
        for (uint32_t v = 0; v < n; v++) l1 += fabs(rs[v] - r[v]);
        free(rs);
    }
    printf("RESULT cs4 variant=%s sched=%s chunk=%d n=%u m=%llu p=%d iters=%d time=%.6f time_per_iter=%.6f "
           "imbalance_R=%.3f static_edge_imbalance=%.3f medges_s=%.2f est_GBs=%.2f l1_vs_serial=%.3e serial_time=%.6f %s\n",
           meth, sch, chunk, n, (unsigned long long)g.m, P, it, secs, secs / it, mean > 0 ? mx / mean : 1.0,
           static_imb, g.m * (double)it / secs / 1e6, bpe * g.m * (double)it / secs / 1e9, l1, tser,
           verify ? (l1 < 1e-6 ? "VERIFY_OK" : "VERIFY_FAIL") : "");
    return 0;
}
