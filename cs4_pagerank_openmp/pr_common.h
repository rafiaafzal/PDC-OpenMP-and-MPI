/* CS4 shared code: graph generation/loading, CSR construction, serial PageRank. */
#ifndef PR_COMMON_H
#define PR_COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

static inline double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static inline uint64_t sm64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

typedef struct {
    uint32_t n;
    uint64_t m;
    uint64_t *out_ptr; uint32_t *out_col;   /* CSR of out-edges (used by push)  */
    uint64_t *in_ptr;  uint32_t *in_col;    /* CSR of in-edges  (used by pull)  */
} graph_t;

static void build_csr(uint32_t n, uint64_t m, const uint32_t *key, const uint32_t *val,
                      uint64_t **ptr_out, uint32_t **col_out) {
    uint64_t *ptr = (uint64_t *)calloc((size_t)n + 1, sizeof(uint64_t));
    uint32_t *col = (uint32_t *)malloc(m * sizeof(uint32_t));
    for (uint64_t e = 0; e < m; e++) ptr[key[e] + 1]++;
    for (uint32_t v = 0; v < n; v++) ptr[v + 1] += ptr[v];
    uint64_t *pos = (uint64_t *)malloc(((size_t)n + 1) * sizeof(uint64_t));
    memcpy(pos, ptr, ((size_t)n + 1) * sizeof(uint64_t));
    for (uint64_t e = 0; e < m; e++) col[pos[key[e]]++] = val[e];
    free(pos);
    *ptr_out = ptr; *col_out = col;
}

static void graph_from_edges(graph_t *g, uint32_t n, uint64_t m, uint32_t *src, uint32_t *dst) {
    g->n = n; g->m = m;
    build_csr(n, m, src, dst, &g->out_ptr, &g->out_col);
    build_csr(n, m, dst, src, &g->in_ptr, &g->in_col);
}

/* R-MAT power-law generator: n = 2^scale, m = n*edge_factor, (a,b,c,d)=(.57,.19,.19,.05).
 * Ids are NOT permuted, so hubs sit at low ids (this makes static scheduling imbalance visible). */
static void graph_rmat(graph_t *g, int scale, int ef, uint64_t seed) {
    uint32_t n = (uint32_t)1 << scale;
    uint64_t m = (uint64_t)n * ef;
    uint32_t *src = (uint32_t *)malloc(m * 4), *dst = (uint32_t *)malloc(m * 4);
    uint64_t st = seed;
    const double a = 0.57, b = 0.19, c = 0.19;
    for (uint64_t e = 0; e < m; e++) {
        uint32_t u = 0, v = 0;
        for (int l = 0; l < scale; l++) {
            double r = (sm64(&st) >> 11) * (1.0 / 9007199254740992.0);
            u <<= 1; v <<= 1;
            if (r < a) { }
            else if (r < a + b) v |= 1;
            else if (r < a + b + c) u |= 1;
            else { u |= 1; v |= 1; }
        }
        src[e] = u; dst[e] = v;
    }
    graph_from_edges(g, n, m, src, dst);
    free(src); free(dst);
}

/* SNAP edge list ("u v" per line, lines starting with # are comments). n = max id + 1. */
static int graph_load(graph_t *g, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    size_t cap = 1 << 20, m = 0;
    uint32_t *src = (uint32_t *)malloc(cap * 4), *dst = (uint32_t *)malloc(cap * 4);
    char line[256]; uint32_t maxid = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 35) continue;
        unsigned u, v;
        if (sscanf(line, "%u %u", &u, &v) != 2) continue;
        if (m == cap) { cap *= 2; src = (uint32_t *)realloc(src, cap * 4); dst = (uint32_t *)realloc(dst, cap * 4); }
        src[m] = u; dst[m] = v; m++;
        if (u > maxid) maxid = u;
        if (v > maxid) maxid = v;
    }
    fclose(f);
    graph_from_edges(g, maxid + 1, m, src, dst);
    free(src); free(dst);
    return 1;
}

/* degree histogram CSV: degree,count_in,count_out  (for the log-log plot) */
static void dump_degrees(const graph_t *g, const char *path) {
    uint32_t maxd = 0;
    for (uint32_t v = 0; v < g->n; v++) {
        uint32_t di = (uint32_t)(g->in_ptr[v + 1] - g->in_ptr[v]);
        uint32_t dout = (uint32_t)(g->out_ptr[v + 1] - g->out_ptr[v]);
        if (di > maxd) maxd = di;
        if (dout > maxd) maxd = dout;
    }
    uint64_t *ci = (uint64_t *)calloc(maxd + 1, 8), *co = (uint64_t *)calloc(maxd + 1, 8);
    for (uint32_t v = 0; v < g->n; v++) {
        ci[g->in_ptr[v + 1] - g->in_ptr[v]]++;
        co[g->out_ptr[v + 1] - g->out_ptr[v]]++;
    }
    FILE *f = fopen(path, "w");
    if (f) {
        fprintf(f, "degree,count_in,count_out\n");
        for (uint32_t d = 0; d <= maxd; d++)
            if (ci[d] || co[d]) fprintf(f, "%u,%llu,%llu\n", d, (unsigned long long)ci[d], (unsigned long long)co[d]);
        fclose(f);
    }
    free(ci); free(co);
}

/* Serial pull-style PageRank. Dangling nodes (outdeg 0) spread their rank uniformly:
 *   r_v = (1-d)/n + d * ( sum_{u->v} r_u/outdeg(u) + dangling_mass/n )
 * Stops when ||r_new - r_old||_1 < eps or after maxit iterations. Returns iterations. */
static int pagerank_serial(const graph_t *g, double d, double eps, int maxit, double *r) {
    uint32_t n = g->n;
    double *nr = (double *)malloc(n * sizeof(double));
    double *contrib = (double *)malloc(n * sizeof(double));
    for (uint32_t v = 0; v < n; v++) r[v] = 1.0 / n;
    int it = 0;
    for (; it < maxit;) {
        double dang = 0;
        for (uint32_t u = 0; u < n; u++) {
            uint64_t od = g->out_ptr[u + 1] - g->out_ptr[u];
            if (od) contrib[u] = r[u] / (double)od; else { contrib[u] = 0; dang += r[u]; }
        }
        double base = (1.0 - d) / n + d * dang / n, diff = 0;
        for (uint32_t v = 0; v < n; v++) {
            double s = 0;
            for (uint64_t e = g->in_ptr[v]; e < g->in_ptr[v + 1]; e++) s += contrib[g->in_col[e]];
            double x = base + d * s;
            diff += fabs(x - r[v]);
            nr[v] = x;
        }
        memcpy(r, nr, n * sizeof(double));
        it++;
        if (diff < eps) break;
    }
    free(nr); free(contrib);
    return it;
}

/* tiny "-key value" argument parser */
static const char *argval(int argc, char **argv, const char *key, const char *def) {
    for (int i = 1; i + 1 < argc; i++) if (!strcmp(argv[i], key)) return argv[i + 1];
    return def;
}
#endif
