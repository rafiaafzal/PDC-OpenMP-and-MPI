/* xoshiro256** with jump() for provably non-overlapping parallel streams.
 * Stream t = seed state advanced by t jumps (each jump skips 2^128 outputs).
 * Also: Box-Muller normals, Black-Scholes closed form, timer. */
#ifndef RNG_H
#define RNG_H
#include <stdint.h>
#include <math.h>
#include <time.h>

typedef struct { uint64_t s[4]; } xo_t;

static inline double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t xo_next(xo_t *g) {
    uint64_t *s = g->s;
    uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return r;
}
static inline double xo_uniform(xo_t *g) { /* (0,1) */
    return ((xo_next(g) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}
static inline uint64_t splitmix(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline void xo_jump(xo_t *g) {
    static const uint64_t J[4] = {0x180ec6d33cfd0abaULL, 0xd5a61266f0c9392cULL,
                                  0xa9582618e03fc9aaULL, 0x39abdc4529b1661cULL};
    uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    for (int i = 0; i < 4; i++)
        for (int b = 0; b < 64; b++) {
            if (J[i] & (1ULL << b)) { s0 ^= g->s[0]; s1 ^= g->s[1]; s2 ^= g->s[2]; s3 ^= g->s[3]; }
            xo_next(g);
        }
    g->s[0] = s0; g->s[1] = s1; g->s[2] = s2; g->s[3] = s3;
}
/* stream number `stream` of the generator family identified by `seed` */
static inline xo_t xo_stream(uint64_t seed, int stream) {
    xo_t g; uint64_t x = seed;
    for (int i = 0; i < 4; i++) g.s[i] = splitmix(&x);
    for (int i = 0; i < stream; i++) xo_jump(&g);
    return g;
}
/* two independent N(0,1) draws via Box-Muller */
static inline void xo_normal2(xo_t *g, double *z0, double *z1) {
    double u1 = xo_uniform(g), u2 = xo_uniform(g);
    double r = sqrt(-2.0 * log(u1)), th = 6.283185307179586 * u2;
    *z0 = r * cos(th); *z1 = r * sin(th);
}
static inline double norm_cdf(double x) { return 0.5 * erfc(-x / 1.4142135623730951); }
static inline double black_scholes_call(double S0, double K, double r, double sig, double T) {
    double d1 = (log(S0 / K) + (r + 0.5 * sig * sig) * T) / (sig * sqrt(T));
    double d2 = d1 - sig * sqrt(T);
    return S0 * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2);
}
/* Parameters used everywhere */
#define S0_ 100.0
#define K_  100.0
#define R_  0.05
#define SIG_ 0.2
#define T_  1.0

/* simulate m paths (terminal price only); accumulates discounted payoff sum / sum of squares */
static inline void mc_batch(xo_t *g, long long m, double *sum, double *sumsq) {
    const double drift = (R_ - 0.5 * SIG_ * SIG_) * T_, vol = SIG_ * sqrt(T_), disc = exp(-R_ * T_);
    double s = 0, s2 = 0;
    for (long long i = 0; i < m; i += 2) {
        double z0, z1; xo_normal2(g, &z0, &z1);
        double p0 = disc * fmax(S0_ * exp(drift + vol * z0) - K_, 0.0);
        s += p0; s2 += p0 * p0;
        if (i + 1 < m) {
            double p1 = disc * fmax(S0_ * exp(drift + vol * z1) - K_, 0.0);
            s += p1; s2 += p1 * p1;
        }
    }
    *sum += s; *sumsq += s2;
}
static inline void mc_report(const char *variant, long long M, int p, double sum, double sumsq,
                             double secs, const char *extra) {
    double mean = sum / M, var = sumsq / M - mean * mean, se = sqrt(var / M);
    double bs = black_scholes_call(S0_, K_, R_, SIG_, T_);
    printf("RESULT cs3 variant=%s M=%lld p=%d time=%.6f price=%.8f stderr=%.8f bs=%.8f abs_err=%.8f err_in_se=%.2f Mpaths_s=%.2f %s\n",
           variant, M, p, secs, mean, se, bs, fabs(mean - bs), fabs(mean - bs) / se, M / secs / 1e6, extra);
}
#endif
