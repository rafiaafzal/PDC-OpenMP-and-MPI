/* CS1 parallel: 2-D heat diffusion with MPI.
 * usage: mpirun -np P ./parallel N T K MODE DIM
 *   N    grid size (N x N)
 *   T    time steps
 *   K    convergence check (MPI_Allreduce MAX) every K steps
 *   MODE 0 = blocking halo exchange (MPI_Sendrecv)
 *        1 = non-blocking (Isend/Irecv, interior computed while messages fly)
 *   DIM  1 = 1-D row-block decomposition  (dims = {P,1})
 *        2 = 2-D block decomposition       (dims from MPI_Dims_create)
 * Both decompositions use MPI_Cart_create; DIM=1 just forces dims={P,1}.
 * Global boundary cells are fixed (never updated); top row = 100, others 0.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static const double ALPHA = 0.2;

static int N, ly, lx, W, oy, ox;
static int istart, iend, jstart, jend;   /* cells this rank must update (local idx) */
static int ii0, ii1, jj0, jj1, interior_nonempty; /* halo-independent interior */
static double *cur, *nxt;

#define IDX(i, j) ((size_t)(i) * W + (j))

/* update rows i0..i1, cols j0..j1 (inclusive); returns max |change| if want_diff */
static double update(int i0, int i1, int j0, int j1, int want_diff) {
    double md = 0.0;
    for (int i = i0; i <= i1; i++) {
        const double *c = cur + IDX(i, 0);
        const double *up = c - W, *dn = c + W;
        double *o = nxt + IDX(i, 0);
        if (want_diff) {
            for (int j = j0; j <= j1; j++) {
                double v = c[j] + ALPHA * (up[j] + dn[j] + c[j - 1] + c[j + 1] - 4.0 * c[j]);
                o[j] = v;
                double d = fabs(v - c[j]);
                if (d > md) md = d;
            }
        } else {
            for (int j = j0; j <= j1; j++)
                o[j] = c[j] + ALPHA * (up[j] + dn[j] + c[j - 1] + c[j + 1] - 4.0 * c[j]);
        }
    }
    return md;
}

/* update the cells not covered by the interior rectangle */
static double update_ring(int want_diff) {
    double md = 0.0, d;
    if (istart > iend || jstart > jend) return 0.0;
    if (!interior_nonempty) return update(istart, iend, jstart, jend, want_diff);
    for (int i = istart; i <= iend; i++) {
        if (i < ii0 || i > ii1) {
            d = update(i, i, jstart, jend, want_diff);
            if (d > md) md = d;
        } else {
            if (jstart < jj0) { d = update(i, i, jstart, jj0 - 1, want_diff); if (d > md) md = d; }
            if (jj1 < jend)   { d = update(i, i, jj1 + 1, jend, want_diff);   if (d > md) md = d; }
        }
    }
    return md;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    N = argc > 1 ? atoi(argv[1]) : 1024;
    int T = argc > 2 ? atoi(argv[2]) : 500;
    int K = argc > 3 ? atoi(argv[3]) : 50;
    int mode = argc > 4 ? atoi(argv[4]) : 0;
    int dim = argc > 5 ? atoi(argv[5]) : 2;

    int dims[2] = {0, 0}, periods[2] = {0, 0}, coords[2];
    if (dim == 1) { dims[0] = p; dims[1] = 1; }
    else MPI_Dims_create(p, 2, dims);
    MPI_Comm cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart);
    MPI_Comm_rank(cart, &rank);
    MPI_Cart_coords(cart, rank, 2, coords);
    int up, down, left, right;
    MPI_Cart_shift(cart, 0, 1, &up, &down);
    MPI_Cart_shift(cart, 1, 1, &left, &right);

    /* block distribution with remainder */
    int by = N / dims[0], ry = N % dims[0];
    int bx = N / dims[1], rx = N % dims[1];
    ly = by + (coords[0] < ry);
    oy = coords[0] * by + (coords[0] < ry ? coords[0] : ry);
    lx = bx + (coords[1] < rx);
    ox = coords[1] * bx + (coords[1] < rx ? coords[1] : rx);
    if (ly < 1 || lx < 1) {
        if (rank == 0) fprintf(stderr, "grid too small for %d ranks\n", p);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    W = lx + 2;
    cur = calloc((size_t)(ly + 2) * W, sizeof(double));
    nxt = calloc((size_t)(ly + 2) * W, sizeof(double));
    if (!cur || !nxt) { fprintf(stderr, "alloc failed\n"); MPI_Abort(MPI_COMM_WORLD, 1); }

    if (oy == 0)                              /* hot top boundary */
        for (int j = 0; j < W; j++) cur[IDX(1, j)] = nxt[IDX(1, j)] = 100.0;

    istart = (oy == 0) ? 2 : 1;
    iend   = (oy + ly == N) ? ly - 1 : ly;
    jstart = (ox == 0) ? 2 : 1;
    jend   = (ox + lx == N) ? lx - 1 : lx;
    ii0 = istart > 2 ? istart : 2;       ii1 = iend < ly - 1 ? iend : ly - 1;
    jj0 = jstart > 2 ? jstart : 2;       jj1 = jend < lx - 1 ? jend : lx - 1;
    interior_nonempty = (ii0 <= ii1 && jj0 <= jj1);

    MPI_Datatype coltype;
    MPI_Type_vector(ly, 1, W, MPI_DOUBLE, &coltype);
    MPI_Type_commit(&coltype);

    double comm_time = 0.0, maxdiff = 0.0;
    MPI_Barrier(cart);
    double t0 = MPI_Wtime();

    for (int t = 1; t <= T; t++) {
        int check = (t % K == 0) || (t == T);
        double md = 0.0, d;
        double c0 = MPI_Wtime();
        if (mode == 0) {
            /* blocking: send first interior row up, receive bottom halo from below, etc. */
            MPI_Sendrecv(&cur[IDX(1, 1)],      lx, MPI_DOUBLE, up,    0,
                         &cur[IDX(ly + 1, 1)], lx, MPI_DOUBLE, down,  0, cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&cur[IDX(ly, 1)],     lx, MPI_DOUBLE, down,  1,
                         &cur[IDX(0, 1)],      lx, MPI_DOUBLE, up,    1, cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&cur[IDX(1, 1)],       1, coltype,    left,  2,
                         &cur[IDX(1, lx + 1)],  1, coltype,    right, 2, cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&cur[IDX(1, lx)],      1, coltype,    right, 3,
                         &cur[IDX(1, 0)],       1, coltype,    left,  3, cart, MPI_STATUS_IGNORE);
            comm_time += MPI_Wtime() - c0;
            md = update(istart, iend, jstart, jend, check);
        } else {
            MPI_Request rq[8];
            MPI_Irecv(&cur[IDX(0, 1)],       lx, MPI_DOUBLE, up,    0, cart, &rq[0]);
            MPI_Irecv(&cur[IDX(ly + 1, 1)],  lx, MPI_DOUBLE, down,  1, cart, &rq[1]);
            MPI_Irecv(&cur[IDX(1, 0)],        1, coltype,    left,  2, cart, &rq[2]);
            MPI_Irecv(&cur[IDX(1, lx + 1)],   1, coltype,    right, 3, cart, &rq[3]);
            MPI_Isend(&cur[IDX(1, 1)],       lx, MPI_DOUBLE, up,    1, cart, &rq[4]);
            MPI_Isend(&cur[IDX(ly, 1)],      lx, MPI_DOUBLE, down,  0, cart, &rq[5]);
            MPI_Isend(&cur[IDX(1, 1)],        1, coltype,    left,  3, cart, &rq[6]);
            MPI_Isend(&cur[IDX(1, lx)],       1, coltype,    right, 2, cart, &rq[7]);
            comm_time += MPI_Wtime() - c0;                    /* posting cost */
            if (interior_nonempty) md = update(ii0, ii1, jj0, jj1, check); /* overlap */
            c0 = MPI_Wtime();
            MPI_Waitall(8, rq, MPI_STATUSES_IGNORE);
            comm_time += MPI_Wtime() - c0;                    /* exposed (un-hidden) wait */
            d = update_ring(check);
            if (d > md) md = d;
        }
        if (check) {
            double c1 = MPI_Wtime();
            MPI_Allreduce(&md, &maxdiff, 1, MPI_DOUBLE, MPI_MAX, cart);
            comm_time += MPI_Wtime() - c1;
        }
        double *tmp = cur; cur = nxt; nxt = tmp;
    }
    double t1 = MPI_Wtime();

    double local = 0.0, sum;
    for (int i = 1; i <= ly; i++)
        for (int j = 1; j <= lx; j++) local += cur[IDX(i, j)];
    MPI_Reduce(&local, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, cart);
    double tloc = t1 - t0, tmax, cmax;
    MPI_Reduce(&tloc, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, cart);
    MPI_Reduce(&comm_time, &cmax, 1, MPI_DOUBLE, MPI_MAX, 0, cart);
    if (rank == 0)
        printf("RESULT cs1 variant=%s%dd N=%d T=%d p=%d K=%d time=%.6f comm=%.6f checksum=%.9e maxdiff=%.6e\n",
               mode ? "nonblocking_" : "blocking_", dim, N, T, p, K, tmax, cmax, sum, maxdiff);
    MPI_Type_free(&coltype);
    MPI_Finalize();
    return 0;
}
