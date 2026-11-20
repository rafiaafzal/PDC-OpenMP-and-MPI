/* CS3 MPI variant (distributed-memory alternative to parallel.c).
 * usage: mpirun -np P ./parallel_mpi M seed
 * Rank r uses xoshiro stream r. Partial sums are gathered and added in rank order,
 * so results are reproducible for a fixed (seed, P). */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "rng.h"
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int r, p; MPI_Comm_rank(MPI_COMM_WORLD, &r); MPI_Comm_size(MPI_COMM_WORLD, &p);
    long long M = argc > 1 ? atoll(argv[1]) : 100000000LL;
    uint64_t seed = argc > 2 ? strtoull(argv[2], 0, 10) : 42;
    long long m = M / p + (r < M % p ? 1 : 0);
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    xo_t g = xo_stream(seed, r);
    double part[2] = {0, 0};
    mc_batch(&g, m, &part[0], &part[1]);
    double *all = r == 0 ? (double *)malloc(2 * p * sizeof(double)) : NULL;
    MPI_Gather(part, 2, MPI_DOUBLE, all, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double tl = t1 - t0, tm;
    MPI_Reduce(&tl, &tm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (r == 0) {
        double sum = 0, sumsq = 0;
        for (int i = 0; i < p; i++) { sum += all[2 * i]; sumsq += all[2 * i + 1]; }
        mc_report("mpi_independent", M, p, sum, sumsq, tm, "");
        free(all);
    }
    MPI_Finalize();
    return 0;
}
