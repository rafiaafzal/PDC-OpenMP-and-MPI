// Ping-pong between ranks 0 and 1 to estimate Tcomm = alpha + beta*m.
// usage: mpirun -np 2 ./pingpong      -> CSV: bytes,one_way_seconds
// Fit alpha/beta with:  python3 ../tools/fit_alpha_beta.py pingpong.csv
#include <mpi.h>
#include <cstdio>
#include <vector>
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int r, p; MPI_Comm_rank(MPI_COMM_WORLD, &r); MPI_Comm_size(MPI_COMM_WORLD, &p);
    if (p < 2) { if (!r) fprintf(stderr, "run with at least 2 ranks\n"); MPI_Finalize(); return 1; }
    if (r == 0) printf("bytes,one_way_seconds\n");
    for (long m = 1; m <= (64L << 20); m *= 4) {
        std::vector<char> buf(m, 1);
        int reps = m < 4096 ? 2000 : (m < (1 << 20) ? 200 : 20);
        MPI_Barrier(MPI_COMM_WORLD);
        double t0 = MPI_Wtime();
        for (int i = 0; i < reps; i++) {
            if (r == 0) { MPI_Send(buf.data(), (int)m, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                          MPI_Recv(buf.data(), (int)m, MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); }
            else if (r == 1) { MPI_Recv(buf.data(), (int)m, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                               MPI_Send(buf.data(), (int)m, MPI_CHAR, 0, 0, MPI_COMM_WORLD); }
        }
        double t = (MPI_Wtime() - t0) / reps / 2.0;
        if (r == 0) printf("%ld,%.9e\n", m, t);
    }
    MPI_Finalize();
    return 0;
}
