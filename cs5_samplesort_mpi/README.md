# CS5 – Distributed sample sort (MPI, C++)

Files: `keygen.h`, `serial.cpp`, `parallel.cpp`, `pingpong.cpp`, `Makefile`, `run_bench.sh`. Build: `make` (`g++`, `mpicxx`).
Each rank generates its own N/p keys with a per-rank seed (nobody holds all data). Phases: local `std::sort` ->
regular sampling (OVERSAMPLE*(p-1) samples/rank, `MPI_Gather`, root picks p-1 splitters, `MPI_Bcast`) -> partition by `upper_bound` ->
`MPI_Alltoall` (counts) + `MPI_Alltoallv` (data) -> p-way merge (tree of `inplace_merge`).
Distributions: 0 uniform, 1 skewed (Zipf-like, u^4), 2 many duplicates (50% equal keys).

Run:
```
mpirun -np P ./parallel N dist OVERSAMPLE [verbose]    # mpirun -np 8 ./parallel 100000000 1 4 1
./serial N dist P                                      # serial std::sort of the same data a P-rank run generates
mpirun -np 2 ./pingpong > pp.csv && python3 ../tools/fit_alpha_beta.py pp.csv N p
```
Verification (printed `VERIFIED`): each rank sorted, max of lower ranks <= min of this rank, total count == N, sum of keys unchanged.
Serial and parallel print the same `checksum` for the same (N, dist, p).
Memory: N=1e8 keys = 400 MB per copy; the run needs about 3-4 copies in total -> lower `N=` if RAM is small.

Analysis hints: `max_over_mean` is the load balance (dist 2 with duplicates stays imbalanced: all equal keys go to one rank;
oversampling fixes skew, not duplicates); phase times `t_sort,t_sample,t_counts,t_alltoallv,t_merge` are max over ranks;
predict Alltoallv with T ~ (p-1)*alpha + beta*(N/p)*4*(p-1)/p using `fit_alpha_beta.py`; sub-linear speedup: sorting is
O((N/p) log(N/p)) per rank but all data (N*4 B) must still move once and be merged, and communication does not shrink like compute.
