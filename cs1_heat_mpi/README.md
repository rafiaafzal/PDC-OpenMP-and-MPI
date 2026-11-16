# CS1 – Heat diffusion (MPI)

Files: `serial.c`, `parallel.c`, `Makefile`, `run_bench.sh`.

Build: `make`   (needs `gcc`, `mpicc`)

Run:
```
./serial N T K                                  # e.g. ./serial 1024 500 50
mpirun -np P ./parallel N T K MODE DIM          # MODE 0=blocking Sendrecv, 1=non-blocking Isend/Irecv
                                                # DIM  1=1-D row blocks,    2=2-D blocks (MPI_Cart_create)
mpirun -np 4 ./parallel 1024 500 50 1 2         # example
```
Inputs used: N in {1024, 4096}, T = 500, K = 50 (K = 1 for the convergence-cost experiment), alpha = 0.2,
top boundary fixed at 100, other boundaries 0.  `./run_bench.sh` runs everything (strong scaling, weak scaling with
N = 1024*sqrt(p), K=1 vs K=50) and writes `results/`.

Verification: serial and parallel print `checksum` (sum of all cells) and `maxdiff`; they must agree to ~1e-9
relative (floating-point summation order differs). Check:
`./serial 512 100 10` vs `mpirun -np 4 ./parallel 512 100 10 0 2`.

Output fields: `time` = max over ranks of the time-step loop; `comm` = max over ranks of time spent in halo exchange
(blocking) or in posting + `Waitall` (non-blocking) + Allreduce.  For the non-blocking version `comm` is the
*exposed* communication, so overlap hidden = comm_blocking - comm_nonblocking.

Analysis hints (use YOUR numbers):
1. Surface-to-volume. Cells per rank = N^2/p in both. Halo values per rank per step: 1-D = 2N (independent of p);
   2-D (sqrt p x sqrt p) = 4N/sqrt(p) (shrinks with p). Halo/cells: 2p/N vs 4 sqrt(p)/N.
2. Overlap: compare `comm` and `time` of blocking vs non-blocking, same N/p/DIM. On one machine MPI uses shared
   memory and messages are small, so little to hide; many MPI implementations only progress messages inside MPI calls.
3. Strong vs weak: `results/cs1_strong_table.md`, `cs1_weak_table.md`. Fit Tp ~ Tcomp/p + Tcomm with Tcomp = T serial and
   Tcomm = measured `comm`.
4. K=1 vs K=50: an Allreduce is a global synchronization every step; compare `comm` and `time` in `results/cs1_K.log`.
