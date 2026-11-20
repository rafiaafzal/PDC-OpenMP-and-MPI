# CS3 – Monte-Carlo European call (OpenMP main, MPI alternative)

Model choice: **OpenMP** for `parallel.c` (single node, shared-memory, almost no communication: one reduction at the end;
simplest model, lowest overhead). `parallel_mpi.c` is the MPI version for comparison (needs `mpicc`).
Parameters: S0=100, K=100, r=0.05, sigma=0.2, T=1 (Black-Scholes price 10.4506), M = 1e8 paths.
RNG: xoshiro256** with `jump()` -> thread/rank t uses stream t (non-overlapping 2^128 blocks); Box-Muller normals.
Partial sums are combined in thread/rank order -> result is bit-reproducible for fixed (seed, p).

Build: `make` (builds `serial`, `parallel`, `parallel_mpi`; if MPI is missing use `make serial parallel`).
Run:
```
./serial M seed                          # ./serial 100000000 42
./parallel M threads seed MODE           # MODE 0 independent streams | 1 shared generator + critical (BAD) | 2 same seed (BAD)
mpirun -np 8 ./parallel_mpi M seed
```
Output: `price`, `stderr` (standard error of the mean), `bs` (closed form), `abs_err`, `err_in_se` (|error|/stderr, should be ~<3).

Analysis hints:
1. RNG: MODE 1 serializes through a lock (time does not drop with p, use M=1e6); MODE 2 gives p identical streams, so the
   effective sample size is M/p but the printed stderr assumes M -> the true error is ~sqrt(p) times larger than stderr claims
   (compare `abs_err` and `err_in_se` over a few seeds).
2. Statistical efficiency: stderr ~ sigma_payoff/sqrt(M). At fixed time p workers do p*M paths -> stderr falls by sqrt(p).
   To compare, run `./parallel` with M=1e8 on 1 thread and with M=8e8 on 8 threads (same wall time) and compare `stderr`.
3. Amdahl: with S = T_serial/T_p, serial fraction f = (1/S - 1/p)/(1 - 1/p). Losses: thread start-up, final reduction,
   turbo/thermal drop with all cores busy, hyper-threads sharing a core (p above physical cores).
4. Reproducibility: same seed + same p -> identical price (run twice); different p -> different stream assignment and a
   different summation order (floating-point addition is not associative) -> price differs within ~1 stderr.
