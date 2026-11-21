# Parallel Case Studies with OpenMP and MPI – Group XX
(Export to PDF as `GroupXX_report.pdf`, 10-14 pages. Replace every <<...>> with YOUR measured data, rewrite in your words.
Every graph: labelled axes, units, number of runs. Medians of >= 3 runs.)

## 1. Introduction
Team: <<names, roll numbers>>. Machine: <<CPU model, physical cores, logical cores, RAM, OS>>. Toolchain: <<gcc version, OpenMPI/MPICH version, flags -O3 -march=native>>.
Timing: wall-clock of the compute phase (data generation excluded), median of 3 runs, idle machine, serial baseline compiled with the same flags.
Oversubscription: <<state if p exceeded physical cores; speedup beyond physical cores is not meaningful>>.

## 2. CS1 – Heat diffusion (MPI)
**Problem / parallelizability:** each cell depends only on 4 neighbours of the previous step -> independent updates per step.
**Decomposition:** figure of 1-D row blocks and 2-D blocks with halos; MPI_Cart_create; blocking Sendrecv vs Isend/Irecv overlapping the interior update.
**Validation:** checksum serial = <<>>, parallel = <<>> (relative diff <<>>).
**Table/graph:** strong N=1024, N=4096; weak (N = 1024*sqrt(p)); <<results/cs1_*_table.md and PNG>>.
**Analysis**
1. Surface-to-volume: cells/rank = N^2/p. Halo values/rank/step: 1-D = 2N (constant in p); 2-D = 4N/sqrt(p). Halo/compute ratio: 2p/N vs 4sqrt(p)/N. 2-D scales better for large p; measured: <<>>.
2. Overlap: blocking comm <<>> s vs non-blocking exposed comm <<>> s -> hid <<>>%. If small: <<messages are N*8 B = ... small; shared-memory transport; MPI progress only inside MPI calls>>.
3. Strong vs weak: efficiency < 70% at p = <<>>. Model Tp ~ Tcomp/p + Tcomm with Tcomp = <<>>, Tcomm = <<>> -> predicted <<>> vs measured <<>>.
4. K=1 vs K=50: time <<>> vs <<>>; Allreduce is a latency-bound global synchronization each step; K=50 amortizes it 50x.

## 3. CS2 – k-mer counting (OpenMP)
**Decomposition:** k-mer start positions split in contiguous ranges; each thread primes the rolling window with k-1 bases, so boundary k-mers are counted exactly once.
**Analysis**
1. Version A total = <<>> vs expected L-k+1 = <<>>: non-atomic ++ = load, add, store; two threads load same value, one update is lost.
2. Atomic vs private: k=8 <<>> vs <<>>; k=12 <<>> vs <<>>. Private memory = threads*4^k*4 B = <<>> MB (k=8: 256 KB/thread, fits L2; k=12: 64 MB/thread, exceeds cache, merge cost grows). Crossover explanation: <<>>.
3. False sharing: packed (16 B per thread histogram, 4 threads per 64 B line) <<>> s vs padded <<>> s -> <<>>x.
4. Bytes/k-mer ~ 1 B (sequence) + 4 B RMW histogram; throughput <<>> Mkmers/s; scaling flattens at p = <<>> because <<memory bandwidth / cache / atomics>>.

## 4. CS3 – Monte Carlo option pricing (<<OpenMP>> + MPI variant)
Parameters S0=100,K=100,r=.05,sigma=.2,T=1,M=1e8; Black-Scholes = 10.4506. Model choice justification: single-node, embarrassingly parallel, only one reduction -> OpenMP simplest; MPI version shows same scaling.
RNG: xoshiro256** + jump() per thread; sums combined in thread order.
**Analysis**
1. Shared generator + critical: <<time>> (no speedup, lock contention). Same seed: stderr printed <<>> but true error <<>> (samples duplicated p times).
2. stderr ~ 1/sqrt(M): p=1,M=1e8 -> <<>>; p=8, M=8e8 same time -> <<>>. "Speedup" = accuracy per second: variance of estimate falls by p at equal time.
3. Amdahl: f = (1/S - 1/p)/(1 - 1/p) = <<>>; losses: reduction, thread start-up, turbo/thermal, hyper-threading.
4. Reproducibility: same (seed,p) identical <<price>>; different p differs by <<>> (< 1 stderr) because streams and summation order differ and FP addition is not associative.

## 5. CS4 – PageRank (OpenMP)
Graph: R-MAT scale <<20>>, <<>> vertices, <<>> edges (or SNAP <<name>>, cite). CSR. d=0.85, eps=1e-8, dangling mass redistributed. Verification: L1 vs serial = <<>>.
**Analysis**
1. Degree distribution (log-log figure): heavy tail, max degree <<>>. Static R = <<>> (edge-based analytic <<>>). Work per vertex ~ in-degree, so equal vertex counts != equal work.
2. Schedule table <<time vs schedule/chunk>>. Best: <<>>: small chunks balance but pay overhead, large chunks the reverse.
3. Pull vs push: push <<>>x slower. Pull = random reads, no write conflicts; push = random atomic read-modify-writes, hub cache lines contended.
4. Pull moves ~12 B/edge, push ~20 B/edge. Achieved <<>> Medges/s = <<>> GB/s vs machine bandwidth <<>> GB/s -> <<>>% (memory-bound if close).

## 6. CS5 – Sample sort (MPI)
Phases: local sort, regular sampling + Gather + splitters + Bcast, partition, Alltoall(counts)+Alltoallv(data), merge. Verified: sorted, boundaries, count, checksum.
**Analysis**
1. Final elements per rank (uniform / skewed / duplicates): max/mean = <<>> / <<>> / <<>>. Oversampling 1,4,16,64 -> <<>> (more samples = better splitters; duplicates need tie-breaking).
2. Phase breakdown (stacked bar): dominant at small p <<local sort>>, at large p <<Alltoallv/merge>>.
3. alpha = <<>> us, beta = <<>> ns/B from ping-pong. Predicted Alltoallv (p-1)alpha + beta*(N/p)*4*(p-1)/p = <<>> s vs measured <<>> s.
4. vs std::sort <<>> s: speedup <<>>. Sub-linear: compute O(N log N / p) but N*4 B must move and be merged; communication and merge do not shrink like the local sort.

## 7. Cross-case comparison
| Case | Model | Best speedup (p) | Dominant bottleneck |
|---|---|---|---|
| CS1 | MPI | <<>> | <<communication / sync>> |
| CS2 | OpenMP | <<>> | <<memory / atomics>> |
| CS3 | OpenMP (+MPI) | <<>> | <<compute>> |
| CS4 | OpenMP | <<>> | <<memory + imbalance>> |
| CS5 | MPI | <<>> | <<communication>> |
Best scaling: <<CS3 – independent work, one reduction>>.

## 8. Individual contributions
<<copy from CONTRIBUTIONS.md>>

## 9. References
Dataset: <<SNAP web-Google, Leskovec et al.>>; xoshiro256** (Blackman & Vigna); R-MAT (Chakrabarti et al.); OpenMP and MPI standards.
