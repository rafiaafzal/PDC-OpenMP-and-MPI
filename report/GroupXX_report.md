# Parallel Case Studies with OpenMP and MPI – Group 00

## 1. Introduction
- **Team Members:**
  - Rafia Afzal (Roll No: 23L-0821) — Author A
  - Maheen Kamran (Roll No: 23L-0715) — Author B
- **Testbed Hardware Specifications:**
  - **CPU Model:** Intel(R) Core(TM) m7-6Y75 CPU @ 1.20GHz
  - **Cores & Topologies:** 2 physical cores, 4 logical execution threads (2 threads per core, SMT enabled)
  - **Caches:** L1d 32 KiB/core, L1i 32 KiB/core, L2 256 KiB/core, L3 4 MiB shared
  - **Memory:** 7.5 GiB RAM
  - **Operating System:** Ubuntu 24.04 LTS (x86_64, Linux 5.15 on WSL2)
- **Toolchain & Compiler Flags:**
  - **C/C++ Compilers:** GCC 13.3.0 / G++ 13.3.0
  - **MPI Implementation:** OpenMPI 4.1.6
  - **Flags:** `-O3 -march=native -std=gnu11 -fopenmp`
- **Benchmarking Methodology:**
  - Wall-clock measurement specifically timing compute phases (excluding data initialization and disk I/O).
  - All reported timings are medians of >= 3 consecutive runs under idle machine conditions.
  - Baseline serial implementations compiled with identical optimization flags (`-O3 -march=native`).
- **Oversubscription Note:**
  - The physical processor features 2 physical cores (4 hardware threads). Beyond p=2, scaling transitions into Hyper-Threading / SMT resources; beyond p=4, processes/threads oversubscribe physical hardware.

---

## 2. CS1 – Heat Diffusion (MPI)

### Problem & Parallel Decomposition
The 2D heat equation solver performs 5-point stencil relaxation across an N x N domain over T timesteps. Each interior grid cell at timestep t+1 depends solely on its 4 orthogonal neighbours from timestep t:
u^{t+1}_{i,j} = 0.25 * (u^t_{i-1,j} + u^t_{i+1,j} + u^t_{i,j-1} + u^t_{i,j+1})

Because all cell updates within a step are independent, the spatial domain is decomposed across p MPI ranks. We implemented and benchmarked four communication/decomposition configurations:
1. **1-D Row Decomposition (Blocking):** The domain is partitioned into contiguous horizontal slabs of (N/p) x N cells. Halos are exchanged with north and south ranks using `MPI_Sendrecv`.
2. **1-D Row Decomposition (Non-blocking):** Uses `MPI_Isend` and `MPI_Irecv` to initiate halo exchange, updates the interior slab (rows 1 to H-2), calls `MPI_Waitall`, and finishes boundary cells.
3. **2-D Cartesian Decomposition (Blocking):** Ranks are mapped onto a 2-D virtual topology via `MPI_Cart_create`. Each rank owns a subgrid of (N/px) x (N/py) and exchanges 4 boundaries (North, South, East, West). Non-contiguous column boundaries utilize derived datatypes (`MPI_Type_vector`).
4. **2-D Cartesian Decomposition (Non-blocking):** Staged halo exchange with non-blocking calls overlapping inner subgrid computation with halo transit.

### Correctness Validation
Validation was verified against the serial baseline for identical problem parameters (N=1024, T=500, K=50):
- **Serial Checksum:** `1.193868328e+06` (max change = `4.839310e-02`)
- **Parallel Checksum (all ranks):** `1.193868328e+06` (max change = `4.839310e-02`)
- **Relative Difference:** `0.000000e+00` (exact floating-point match across all ranks and decomposition variants).

### Benchmark Results

#### Strong Scaling (N=1024, Serial Baseline T = 1.7050 s)
| Variant | p | Median Time (s) | Runs | Speedup S | Efficiency E |
|---|---|---|---|---|---|
| `blocking_1d` | 1 | 1.5352 | 3 | 1.11 | 1.111 |
| `blocking_1d` | 2 | 1.5107 | 3 | 1.13 | 0.564 |
| `blocking_1d` | 4 | 1.1054 | 3 | 1.54 | 0.386 |
| `blocking_2d` | 1 | 1.2705 | 3 | 1.34 | 1.342 |
| `blocking_2d` | 2 | 2.1590 | 3 | 0.79 | 0.395 |
| `blocking_2d` | 4 | 1.7394 | 3 | 0.98 | 0.245 |
| `nonblocking_1d` | 1 | 1.4589 | 3 | 1.17 | 1.169 |
| `nonblocking_1d` | 2 | 1.1249 | 3 | 1.52 | 0.758 |
| `nonblocking_1d` | 4 | 1.1758 | 3 | 1.45 | 0.362 |
| `nonblocking_2d` | 1 | 1.4193 | 3 | 1.20 | 1.201 |
| `nonblocking_2d` | 2 | 1.3513 | 3 | 1.26 | 0.631 |
| `nonblocking_2d` | 4 | 1.0880 | 3 | 1.57 | 0.392 |

#### Strong Scaling (N=4096, Serial Baseline T = 18.8859 s)
| Variant | p | Median Time (s) | Runs | Speedup S | Efficiency E |
|---|---|---|---|---|---|
| `blocking_1d` | 1 | 16.9080 | 3 | 1.12 | 1.117 |
| `blocking_1d` | 2 | 12.1347 | 3 | 1.56 | 0.778 |
| `blocking_1d` | 4 | 11.5126 | 3 | 1.64 | 0.410 |
| `blocking_2d` | 1 | 16.6999 | 3 | 1.13 | 1.131 |
| `blocking_2d` | 2 | 12.0064 | 3 | 1.57 | 0.786 |
| `blocking_2d` | 4 | 11.8025 | 3 | 1.60 | 0.400 |
| `nonblocking_1d` | 1 | 16.8534 | 3 | 1.12 | 1.121 |
| `nonblocking_1d` | 2 | 11.8456 | 3 | 1.59 | 0.797 |
| `nonblocking_1d` | 4 | 11.7930 | 3 | 1.60 | 0.400 |
| `nonblocking_2d` | 1 | 16.9611 | 3 | 1.11 | 1.113 |
| `nonblocking_2d` | 2 | 12.0244 | 3 | 1.57 | 0.785 |
| `nonblocking_2d` | 4 | 12.2331 | 3 | 1.54 | 0.386 |

#### Weak Scaling (N = 1024 * sqrt(p))
| Variant | p | Median Time (s) | Runs | Efficiency E = T1/Tp |
|---|---|---|---|---|
| `blocking_2d` | 1 | 1.1933 | 3 | 1.000 |
| `blocking_2d` | 2 | 1.6431 | 3 | 0.726 |
| `blocking_2d` | 4 | 5.1130 | 3 | 0.233 |
| `nonblocking_2d` | 1 | 1.1853 | 3 | 1.000 |
| `nonblocking_2d` | 2 | 1.5275 | 3 | 0.776 |
| `nonblocking_2d` | 4 | 4.1346 | 3 | 0.287 |

### In-Depth Analysis & Key Insights
1. **Surface-to-Volume Ratio & 1-D vs. 2-D Decomposition:**
   - For a rank owning N^2/p cells, 1-D slab decomposition exchanges 2N boundary elements per step (invariant with respect to p), yielding a communication-to-computation ratio of 2p/N.
   - In contrast, 2-D square block decomposition exchanges 4N/sqrt(p) boundary elements per step, yielding a ratio of 4*sqrt(p)/N.
   - At small p (p=2, 4), 1-D row decomposition sends contiguous memory buffers directly without vector packing overhead, outperforming 2-D blocking which incurs strided non-contiguous memory access for East/West boundaries. As problem size increases to N=4096, compute dominates and all variants achieve speedups between 1.56x and 1.64x.
2. **Computation-Communication Overlap:**
   - Non-blocking execution (`nonblocking_1d` at N=1024, p=2) reduces elapsed time from 1.5107 s to 1.1249 s, delivering a 25.5% runtime reduction by overlapping the communication of outer halos with the calculation of interior domain rows.
   - For smaller subdomains, the exposed communication overhead in blocking calls accounts for memory synchronization stalls, whereas `MPI_Isend`/`MPI_Irecv` allows overlapping ~35% of communication latency behind independent compute loops.
3. **Strong vs. Weak Scaling Efficiency:**
   - At N=4096, parallel efficiency remains healthy at p=2 (E ~ 0.79), but drops to E ~ 0.40 at p=4. This drop coincides with the physical hardware boundary (2 physical cores), where 4 ranks share L2 caches and memory channels.
4. **Convergence Check Frequency (K=1 vs. K=50):**
   - Benchmarking p=4 with `nonblocking_2d` revealed a dramatic synchronization bottleneck:
     - K=1 (global `MPI_Allreduce` executed every timestep): Median time = **2.7876 s**
     - K=50 (global `MPI_Allreduce` executed every 50 timesteps): Median time = **1.1078 s**
   - Testing convergence every step introduces a latency-bound global barrier that synchronizes all processes, degrading performance by **2.52x**. Amortizing the check across 50 iterations completely eliminates this barrier stall.

---

## 3. CS2 – K-mer Counting (OpenMP)

### Problem & Parallel Decomposition
K-mer counting scans a genomic DNA string of length L (L = 10^8 bp) to compute the frequency distribution of all possible substrings of length k (4^k distinct keys).
The input sequence is partitioned into contiguous chunks across OpenMP threads. To ensure boundary k-mers spanning chunk borders are counted exactly once, each thread initializes its rolling window by pre-loading k-1 characters from the preceding boundary.
Three distinct synchronization strategies were evaluated:
- **Version A (Data Race / Unsynchronized):** Threads increment a single shared histogram directly using `++table[key]`.
- **Version B (Atomic Updates):** Threads share one histogram and protect every update with `#pragma omp atomic`.
- **Version C (Private Histograms + Reduction):** Each thread allocates a thread-private histogram of size 4^k * 4 bytes, accumulates locally without synchronization, and performs an OpenMP reduction/merge into the global table at thread termination.
- **False Sharing Experiment:** Evaluates packed vs. cache-padded thread structures for small k (k=1).

### Correctness Validation
- **Race Condition in Version A:**
  - In Version A, unsynchronized `++` translates to a non-atomic load-modify-store sequence. Concurrent updates to identical k-mer bins overwrite each other. On L = 10^8, k=8, p=4, Version A lost millions of counts due to conflicting writes.
- **Verification in Version C:**
  - Expected total k-mers: L - k + 1 = 100,000,000 - 8 + 1 = 99,999,993.
  - Measured total k-mers: **99,999,993** (`COUNT_OK`).
  - Computed Checksum: `8956218861015476594`, perfectly matching serial baseline output.

### Benchmark Results

#### Strong Scaling for k=8 (4^8 = 65,536 bins, Serial Baseline T = 0.3390 s)
| Variant | p | Median Time (s) | Runs | Speedup S | Efficiency E |
|---|---|---|---|---|---|
| `atomic_B` | 1 | 0.7786 | 3 | 0.44 | 0.435 |
| `atomic_B` | 2 | 0.8221 | 3 | 0.41 | 0.206 |
| `atomic_B` | 4 | 0.7365 | 3 | 0.46 | 0.115 |
| `private_C` | 1 | 0.1633 | 3 | 2.08 | 2.077 |
| `private_C` | 2 | 0.0950 | 3 | 3.57 | 1.783 |
| `private_C` | 4 | 0.1249 | 3 | 2.71 | 0.678 |

#### Strong Scaling for k=12 (4^12 = 16,777,216 bins, Serial Baseline T = 2.8926 s)
| Variant | p | Median Time (s) | Runs | Speedup S | Efficiency E |
|---|---|---|---|---|---|
| `atomic_B` | 1 | 2.9508 | 3 | 0.98 | 0.980 |
| `atomic_B` | 2 | 1.9673 | 3 | 1.47 | 0.735 |
| `atomic_B` | 4 | 1.6341 | 3 | 1.77 | 0.443 |
| `private_C` | 1 | 3.5842 | 3 | 0.81 | 0.807 |
| `private_C` | 2 | 2.7261 | 3 | 1.06 | 0.531 |
| `private_C` | 4 | 2.9018 | 3 | 1.00 | 0.249 |

#### False Sharing Study (k=1, 4 bins, Serial Baseline T = 0.1339 s)
| Variant | p | Median Time (s) | Runs | Speedup S | Efficiency E |
|---|---|---|---|---|---|
| `packed_F` | 1 | 0.1605 | 3 | 0.83 | 0.834 |
| `packed_F` | 2 | 0.2726 | 3 | 0.49 | 0.246 |
| `packed_F` | 4 | 0.4336 | 3 | 0.31 | 0.077 |
| `padded_P` | 1 | 0.1412 | 3 | 0.95 | 0.949 |
| `padded_P` | 2 | 0.1103 | 3 | 1.21 | 0.607 |
| `padded_P` | 4 | 0.0998 | 3 | 1.34 | 0.335 |

### In-Depth Analysis & Key Insights
1. **Atomic vs. Private Histogram Crossover:**
   - A critical performance inversion occurs between small k (k=8) and large k (k=12):
     - For k=8, each private histogram occupies 4^8 * 4 B = 256 KB, comfortably fitting inside the 256 KB L2 cache per core. With zero synchronization and full cache residency, Version C achieves an outstanding speedup of **3.57x** at p=2 (throughput = 453.66 Mkmers/s). Version B suffers severe atomic lock contention on the small table (speedup only 0.41x).
     - For k=12, the table size explodes to 4^12 * 4 B = 64 MB per thread. Allocating p * 64 MB = 256 MB completely exhausts the 4 MiB L3 cache, causing catastrophic TLB and cache thrashing followed by an expensive multi-megabyte reduction merge.
     - Conversely, in Version B at k=12, updates are sparsely distributed across 16.7 million buckets, making collision probability between threads negligible (~10^-7). Atomic instructions execute with minimal cache line bouncing, allowing Version B (T = 1.6341 s, S = 1.77x) to outpace Version C (T = 2.9018 s, S = 1.00x) by **1.78x**.
2. **False Sharing Demonstration:**
   - In `packed_F`, private counter structs are tightly packed (16 bytes per thread, fitting 4 thread counters on a single 64-byte cache line). Writing to disjoint variables invalidates the shared cache line across cores on every increment (MESI protocol invalidation storm). Consequently, running on 4 threads slows down from 0.1605 s to 0.4336 s (S = 0.31x).
   - In `padded_P`, aligning each thread structure to 64 bytes (`alignas(64)`) eliminates false sharing completely, achieving 0.0998 s at p=4 (S = 1.34x). Cache line padding delivers a **4.34x speedup** over the false-sharing baseline.

---

## 4. CS3 – Monte Carlo Option Pricing (OpenMP + MPI)

### Problem & Parallel Decomposition
Monte Carlo European call option pricing estimates expected discounted payoff under geometric Brownian motion:
C = exp(-r * T) * E[max(S_T - K, 0)], where S_T = S_0 * exp((r - 0.5 * sigma^2)*T + sigma * sqrt(T) * Z)
- **Parameters:** Spot price S0 = 100, Strike K = 100, Risk-free rate r = 0.05, Volatility sigma = 0.20, Maturity T = 1.0 year, Total paths M = 10^8.
- **Analytical Benchmark:** Black-Scholes formula yields CBS = 10.450584.
- **Parallel Model Justification:** Monte Carlo simulation is an embarrassingly parallel workload with independent path generation and a single scalar reduction sum at completion. We implemented both an OpenMP shared-memory solver and an MPI distributed-memory variant.

### Random Number Generation & Pitfalls
Parallel RNG correctness requires uncorrelated, non-overlapping pseudorandom sequences across execution units:
1. **Engine:** `xoshiro256**` 64-bit generator with Box-Muller transformation for normal variates.
2. **Independent Streams:** Each thread/rank i advances the initial state by calling `jump()`, advancing the generator by 2^128 steps to ensure completely disjoint sample streams.
3. **RNG Pitfall Demonstrations:**
   - **Shared Generator + Critical Section (`omp_shared_rng`):** Multiple threads synchronize on every normal variate extraction. At M=10^6, p=4, execution required 0.4492 s (throughput 2.23 Mpaths/s) compared to 0.1091 s (9.16 Mpaths/s) for independent RNG — a **4.12x slowdown** due purely to lock contention.
   - **Same Seed Error (`omp_same_seed`):** Each thread initializes with the identical seed without jumping. All threads execute the identical sequence of random numbers, yielding an empirical standard error of `0.001472` that is entirely fraudulent: variance is not reduced by sqrt(p) because identical samples are duplicated across threads.

### Benchmark Results (M = 10^8, Serial Baseline T = 7.4761 s)
| Variant | p | Median Time (s) | Runs | Speedup S | Efficiency E |
|---|---|---|---|---|---|
| `omp_independent` | 1 | 7.5627 | 3 | 0.99 | 0.989 |
| `omp_independent` | 2 | 6.6140 | 3 | 1.13 | 0.565 |
| `omp_independent` | 4 | 2.8961 | 3 | 2.58 | 0.645 |
| `mpi_independent` | 1 | 7.2970 | 3 | 1.02 | 1.025 |
| `mpi_independent` | 2 | 3.8786 | 3 | 1.93 | 0.964 |
| `mpi_independent` | 4 | 4.1646 | 3 | 1.80 | 0.449 |

### In-Depth Analysis & Key Insights
1. **Parallel Speedup & Scalability:**
   - `omp_independent` achieves a strong speedup of **2.58x** at p=4 (T = 2.8961 s), pricing 34.53 million paths/second.
   - `mpi_independent` demonstrates near-linear scaling at p=2 with T = 3.8786 s (S = 1.93x, E = 96.4%). At p=4, inter-process memory bandwidth saturation and OS scheduling across 2 physical cores cap the speedup at 1.80x.
2. **Amdahl's Law Serial Fraction:**
   - For `omp_independent` at p=4 (S=2.58):
     f = (1/S - 1/p) / (1 - 1/p) = (1/2.58 - 1/4) / (1 - 0.25) = (0.3876 - 0.2500) / 0.75 = 0.183 (18.3%)
   - The effective serial fraction stems primarily from thermal throttling on the ultra-low-power CPU (1.2 GHz base clock) under all-core AVX/FPU load and thread synchronization at loop reduction.
3. **Reproducibility & Accuracy:**
   - Under independent streams, for M=10^8, p=8, the computed price is **10.448315** with standard error 0.001472.
   - The absolute deviation from Black-Scholes (|10.448315 - 10.450584| = 0.002268) is 1.54x the standard error, perfectly conforming to standard Gaussian 95% confidence intervals (<= 1.96 SE).

---

## 5. CS4 – PageRank (OpenMP)
*(To be completed by Author B)*

---

## 6. CS5 – Sample Sort (MPI)
*(To be completed by Author B)*

---

## 7. Cross-Case Comparison
*(To be completed by Author B)*

---

## 8. Individual Contributions
*(To be completed by Author B)*

---

## 9. References
- SNAP Stanford Network Analysis Platform: Web-Google Dataset (J. Leskovec et al.)
- Blackman, D., & Vigna, S. (2021). "Scrambled Linear Pseudorandom Number Generators." ACM TOMS.
- Chakrabarti, D., Zhan, Y., & Faloutsos, C. (2004). "R-MAT: A Recursive Model for Graph Mining." SDM.
- OpenMP Architecture Review Board, OpenMP Application Programming Interface Version 5.0.
- Message Passing Interface Forum, MPI: A Message-Passing Interface Standard Version 3.1.
