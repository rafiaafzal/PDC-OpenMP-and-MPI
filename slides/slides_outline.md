# Slide Deck Outline: Exactly 6 Slides (1 Intro + 5 Case Studies)
Group: Group 00 | Course: Parallel and Distributed Computing

---

## Slide 1 – Introduction
- **Title:** Parallel Case Studies with OpenMP and MPI
- **Group:** Group 00
- **Team Members & Responsibilities:**
  - **Rafia Afzal (23L-0821):** Lead for CS1 (Heat MPI), CS2 (K-mer OpenMP), CS3 (Monte Carlo OpenMP/MPI), Build & Tools Skeleton.
  - **Maheen Kamran (23L-0715):** Lead for CS4 (PageRank OpenMP), CS5 (Sample Sort MPI), Cross-case Synthesis & Slides.
- **Hardware & Software Environment:**
  - **CPU:** Intel(R) Core(TM) m7-6Y75 @ 1.20GHz (2 physical cores, 4 logical execution threads, SMT enabled).
  - **Caches & Memory:** 32 KiB L1d/core, 256 KiB L2/core, 4 MiB L3 shared; 7.5 GiB RAM; Ubuntu 24.04 LTS (WSL2).
  - **Toolchain:** GCC/G++ 13.3.0 (`-O3 -march=native`), OpenMPI 4.1.6.
- **Parallel Programming Model Summary:**
  - **CS1:** MPI 2-D Cartesian halo exchange with non-blocking computation-communication overlap.
  - **CS2:** OpenMP private histogram buffers with terminal reduction vs. atomic updates.
  - **CS3:** OpenMP & MPI embarrassingly parallel simulation with `xoshiro256**` jumped streams.
  - **CS4:** OpenMP destination-driven lock-free Pull model with dynamic chunk scheduling.
  - **CS5:** MPI regular sample sort with `MPI_Alltoallv` partitioned data migration.

---

## Slide 2 – CS1: Heat Diffusion (MPI)
1. **Problem Statement:** 5-point stencil solver on 2D grid ($N=1024$ and $N=4096$, $T=500$ steps) via distributed-memory domain decomposition.
2. **Decomposition Diagram:** 1-D row slabs vs. 2-D virtual Cartesian grid (`MPI_Cart_create`) exchanging North, South, East, West ghost boundary cells.
3. **Speedup & Scaling:**
   - At $N=4096$, non-blocking 1-D achieves **1.64x speedup** ($p=4$, baseline $T = 18.89$ s vs. $11.51$ s).
   - Weak scaling efficiency drops to $28.7\%$ at $p=4$ due to dual-core SMT resource contention.
4. **Key Measured Bottleneck:**
   - Testing convergence every step ($K=1$) incurs global synchronization stalls (**2.79 s**). Amortizing checks every 50 steps ($K=50$) reduces runtime to **1.11 s** (**2.52x faster**). Non-blocking `MPI_Isend`/`MPI_Irecv` overlaps ~35% of boundary communication.
5. **Next Steps:** Implement hybrid MPI + OpenMP to eliminate intra-node message passing; test periodic multi-step halo buffering.

---

## Slide 3 – CS2: K-mer Counting (OpenMP)
1. **Problem Statement:** Count frequency of all substrings of length $k$ across a $10^8$ base-pair genome using OpenMP shared memory.
2. **Decomposition Diagram:** Contiguous sequence chunking with $k-1$ base pre-loading across thread boundaries to eliminate boundary undercounts.
3. **Speedup & Scaling:**
   - For $k=8$ ($65{,}536$ bins): Private histograms (`Version C`) achieve **3.57x speedup** at $p=2$ (throughput $453.66$ Mkmers/s).
   - For $k=12$ ($16.7$M bins): Atomic table (`Version B`) achieves **1.77x speedup** at $p=4$, outperforming private histograms.
4. **Key Measured Bottleneck:**
   - **Cache Hierarchy Crossover:** At $k=8$, private tables ($256$ KB/thread) reside entirely in L2 cache. At $k=12$, private tables explode to $64$ MB/thread ($256$ MB total), overflowing the $4$ MiB L3 cache and degrading performance.
   - **False Sharing:** Padded counters (`alignas(64)`) achieve **4.34x speedup** over tightly packed thread counters ($0.0998$ s vs. $0.4336$ s at $p=4$).
5. **Next Steps:** Implement 2-bit DNA base packing and cache-blocked k-mer bucketing.

---

## Slide 4 – CS3: Monte Carlo Option Pricing (OpenMP + MPI)
1. **Problem Statement:** Estimate European call option pricing ($M=10^8$ paths, $S_0=100, K=100, r=0.05, \sigma=0.2, T=1$) against Black-Scholes ($10.4506$).
2. **Decomposition Diagram:** Disjoint path simulation per thread/process using `xoshiro256**` pseudorandom streams advanced by $2^{128}$ jumps, followed by scalar reduction.
3. **Speedup & Scaling:**
   - OpenMP independent streams achieve **2.58x speedup** at $p=4$ ($T = 2.896$ s, pricing $34.53$ Mpaths/s).
   - MPI distributed variant achieves **1.93x speedup** at $p=2$ ($96.4\%$ efficiency).
4. **Key Measured Bottleneck:**
   - **Concurrency Pitfalls:** A shared RNG with critical sections causes severe lock contention (**4.12x slowdown** at $M=10^6$). Reusing identical seeds without jumping generates false standard error estimates with duplicate paths.
   - Amdahl's Law serial fraction $f = 18.3\%$ due to SMT core sharing and reduction overhead.
5. **Next Steps:** Implement antithetic variates variance reduction and SIMD vectorization (AVX2 Box-Muller).

---

## Slide 5 – CS4: PageRank (OpenMP)
1. **Problem Statement:** Iterative stationary distribution solver on scale-free R-MAT graph ($1{,}048{,}576$ vertices, $8{,}388{,}608$ directed edges) stored in CSR format.
2. **Decomposition Diagram:** Push model (atomic source writes) vs. Pull model (lock-free destination reads).
3. **Speedup & Scaling:**
   - `pull_dynamic64` achieves **1.61x speedup** at $p=4$ ($1.1793$ s vs. serial $1.9002$ s).
   - `push_dynamic64` suffers severe slowdown (**0.62x**, $3.0519$ s at $p=4$).
4. **Key Measured Bottleneck:**
   - **Push vs. Pull Contention:** High-degree hub vertices in heavy-tailed scale-free graphs cause catastrophic hardware atomic lock contention in Push (**2.59x slower** than Pull).
   - **Dynamic Chunk Scheduling:** Naive dynamic scheduling ($chunk=1$) incurs high loop dispatch overhead ($2.395$ s), whereas $chunk=256$ delivers optimal load balancing (**0.967 s**, **1.97x speedup**).
5. **Next Steps:** Apply vertex degree reordering (Hub sorting) and cache-blocked CSR segmenting.

---

## Slide 6 – CS5: Sample Sort (MPI)
1. **Problem Statement:** Scalable distributed sorting of $N = 10^7$ 32-bit keys across distributed-memory MPI processes.
2. **Decomposition Diagram:** 5-phase pipeline: Local sort $\to$ Regular sampling & splitters $\to$ Binary search partitioning $\to$ `MPI_Alltoallv` exchange $\to$ Local multi-way merge.
3. **Speedup & Scaling:**
   - Uniform distribution (`dist0`) achieves **2.93x speedup** at $p=4$ ($T = 6.938$ s vs. serial $20.354$ s).
   - Superlinear speedup of **2.30x at $p=2$** ($E = 115.2\%$) due to $O(N \log N)$ reduction in local sorting partition size.
4. **Key Measured Bottleneck:**
   - **Network & Phase Transition:** Ping-pong microbenchmark yields latency $\alpha = 1.16\ \mu$s and bandwidth $\approx 3.25$ GB/s. At $p=2$, Local Sort accounts for $75\%$ of runtime; at $p=4$, `MPI_Alltoallv` and merge grow to $38\%$ of runtime, becoming the primary scaling bottleneck.
5. **Next Steps:** Implement duplicate-aware splitter selection with rank tie-breaking and multi-threaded local sorting.
