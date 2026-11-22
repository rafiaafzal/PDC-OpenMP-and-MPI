# CONTRIBUTIONS

**Group:** Group 00  
**Members:**
- **Member A (Lead CS1, CS2, CS3):** Rafia Afzal (Roll No: 23L-0821, Git: `rafiaafzal`)
- **Member B (Lead CS4, CS5):** Maheen Kamran (Roll No: 23L-0715, Git: `maheenkamran`)

---

## Lead / Reviewer Assignment
| Case Study | Topic & Framework | Lead | Reviewer |
|---|---|---|---|
| **CS1** | Heat Diffusion (MPI 1-D / 2-D Stencil) | Rafia Afzal (A) | Maheen Kamran (B) |
| **CS2** | K-mer Counting (OpenMP Atomics vs Private) | Rafia Afzal (A) | Maheen Kamran (B) |
| **CS3** | Monte Carlo Option Pricing (OpenMP + MPI) | Rafia Afzal (A) | Maheen Kamran (B) |
| **CS4** | PageRank (OpenMP Push vs Pull & Sched) | Maheen Kamran (B) | Rafia Afzal (A) |
| **CS5** | Sample Sort (MPI Regular Sampling) | Maheen Kamran (B) | Rafia Afzal (A) |

---

## Member A — Rafia Afzal (`rafiaafzal`)
- **Roles & Ownership:** Lead on CS1, CS2, CS3; Project Architecture; Cross-case tools.
- **Tasks Completed:**
  - Setup repository skeleton, top-level `Makefile`, `.gitignore`, and shared analysis tools (`tools/lib.sh`, `tools/analyze.py`, etc.).
  - Implemented serial, 1-D blocking/non-blocking, and 2-D Cartesian blocking/non-blocking MPI solvers for CS1 (`cs1_heat_mpi/`).
  - Implemented strong and weak scaling benchmarks for CS1, testing convergence frequency ($K=1$ vs. $K=50$).
  - Implemented OpenMP K-mer counting versions A (data race demo), B (atomic updates), and C (private histograms + reduction) for CS2 (`cs2_kmer_openmp/`).
  - Implemented false sharing evaluation (packed vs. 64-byte padded counters) for CS2.
  - Implemented parallel Monte Carlo European call option pricing with `xoshiro256**` jump streams in OpenMP and MPI (`cs3_montecarlo/`).
  - Implemented RNG pitfalls demonstrating lock contention in shared generators and variance distortion with identical seeds.
  - Authored full methodology, performance analysis, and tables for CS1, CS2, and CS3 in `report/GroupXX_report.md`.
- **Approx. Hours:** ~28 hours.
- **Commits Authored:** Commits 1, 2, 3, 5, 7, 9, 11, 12.

---

## Member B — Maheen Kamran (`maheenkamran`)
- **Roles & Ownership:** Lead on CS4, CS5; Report synthesis and presentation deck.
- **Tasks Completed:**
  - Implemented serial and OpenMP PageRank solver supporting both Pull (lock-free) and Push (atomic) traversal models (`cs4_pagerank_openmp/`).
  - Benchmarked thread scheduling strategies (`static` vs. `dynamic` with chunk sizes 1 to 1024) across scale-free R-MAT graph degree distributions.
  - Implemented serial and distributed-memory MPI sample sort with regular sampling, all-to-all communications (`MPI_Alltoallv`), and multi-way local merge (`cs5_samplesort_mpi/`).
  - Implemented MPI ping-pong latency and bandwidth microbenchmark to model `MPI_Alltoallv` communication performance.
  - Executed benchmark suites for CS4 and CS5 across uniform, skewed, and duplicate key distributions.
  - Authored full methodology, performance analysis, and tables for CS4 and CS5 in `report/GroupXX_report.md`.
  - Authored Section 7 (Cross-case comparison), Section 8 (Contributions), and structured presentation deck outline (`slides/slides_outline.md`).
- **Approx. Hours:** ~26 hours.
- **Commits Authored:** Commits 4, 6, 8, 10, 13, 14.

---

## Joint Work
- Cross-case performance comparison and architectural synthesis across shared-memory (OpenMP) and distributed-memory (MPI) models.
- End-to-end clean build integration verification and test script execution on Ubuntu 24.04 LTS (WSL2).
- Presentation slide outline alignment and data validation against measured benchmark metrics.
