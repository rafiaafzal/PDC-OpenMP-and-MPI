# CS3006 Assignment 3 – Parallel Case Studies (OpenMP + MPI)

Model choices (rule: >=2 OpenMP, >=2 MPI):

| Case | Problem | Model | Directory |
|---|---|---|---|
| 1 | Heat diffusion (5-point stencil) | **MPI** (1-D and 2-D decomposition, blocking + non-blocking halo exchange) | `cs1_heat_mpi` |
| 2 | k-mer counting | **OpenMP** (race / atomic / private histograms / false sharing) | `cs2_kmer_openmp` |
| 3 | Monte-Carlo option pricing | **OpenMP** (main) + **MPI** variant | `cs3_montecarlo` |
| 4 | PageRank | **OpenMP** (pull vs push, schedules) | `cs4_pagerank_openmp` |
| 5 | Sample sort | **MPI** | `cs5_samplesort_mpi` |

Replace `Group00` with your group number when you rename the folder / zip.

## 1. Setup (Linux / WSL / macOS)
```
sudo apt install build-essential openmpi-bin libopenmpi-dev python3 python3-matplotlib   # Ubuntu/WSL
mpicc --version; mpirun --version; gcc --version; lscpu | head -20   # record CPU model/cores/RAM for the report
```
macOS: `brew install gcc open-mpi` and use `gcc-14` (set `CC=gcc-14`); Apple clang has no `-fopenmp`.
Native Windows: use WSL (recommended) or MS-MPI.

## 2. Build everything
```
make            # from this folder; builds serial + parallel in all five case studies
```
(Each case study also has its own Makefile: `cd cs1_heat_mpi && make`.)

## 3. Run the experiments (writes CSV + tables + PNG graphs into each `results/`)
Close other programs first (idle machine). Defaults: 3 runs per point, p = 1,2,4,... up to `nproc`.
```
cd cs1_heat_mpi         && ./run_bench.sh
cd ../cs2_kmer_openmp   && ./run_bench.sh
cd ../cs3_montecarlo    && ./run_bench.sh
cd ../cs4_pagerank_openmp && ./run_bench.sh
cd ../cs5_samplesort_mpi  && ./run_bench.sh
```
Environment knobs (all scripts): `RUNS=5`, `MAXP=8` (max threads/ranks), `MPIRUN="mpirun --oversubscribe"`
(only if you deliberately use more ranks than slots; say so in the report). Smaller/faster test run:
`RUNS=1 MAXP=2 ./run_bench.sh`. Memory-hungry defaults can be reduced with `N=`, `L=`, `M=`, `SCALE=`.

Logical vs physical cores: speedup beyond the physical core count is not meaningful. Set `MAXP` to the
number of **physical** cores for the main strong-scaling graphs.

## 4. What each folder contains
`serial.*` baseline, `parallel.*` parallel version(s), `Makefile`, `README.md` (how to run + what to measure
for the analysis questions), `run_bench.sh`. Every program prints one machine-readable `RESULT ...` line.
`tools/` holds `analyze.py` (speedup/efficiency tables + graphs), `plot_phases.py`, `plot_degree.py`,
`fit_alpha_beta.py`.

## 5. What you still have to do yourselves
1. Run everything on your machine and look at the numbers (graphs are produced automatically).
2. Fill in `report/GroupXX_report_TEMPLATE.md` with your measured numbers, export to PDF
   (Word/Google Docs/`pandoc report.md -o GroupXX_report.pdf`). The template already contains the derivations
   and the explanations the questions ask for; confirm each claim against YOUR data and rewrite in your words.
3. Build the 6-slide deck from `slides/slides_outline.md` (PowerPoint/Google Slides -> PDF).
4. Fill in `CONTRIBUTIONS.md`, make real Git commits (both members, over the two weeks), submit peer forms.
5. `zip -r GroupXX_Assignment3.zip GroupXX_Assignment3 -x "*/results/*.log" "*.o"` (keep results/ CSV + PNG).
