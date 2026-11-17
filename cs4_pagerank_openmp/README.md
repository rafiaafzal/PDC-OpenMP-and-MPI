# CS4 – PageRank (OpenMP)

Files: `pr_common.h` (graph generator/loader, CSR, serial PageRank), `serial.c`, `parallel.c`, `Makefile`, `run_bench.sh`.
Graph: synthetic R-MAT power-law graph (a,b,c,d = .57,.19,.19,.05, seed fixed), default scale 20 -> 1,048,576 vertices and
8,388,608 edges (>1e6). Hub ids are not scrambled (hubs have low ids), which makes `schedule(static)` imbalance clearly
visible. Or load a SNAP file: `-f web-Google.txt` (https://snap.stanford.edu/data/web-Google.html, cite it). CSR for out-edges (push)
and in-edges (pull) are both built. d=0.85, eps=1e-8 (L1), dangling nodes spread uniformly each iteration.

Build: `make`.  Run:
```
./serial [-s 20] [-e 8] [-f file] [-d results/deg.csv]
./parallel -m pull|push -S static|dynamic|guided -c CHUNK -t THREADS -v 1 [-s 20 -e 8 | -f file]
./parallel -m pull -S dynamic -c 64 -t 8 -v 1      # example (-v 1 also checks against serial in-process)
```
Verification: `l1_vs_serial` must be < 1e-6 (prints `VERIFY_OK`).

Output fields: `imbalance_R` = slowest thread time / mean thread time in the main vertex loop (measured);
`static_edge_imbalance` = same ratio computed analytically for contiguous equal-size blocks (edges per thread);
`medges_s` = edges processed per second; `est_GBs` = estimated traffic using 12 B/edge for pull
(4 B column index + 8 B random read of contrib[u]) and 20 B/edge for push (4 B index + 8 B read + 8 B write of the atomic RMW target),
ignoring sequential per-vertex arrays. Compare with your measured memory bandwidth (e.g. STREAM) for the memory-bound argument.

Analysis hints: degree plot -> `python3 ../tools/plot_degree.py results/deg.csv` (made by run_bench.sh);
schedule table -> `results/cs4_sched.csv/log`; pull vs push -> `results/cs4_strong_table.md`.
Equal vertex counts != equal work because work per vertex ~ in-degree, which is power-law distributed.
Small dynamic chunks balance well but pay scheduling overhead; large chunks reduce overhead but reintroduce imbalance;
`guided` starts large and shrinks. Push adds atomic RMW on random addresses (contention on hubs, cache-line ping-pong).
