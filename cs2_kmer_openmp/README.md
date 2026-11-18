# CS2 – k-mer counting (OpenMP)

Files: `common.h` (timer, RNG, DNA generator), `serial.c`, `parallel.c`, `Makefile`, `run_bench.sh`.
Build: `make`  (`gcc -fopenmp`).

Run:
```
./serial L k skew                       # L=100000000, k in {8,12}, skew 0 uniform / 1 skewed bases
./parallel L k skew MODE threads        # MODE: A race(wrong) | B atomic | C private histograms | F packed | P padded
./parallel 100000000 8 0 A 8            # shows COUNT_WRONG
```
The sequence is generated in memory (1 byte per base value 0..3, fixed seed) and generation is not timed.
Verification: `total` must equal `expected` = L-k+1 and `checksum` must equal the serial checksum (A fails, B/C/F/P pass).
Chunk boundaries: k-mer *start* positions are split among threads; each thread primes its rolling window with the
k-1 bases before its first start position, so a k-mer straddling two chunks is counted once, by the thread that owns its first base.

Analysis hints:
1. Race: A increments a shared counter with a non-atomic load-add-store; two threads read the same old value and one
   increment is lost. Report `total` vs `expected`.
2. Atomic vs private: C extra memory = threads * 4^k * 4 bytes (printed as `extraMB`): k=8 -> 256 KB/thread (fits L2,
   C wins easily), k=12 -> 64 MB/thread (exceeds cache, merge cost grows, atomic contention shrinks because 16M bins are
   rarely hit by two threads at once) -> crossover. Look at `results/cs2_strong_table.md`.
3. False sharing: `F` vs `P` at k=1 (each thread histogram = 16 B, four threads share one 64 B line). `P` pads to 64 B.
4. Bytes per k-mer: 1 B streamed from the sequence + a 4 B read-modify-write on the histogram (cache-resident for small k,
   ~64 B cache line traffic when random over a large histogram). Throughput printed as `Mkmers_s`; bandwidth ~ Mkmers_s * bytes.
   Compare with your machine memory bandwidth to see where scaling flattens.
