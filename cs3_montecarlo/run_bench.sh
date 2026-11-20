#!/bin/bash
# CS3 benchmarks. Usage: ./run_bench.sh   (env: RUNS, MAXP, MPIRUN, M)
cd "$(dirname "$0")"; source ../tools/lib.sh
mkdir -p results; M=${M:-100000000}; S=results/cs3_strong.csv
rm -f results/cs3_*.csv results/cs3_*.log
bench M$M serial 1 $S ./serial $M 42
for p in $(plist); do
  bench M$M omp_independent $p $S ./parallel $M $p 42 0
  [ -x ./parallel_mpi ] && bench M$M mpi_independent $p $S $MPIRUN -np $p ./parallel_mpi $M 42
done
echo "== RNG pitfalls (small M) =="
./parallel 1000000 $MAXP 42 0 | tee results/cs3_rng_good.txt
./parallel 1000000 $MAXP 42 1 | tee results/cs3_rng_shared.txt
./parallel 100000000 $MAXP 42 2 | tee results/cs3_rng_sameseed.txt
echo "== reproducibility: same p twice, then a different p =="
./parallel $M 4 42 0; ./parallel $M 4 42 0; ./parallel $M 8 42 0 | tee results/cs3_repro.txt
python3 ../tools/analyze.py $S
