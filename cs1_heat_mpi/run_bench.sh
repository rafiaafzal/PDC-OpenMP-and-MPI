#!/bin/bash
# CS1 benchmarks. Usage: ./run_bench.sh   (env: RUNS, MAXP, MPIRUN, T, SIZES)
cd "$(dirname "$0")"; source ../tools/lib.sh
mkdir -p results; T=${T:-500}; SIZES=${SIZES:-"1024 4096"}
S=results/cs1_strong.csv; W=results/cs1_weak.csv; Kc=results/cs1_K.csv
rm -f results/cs1_*.csv results/cs1_*.log
for N in $SIZES; do
  bench N$N serial 1 $S ./serial $N $T 50
  for p in $(plist); do
    for dim in 1 2; do
      bench N$N blocking_${dim}d    $p $S $MPIRUN -np $p ./parallel $N $T 50 0 $dim
      bench N$N nonblocking_${dim}d $p $S $MPIRUN -np $p ./parallel $N $T 50 1 $dim
    done
  done
done
# weak scaling: ~1024x1024 cells per rank  ->  N = 1024*sqrt(p)
for p in $(plist); do
  N=$(awk -v p=$p "BEGIN{printf \"%d\", 1024*sqrt(p)+0.5}")
  bench weak blocking_2d    $p $W $MPIRUN -np $p ./parallel $N $T 50 0 2
  bench weak nonblocking_2d $p $W $MPIRUN -np $p ./parallel $N $T 50 1 2
done
# cost of convergence checks: K=1 vs K=50 at the largest p
N=$(echo $SIZES | awk "{print \$1}")
for K in 1 50; do bench K$K nonblocking_2d $MAXP $Kc $MPIRUN -np $MAXP ./parallel $N $T $K 1 2; done
python3 ../tools/analyze.py $S; python3 ../tools/analyze.py $W --weak
echo "K-check rows (compare time and comm fields):"; grep -h "K=" results/cs1_K.log | head -6
