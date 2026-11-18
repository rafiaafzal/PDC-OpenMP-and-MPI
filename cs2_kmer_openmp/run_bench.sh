#!/bin/bash
# CS2 benchmarks. Usage: ./run_bench.sh   (env: RUNS, MAXP, L)
cd "$(dirname "$0")"; source ../tools/lib.sh
mkdir -p results; L=${L:-100000000}; S=results/cs2_strong.csv; F=results/cs2_falsesharing.csv
rm -f results/cs2_*.csv results/cs2_*.log
echo "== correctness (Version A should say COUNT_WRONG on a multi-core machine) =="
./serial $L 8 0; ./parallel $L 8 0 A $MAXP; ./parallel $L 8 0 B $MAXP; ./parallel $L 8 0 C $MAXP | tee results/cs2_correctness.txt
for k in 8 12; do
  bench k$k serial 1 $S ./serial $L $k 0
  for p in $(plist); do
    bench k$k atomic_B  $p $S ./parallel $L $k 0 B $p
    bench k$k private_C $p $S ./parallel $L $k 0 C $p
  done
done
# false sharing: k=1 (4 bins x 4B = 16 B per thread) packed vs padded
bench k1 serial 1 $F ./serial $L 1 0
for p in $(plist); do
  bench k1 packed_F $p $F ./parallel $L 1 0 F $p
  bench k1 padded_P $p $F ./parallel $L 1 0 P $p
done
python3 ../tools/analyze.py $S; python3 ../tools/analyze.py $F
