#!/bin/bash
# CS5 benchmarks. Usage: ./run_bench.sh   (env: RUNS, MAXP, MPIRUN, N)   N=1e8 needs ~1.5 GB RAM in total
cd "$(dirname "$0")"; source ../tools/lib.sh
mkdir -p results; N=${N:-100000000}; S=results/cs5_strong.csv
rm -f results/cs5_*.csv results/cs5_*.log results/cs5_phases.log
for dist in 0 1 2; do
  bench dist$dist serial 1 $S ./serial $N $dist 1
  for p in $(plist); do
    [ $p -eq 1 ] && continue
    bench dist$dist samplesort_os1 $p $S $MPIRUN -np $p ./parallel $N $dist 1
  done
done
cp results/cs5_strong.log results/cs5_phases.log
echo "== oversampling effect at p=$MAXP =="
for os in 1 4 16 64; do $MPIRUN -np $MAXP ./parallel $N 1 $os 1 | tee -a results/cs5_oversample.log; done
echo "== checksum cross-check against serial (same data) =="
./serial 2000000 1 $MAXP; $MPIRUN -np $MAXP ./parallel 2000000 1 1
echo "== ping-pong for alpha/beta =="
$MPIRUN -np 2 ./pingpong > results/pingpong.csv && python3 ../tools/fit_alpha_beta.py results/pingpong.csv $N $MAXP
python3 ../tools/analyze.py $S
python3 ../tools/plot_phases.py results/cs5_phases.log 0
