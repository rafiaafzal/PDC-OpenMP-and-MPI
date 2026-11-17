#!/bin/bash
# CS4 benchmarks. Usage: ./run_bench.sh   (env: RUNS, MAXP, SCALE, EF, GRAPHARGS)
# GRAPHARGS example for a SNAP file:  GRAPHARGS="-f web-Google.txt" ./run_bench.sh
cd "$(dirname "$0")"; source ../tools/lib.sh
mkdir -p results; G=${GRAPHARGS:-"-s ${SCALE:-20} -e ${EF:-8}"}
S=results/cs4_strong.csv; C=results/cs4_sched.csv
rm -f results/cs4_*.csv results/cs4_*.log
./serial $G -d results/deg.csv | tee results/cs4_serial.txt
python3 ../tools/plot_degree.py results/deg.csv
bench g serial 1 $S ./serial $G
for p in $(plist); do
  bench g pull_static  $p $S ./parallel $G -m pull -S static -t $p -v 0
  bench g pull_dynamic64 $p $S ./parallel $G -m pull -S dynamic -c 64 -t $p -v 0
  bench g push_dynamic64 $p $S ./parallel $G -m push -S dynamic -c 64 -t $p -v 0
done
echo "== correctness check (-v 1) =="; ./parallel $G -m pull -t $MAXP -v 1; ./parallel $G -m push -t $MAXP -v 1
echo "== schedule/chunk sweep at p=$MAXP =="
bench sched pull_static $MAXP $C ./parallel $G -m pull -S static -t $MAXP -v 0
for c in 1 16 64 256 1024 4096; do bench sched pull_dynamic_c$c $MAXP $C ./parallel $G -m pull -S dynamic -c $c -t $MAXP -v 0; done
for c in 1 64 1024; do bench sched pull_guided_c$c $MAXP $C ./parallel $G -m pull -S guided -c $c -t $MAXP -v 0; done
bench sched push_static  $MAXP $C ./parallel $G -m push -S static -t $MAXP -v 0
bench sched push_dynamic64 $MAXP $C ./parallel $G -m push -S dynamic -c 64 -t $MAXP -v 0
python3 ../tools/analyze.py $S
echo "imbalance_R per schedule:"; grep -h -o "variant=[a-z]* sched=[a-z]* chunk=[0-9]*\|imbalance_R=[0-9.]*" results/cs4_sched.log | paste - - | sort -u
