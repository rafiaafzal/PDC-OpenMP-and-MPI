# Shared helpers for the run_bench.sh scripts (source this file).
# Env knobs:  RUNS=3   MAXP=$(nproc)   MPIRUN="mpirun" (use "mpirun --oversubscribe" if you need more ranks than slots)
RUNS=${RUNS:-3}
MAXP=${MAXP:-$(nproc)}
MPIRUN=${MPIRUN:-mpirun}

# p list: 1,2,4,8,... up to MAXP (and MAXP itself if it is not a power of two)
plist() { local p=1; while [ $p -le $MAXP ]; do echo $p; p=$((p*2)); done
          if [ $((p/2)) -ne $MAXP ]; then echo $MAXP; fi; }

# bench GROUP LABEL P OUTCSV cmd...   -> runs cmd RUNS times, appends "group,label,p,run,time" rows
bench() {
  local group=$1 label=$2 p=$3 out=$4; shift 4
  [ -f "$out" ] || echo "group,label,p,run,time" > "$out"
  for r in $(seq 1 $RUNS); do
    local line; line=$("$@" 2>&1 | grep RESULT | head -1)
    local t; t=$(echo "$line" | grep -o " time=[0-9.eE+-]*" | head -1 | cut -d= -f2)
    if [ -z "$t" ]; then echo "!! no RESULT from: $*" >&2; continue; fi
    echo "$group,$label,$p,$r,$t" >> "$out"
    echo "$line" >> "${out%.csv}.log"
  done
  echo "done: $group $label p=$p" >&2
}
