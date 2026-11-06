N=nucleation.lammpstrj
# The sweep moves the molecules near the box centre of frame FRAME_A to their
# positions in frame FRAME_B, on the graph of cutoff CUTOFF angstrom over the
# atoms of type ATOM_TYPE, and times CALLS calls of each update.
FRAME_A=147
FRAME_B=148
CUTOFF=3.5
ATOM_TYPE=1
CALLS=31
taskset -c $CORE build/tests/bench_churn_sweep $N $FRAME_A $FRAME_B \
  $CUTOFF $ATOM_TYPE $CALLS > $D/churn-sweep.txt
python3 paper/tex/figures/nucleation_churn.py $N > $D/nucleation-churn.txt
# The nucleation runs read frames FIRST to LAST, and each is timed RUNS times
# after WARMUP runs.
# The clock date +%s%N counts nanoseconds.
FIRST=1
LAST=148
WARMUP=2
RUNS=7
NS_PER_US=1000
US_PER_MS=1000
for graph in cutoff seeded; do
  taskset -c $CORE build/seams cages --graph $graph --frame $FIRST \
    --last $LAST $N > $D/nucleation-$graph-carry.txt
  taskset -c $CORE build/seams cages --graph $graph --rebuild \
    --frame $FIRST --last $LAST $N > $D/nucleation-$graph-rebuild.txt
  cmp $D/nucleation-$graph-carry.txt $D/nucleation-$graph-rebuild.txt
done
{
  echo "columns graph mode rep wall_ms"
  for graph in cutoff seeded; do
    for mode in carry rebuild; do
      flag=; test $mode = rebuild && flag=--rebuild
      for rep in $(seq 1 $((WARMUP + RUNS))); do
        t0=$(date +%s%N)
        taskset -c $CORE build/seams cages --graph $graph $flag \
          --frame $FIRST --last $LAST $N > /dev/null
        t1=$(date +%s%N)
        if [ $rep -gt $WARMUP ]; then
          us=$(( (t1 - t0) / NS_PER_US ))
          echo "$graph $mode $((rep - WARMUP)) $us" | awk -v d=$US_PER_MS \
            '{printf "%s %s %s %.3f\n", $1, $2, $3, $4 / d}'
        fi
      done
    done
  done
} > $D/nucleation-timings.txt
