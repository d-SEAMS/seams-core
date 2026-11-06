# Every run reads frame FRAME of its file, and d-SEAMS makes WARMUP warm-up
# calls and REPEATS timed calls.
# The ladder frame of m has m^3 cubic cells of SITES_PER_CELL molecules, for
# each m in CELLS.
# Each SPECS entry gives the molecules, warm-up calls and timed calls of
# d-SEAMS 1 on a ladder frame, read with the atoms of type ATOM_TYPE.
FRAME=1
WARMUP=5
REPEATS=30
CELLS="3 4 5 6 7 8"
SITES_PER_CELL=8
ATOM_TYPE=1
SPECS="216:5:30 512:2:10 1000:1:6 1728:1:6"
python3 paper/tex/figures/ladder_frames.py ladder $CELLS
for m in $CELLS; do
  n=$((SITES_PER_CELL * m * m * m))
  taskset -c $CORE build/tests/bench_paper ladder/ic-$n.lammpstrj $FRAME \
    $WARMUP $REPEATS cages > $D/dseams2-ic-$n.txt 2> $D/dseams2-ic-$n.err
done
taskset -c $CORE build/tests/bench_paper input/traj/mW_cubic.lammpstrj \
  $FRAME $WARMUP $REPEATS cages > $D/dseams2-mw.txt 2> $D/dseams2-mw.err
LD_LIBRARY_PATH=$GSL_LIB bash scripts/baseline_dseams1.sh
V1=build/baseline-v1.0.1/driver-build/baseline_dseams1
for spec in $SPECS; do
  n=${spec%%:*}; rest=${spec#*:}; warm=${rest%%:*}; reps=${rest#*:}
  LD_LIBRARY_PATH=$GSL_LIB taskset -c $CORE $V1 ladder/ic-$n.lammpstrj \
    $FRAME $reps $ATOM_TYPE 6e51956f3b9059729f36bcc208f6b6d18c3fec22 $warm \
    > $D/dseams1-ic-$n.txt 2> $D/dseams1-ic-$n.err
done
