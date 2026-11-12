#!/bin/sh
# Copyright (c) d-SEAMS developers
# Time findDDC on the ideal-ice ladder.
# Release build, one thread on core 0, one warm-up, six timed calls.
# Usage: ddc_ladder.sh BENCH OUTDIR
set -eu
BENCH=$1
OUT=$2
WARMUP=1
REPEATS=6
CELLS="3 4 5 6 7 8"
mkdir -p "$OUT/frames"
python3 paper/tex/figures/ladder_frames.py "$OUT/frames" $CELLS
export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1
for m in $CELLS; do
  n=$((8 * m * m * m))
  taskset -c 0 "$BENCH" "$OUT/frames/ic-$n.lammpstrj" $WARMUP $REPEATS \
    > "$OUT/ic-$n.log" 2> "$OUT/ic-$n.err"
done
