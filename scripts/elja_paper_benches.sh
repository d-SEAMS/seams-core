#!/usr/bin/env bash
# Exclusive-node remeasure of every d-SEAMS 2.0 paper timing row.
#
# Runs as an sbatch body on one exclusive node. Builds the tip tree, gates on
# a green test suite, builds the baseline tree in a git worktree, compiles the
# two baseline drivers directly against the baseline library, then runs the
# repo benches with pinned single-thread OpenMP (plus the named 1/2/4 sweep)
# and assembles one JSON object next to the logs.
#
# Environment:
#   SEAMS_ROOT  repo checkout to run from (default: directory of this script's
#               parent)
#   BASE_SHA    baseline commit (default: merge-base with main is NOT computed
#               here; the audited baseline is pinned)
#   OUT_DIR     output directory (default: $SEAMS_ROOT/bench-elja-$SLURM_JOB_ID)
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
SEAMS_ROOT=${SEAMS_ROOT:-$(dirname "$SCRIPT_DIR")}
BASE_SHA=${BASE_SHA:-1c8efa382c4ae88ebc68e40a991260c4d57edbd5}
JOB_ID=${SLURM_JOB_ID:-nojob}
OUT_DIR=${OUT_DIR:-$SEAMS_ROOT/bench-elja-$JOB_ID}
BASE_TREE=$SEAMS_ROOT/../seams-base-$JOB_ID

mkdir -p "$OUT_DIR"
cd "$SEAMS_ROOT"

log() { printf '%s %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$OUT_DIR/driver.log"; }

# ---------------------------------------------------------------- conditions
TIP_SHA=$(git rev-parse HEAD)
log "tip $TIP_SHA base $BASE_SHA job $JOB_ID node $(hostname)"
{
  echo "hostname: $(hostname)"
  echo "jobid: $JOB_ID"
  echo "partition: ${SLURM_JOB_PARTITION:-unknown}"
  echo "tip_sha: $TIP_SHA"
  echo "base_sha: $BASE_SHA"
  lscpu | grep -E 'Model name|Socket|Core|Thread' || true
  echo "loadavg_at_start: $(cut -d' ' -f1 /proc/loadavg)"
} > "$OUT_DIR/conditions.txt"

# ------------------------------------------------------------------ tip build
export PATH=$HOME/.pixi/bin:$PATH
pixi install >> "$OUT_DIR/driver.log" 2>&1
run() { pixi run -- "$@"; }

log "configuring tip"
run meson setup build-tip --buildtype=release \
  -Dwith_python=false -Dwith_tests=true \
  > "$OUT_DIR/tip-setup.log" 2>&1
run meson compile -C build-tip >> "$OUT_DIR/tip-setup.log" 2>&1
run meson configure build-tip > "$OUT_DIR/tip-config.txt" 2>&1 || true
grep -A40 'Optional' "$OUT_DIR/tip-setup.log" | tail -30 \
  > "$OUT_DIR/tip-backends.txt" || true

# Identity gate: a red suite aborts; it does not emit a paper table
log "identity gate: meson test on tip"
if ! run meson test -C build-tip > "$OUT_DIR/tip-test.log" 2>&1; then
  log "TIP SUITE RED -- aborting without a table"
  tail -30 "$OUT_DIR/tip-test.log"
  exit 1
fi
log "tip suite green"

# ----------------------------------------------------------------- base build
log "configuring base at $BASE_SHA"
git worktree add --detach "$BASE_TREE" "$BASE_SHA" >> "$OUT_DIR/driver.log" 2>&1
(
  cd "$BASE_TREE"
  run meson setup build-base --buildtype=release \
    -Dwith_python=false -Dwith_tests=false \
    > "$OUT_DIR/base-setup.log" 2>&1
  run meson compile -C build-base >> "$OUT_DIR/base-setup.log" 2>&1
)

# Baseline drivers: compiled straight against the baseline library, since the
# baseline predates the bench targets
EIGEN_INC=$(run bash -c 'echo $CONDA_PREFIX/include/eigen3')
BASE_INC="-I$BASE_TREE/src/include/internal -I$BASE_TREE/src/include/external -I$EIGEN_INC"
BASE_LINK="-L$BASE_TREE/build-base/src -lyodaLib -Wl,-rpath,$BASE_TREE/build-base/src -Wl,--enable-new-dtags -Wl,--allow-shlib-undefined"
log "compiling baseline drivers"
run g++ -std=c++17 -O2 -o "$OUT_DIR/base_bench_scaling" \
  tests/bench_scaling.cpp $BASE_INC $BASE_LINK \
  >> "$OUT_DIR/driver.log" 2>&1
run g++ -std=c++17 -O2 -o "$OUT_DIR/base_bench_cages" \
  tests/bench_cages_base.cpp $BASE_INC $BASE_LINK \
  >> "$OUT_DIR/driver.log" 2>&1

# ------------------------------------------------------------------ run gates
export OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores
LOAD=$(cut -d' ' -f1 /proc/loadavg)
echo "loadavg_before_timing: $LOAD" >> "$OUT_DIR/conditions.txt"
if python3 -c "import sys; sys.exit(0 if float('$LOAD') < 2.0 else 1)"; then
  log "load $LOAD ok"
else
  log "LOAD $LOAD >= 2.0 -- refusing to emit rows"
  exit 1
fi

# ---------------------------------------------------------------- tip benches
cd "$SEAMS_ROOT/input"
TIPB=$SEAMS_ROOT/build-tip/tests
log "tip bench_scaling"
"$TIPB/bench_scaling" 1000 2000 4000 8000 16000 32000 \
  > "$OUT_DIR/tip-scaling.txt" 2>&1
log "tip bench_cages"
"$TIPB/bench_cages" traj/mW_cubic.lammpstrj 1 5 \
  > "$OUT_DIR/tip-cages.txt" 2>&1
log "tip bench_overhead"
"$TIPB/bench_overhead" 4096 16000 > "$OUT_DIR/tip-overhead.txt" 2>&1
log "tip bench_strong sweep"
for t in 1 2 4; do
  OMP_NUM_THREADS=$t "$TIPB/bench_strong" 65536 \
    > "$OUT_DIR/tip-strong-t$t.txt" 2>&1
done

# --------------------------------------------------------------- base benches
log "base bench_scaling"
"$OUT_DIR/base_bench_scaling" 1000 2000 4000 8000 16000 32000 \
  > "$OUT_DIR/base-scaling.txt" 2>&1
log "base bench_cages"
"$OUT_DIR/base_bench_cages" traj/mW_cubic.lammpstrj 1 5 \
  > "$OUT_DIR/base-cages.txt" 2>&1

# ------------------------------------------------------------ identity checks
for f in tip-cages base-cages; do
  if ! grep -q 'six-rings  8192' "$OUT_DIR/$f.txt"; then
    log "WARNING: $f does not report 8192 six-rings"
  fi
done

echo "loadavg_at_end: $(cut -d' ' -f1 /proc/loadavg)" >> "$OUT_DIR/conditions.txt"
log "assembling bench.json"
python3 "$SCRIPT_DIR/elja_bench_json.py" "$OUT_DIR" > "$OUT_DIR/bench.json"
log "done: $OUT_DIR/bench.json"
