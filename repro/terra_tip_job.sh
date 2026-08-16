#!/usr/bin/env bash
#SBATCH --job-name=seams-repro
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=repro/results/terra-%j.out
# Tip-only campaign: pipeline, incremental, stages. No baseline worktree.
set -euo pipefail
export PATH=$HOME/.pixi/bin:$PATH
export SLURM_CONF=${SLURM_CONF:-/etc/slurm-llnl/slurm.conf}
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
OUT=$ROOT/repro/results
BUILD=/tmp/seams-repro-${SLURM_JOB_ID:-manual}
mkdir -p "$OUT" "$BUILD"
cd "$ROOT"
{
  echo "hostname: $(hostname)"
  echo "jobid: ${SLURM_JOB_ID:-none}"
  echo "partition: ${SLURM_JOB_PARTITION:-none}"
  echo "tip_note: tip-only; no baseline worktree on this host"
  lscpu | grep -E 'Model name|Socket|Core|Thread' || true
  echo "loadavg_at_start: $(cut -d' ' -f1 /proc/loadavg)"
} | tee "$OUT/conditions.txt"

pixi run -- meson setup "$BUILD" --buildtype=release -Dwith_tests=true
pixi run -- meson compile -C "$BUILD"
pixi run -- meson test -C "$BUILD" --print-errorlogs | tee "$OUT/tip-test.log"

cd "$ROOT/input"
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=$BUILD/src:${LD_LIBRARY_PATH:-}

"$BUILD/tests/bench_pipeline" 1000 2000 4000 8000 16000 32000 | tee "$OUT/tip-pipeline.txt"
"$BUILD/tests/bench_scaling" 1000 2000 4000 8000 16000 32000 | tee "$OUT/tip-scaling.txt"
for n in 1000 2000 4000 8000 16000 32000; do
  "$BUILD/tests/bench_rings_incremental" "$n" 5 3 | tee "$OUT/tip-incremental-$n.txt"
done
"$BUILD/tests/bench_stages" traj/mW_cubic.lammpstrj 1 1 5 | tee "$OUT/tip-stages-cubic.txt"
"$BUILD/tests/bench_cages" traj/mW_cubic.lammpstrj 1 5 | tee "$OUT/tip-cages.txt"
"$BUILD/tests/bench_trajectory_incremental" traj/mW_cubic.lammpstrj 11 | tee "$OUT/trajectory-incremental.txt"
if [ -f "$ROOT/repro/figshare/nucleation.lammpstrj" ]; then
  "$BUILD/tests/bench_stages" "$ROOT/repro/figshare/nucleation.lammpstrj" 1 1 5 \
    | tee "$OUT/tip-stages-nucleation.txt"
fi
echo DONE
