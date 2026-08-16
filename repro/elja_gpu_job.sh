#!/usr/bin/env bash
#SBATCH --job-name=seams-gpu
#SBATCH --partition=gpu-2xA100
#SBATCH --account=chem-ui
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=repro/results/elja-gpu-%j.out
# Device-resident N-frame batch on Elja. Not Terra.
set -euo pipefail
export PATH=$HOME/.pixi/bin:$PATH
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
OUT=$ROOT/repro/results
BUILD=/tmp/seams-gpu-${SLURM_JOB_ID:-manual}
mkdir -p "$OUT" "$BUILD"
cd "$ROOT"
{
  echo "hostname: $(hostname)"
  echo "jobid: ${SLURM_JOB_ID:-none}"
  echo "partition: ${SLURM_JOB_PARTITION:-none}"
  nvidia-smi -L 2>/dev/null || true
  echo "loadavg_at_start: $(cut -d' ' -f1 /proc/loadavg)"
} | tee "$OUT/gpu-conditions.txt"

pixi run -- meson setup "$BUILD" --buildtype=release -Dwith_tests=true \
  -Dwith_gpulite=enabled
pixi run -- meson compile -C "$BUILD"
pixi run -- meson test -C "$BUILD" --print-errorlogs | tee "$OUT/gpu-test.log"

cd "$ROOT/input"
export LD_LIBRARY_PATH=$BUILD/src:${LD_LIBRARY_PATH:-}
"$BUILD/tests/bench_gpu_batch" traj/mW_cubic.lammpstrj 11 1 \
  | tee "$OUT/tip-gpu-batch.txt"
echo DONE
