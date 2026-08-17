#!/usr/bin/env bash
#SBATCH --job-name=seams-gpu
#SBATCH --partition=gpu-2xA100
#SBATCH --account=chem-ui
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --time=01:00:00
#SBATCH --output=repro/results/elja-gpu-%j.out
# Device-resident N-frame batch on Elja. Not Terra.
# gpu-2xA100 advertises GRES; gpu-1xA100 does not. Slurm copies the
# batch script into the spool, so BASH_SOURCE is not the repo path.
set -euo pipefail
# cargo builds the linkcell wrap. libcudart/libnvrtc live under the
# OpenHPC CUDA prefix and are not on the default linker path.
export PATH=$HOME/.pixi/bin:$HOME/.cargo/bin:/opt/ohpc/pub/compiler/cuda/12.2/bin:$PATH
CUDA_LIB=
for d in \
  /opt/ohpc/pub/compiler/cuda/12.2/targets/x86_64-linux/lib \
  /opt/ohpc/pub/compiler/cuda/12.2/lib64; do
  if [[ -e $d/libcudart.so.12 || -e $d/libcudart.so ]]; then
    CUDA_LIB=$d
    break
  fi
done
if [[ -z $CUDA_LIB ]]; then
  echo "no libcudart.so on this node" >&2
  exit 1
fi
export LD_LIBRARY_PATH=$CUDA_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
ROOT=${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}
OUT=$ROOT/repro/results
BUILD=/tmp/seams-gpu-${SLURM_JOB_ID:-manual}
mkdir -p "$OUT" "$BUILD"
: > "$OUT/tip-gpu-batch.txt"
: > "$OUT/gpu-test.log"
cd "$ROOT"
{
  echo "hostname: $(hostname)"
  echo "jobid: ${SLURM_JOB_ID:-none}"
  echo "partition: ${SLURM_JOB_PARTITION:-none}"
  echo "cuda_lib: $CUDA_LIB"
  echo "nsys: $(command -v nsys || echo missing)"
  echo "cargo: $(command -v cargo || echo missing)"
  echo "rustc: $(command -v rustc || echo missing)"
  nvidia-smi -L 2>/dev/null || true
  echo "loadavg_at_start: $(cut -d' ' -f1 /proc/loadavg)"
} | tee "$OUT/gpu-conditions.txt"

# cargo is already on PATH. Do not wrap pixi run in `env PATH=...`:
# that expansion is the pre-pixi PATH and drops meson.
pixi run -- meson setup "$BUILD" --buildtype=release -Dwith_tests=true \
  -Dwith_gpulite=enabled
pixi run -- meson compile -C "$BUILD"
pixi run -- meson test -C "$BUILD" --print-errorlogs | tee "$OUT/gpu-test.log"

cd "$ROOT/input"
export LD_LIBRARY_PATH=$BUILD/src:$CUDA_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
NSYS=$(command -v nsys || true)
if [[ -z $NSYS ]]; then
  echo "nsys not on PATH" >&2
  exit 1
fi
# One process: cold NVRTC plus warm repeats. nsys owns the timeline.
$NSYS profile --trace=cuda,nvtx,osrt --sample=none --cpuctxsw=none \
  --force-overwrite=true --stats=true \
  -o "$OUT/tip-gpu-nsys" \
  "$BUILD/tests/bench_gpu_batch" traj/mW_cubic.lammpstrj 11 1 5 \
  | tee "$OUT/tip-gpu-batch.txt"
if ! grep -q '^resident yes$' "$OUT/tip-gpu-batch.txt"; then
  echo "device batch did not reside" >&2
  exit 1
fi
$NSYS stats --force-export=true --report cuda_gpu_kern_sum --report cuda_api_sum \
  --report cuda_gpu_mem_time_sum \
  "$OUT/tip-gpu-nsys.nsys-rep" | tee "$OUT/tip-gpu-nsys-stats.txt"
echo DONE
