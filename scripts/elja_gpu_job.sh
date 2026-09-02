#!/usr/bin/env bash
#SBATCH --job-name=seams-offload
#SBATCH --partition=gpu-1xA100
#SBATCH --account=chem-ui
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=offload-elja-%j.out
#
# OpenMP target offload build, Catch2 suite, and nsys profile on an
# Elja A100 node. Not the host-only pixi environment. Slurm copies
# the batch script into the spool, so BASH_SOURCE is not the repo
# path; SLURM_SUBMIT_DIR is.
#
# Hierarchical EasyBuild modules live under
# /hpcapps/lib-edda/modules/all. GCC/13.3.0 and nvidia/nvhpc/22.3
# both claim the compiler family, so this script loads nvhpc (or
# uses clang++ from a prefix) and puts Meson, Ninja, Eigen, and
# FlexiBLAS on PATH / PKG_CONFIG_PATH from the GCCcore-13.3.0
# prefixes. Cluster Catch2 is 2.x; the tree carries a Catch2 3 wrap.
set -euo pipefail

if [ -f /etc/profile.d/lmod.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/lmod.sh
fi
module use /hpcapps/lib-edda/modules/all
module load nvidia/nvhpc/22.3

EB=/hpcapps/lib-edda/easybuild/software
MESON_PRE=$EB/Meson/1.4.0-GCCcore-13.3.0
NINJA_PRE=$EB/Ninja/1.12.1-GCCcore-13.3.0
EIGEN_PRE=$EB/Eigen/3.4.0-GCCcore-13.3.0
FLEXI_PRE=$EB/FlexiBLAS/3.4.4-GCC-13.3.0
CUDA124=$EB/CUDA/12.4.0
CLANG17=$EB/Clang/17.0.6-GCCcore-13.2.0

export PATH=$MESON_PRE/bin:$NINJA_PRE/bin:$CUDA124/bin:$PATH
export PKG_CONFIG_PATH=$EIGEN_PRE/share/pkgconfig:$FLEXI_PRE/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=$FLEXI_PRE/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

ROOT=${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}
JOB_ID=${SLURM_JOB_ID:-manual}
OUT=${OUT_DIR:-$ROOT/offload-elja-$JOB_ID}
BUILD=$OUT/build
mkdir -p "$OUT"
cd "$ROOT"

export SEAMS_OFFLOAD=${SEAMS_OFFLOAD:-1}

log_env() {
  {
    echo "hostname: $(hostname)"
    echo "jobid: $JOB_ID"
    echo "partition: ${SLURM_JOB_PARTITION:-none}"
    echo "CC: ${CC:-} ($(command -v "${CC:-nvc}" || echo missing))"
    echo "CXX: ${CXX:-} ($(command -v "${CXX:-nvc++}" || echo missing))"
    echo "nvc++: $(command -v nvc++ || echo missing)"
    echo "clang++: $(command -v clang++ || echo missing)"
    echo "nsys: $(command -v nsys || echo missing)"
    echo "meson: $(command -v meson || echo missing)"
    echo "ninja: $(command -v ninja || echo missing)"
    echo "PKG_CONFIG_PATH=$PKG_CONFIG_PATH"
    echo "modules:"
    module list 2>&1 || true
    nvidia-smi -L 2>/dev/null || true
    echo "loadavg_at_start: $(cut -d' ' -f1 /proc/loadavg)"
  } | tee "$OUT/gpu-conditions.txt"
}

if ! command -v meson >/dev/null || ! command -v ninja >/dev/null; then
  echo "meson or ninja missing after prefix setup" >&2
  exit 1
fi

configure_offload() {
  local compiler=$1
  local cxx=$2
  local cc=$3
  rm -rf "$BUILD"
  mkdir -p "$BUILD"
  export CC=$cc
  export CXX=$cxx
  echo "configuring with $compiler CC=$cc CXX=$cxx" | tee -a "$OUT/setup.log"
  if meson setup "$BUILD" --buildtype=debugoptimized \
    -Dwith_openmp_offload=enabled \
    -Dwith_python=false \
    -Dwith_lua=disabled \
    -Dwith_mpi=disabled \
    -Dwith_gpulite=disabled \
    >> "$OUT/setup.log" 2>&1; then
    meson configure "$BUILD" | tee "$OUT/meson-config.txt"
    if grep -E 'OpenMP target offload.*NO|cannot link a target region' \
      "$OUT/setup.log" "$OUT/meson-config.txt"; then
      return 1
    fi
    if grep -q 'OpenMP target offload' "$OUT/meson-config.txt" \
      || grep -q 'SEAMS_HAS_OFFLOAD' "$OUT/setup.log"; then
      return 0
    fi
    # meson configure prints the summary bool as YES
    if meson configure "$BUILD" | grep -i offload | grep -qi yes; then
      return 0
    fi
  fi
  return 1
}

: > "$OUT/setup.log"
log_env

USED_COMPILER=
if configure_offload nvc++ nvc++ nvc; then
  USED_COMPILER=nvc++
else
  echo "nvc++ offload configure failed; trying clang++ 17 libomptarget" \
    | tee -a "$OUT/setup.log"
  if [[ ! -x $CLANG17/bin/clang++ ]]; then
    echo "clang++ 17 prefix missing: $CLANG17" >&2
    exit 1
  fi
  export PATH=$CLANG17/bin:$PATH
  export LD_LIBRARY_PATH=$CLANG17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
  if configure_offload clang++ "$CLANG17/bin/clang++" "$CLANG17/bin/clang"; then
    USED_COMPILER=clang++
  else
    echo "offload probe failed for nvc++ and clang++; see $OUT/setup.log" >&2
    tail -80 "$OUT/setup.log" >&2 || true
    exit 1
  fi
fi
echo "used_compiler: $USED_COMPILER" | tee -a "$OUT/gpu-conditions.txt"

meson compile -C "$BUILD" | tee "$OUT/compile.log"
meson test -C "$BUILD" --print-errorlogs | tee "$OUT/gpu-test.log"

NSYS=$(command -v nsys || true)
if [[ -z $NSYS ]]; then
  echo "nsys not on PATH" >&2
  exit 1
fi

export LD_LIBRARY_PATH=$BUILD/src${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
cd "$ROOT/input"
$NSYS profile --trace=cuda,nvtx,osrt,openmp --sample=none --cpuctxsw=none \
  --force-overwrite=true \
  -o "$OUT/tip-gpu-nsys" \
  "$BUILD/tests/bench_strong" 4096 2 \
  | tee "$OUT/tip-gpu-batch.txt"

$NSYS stats --force-export=true --report cuda_gpu_kern_sum \
  --report cuda_api_sum --report cuda_gpu_mem_time_sum \
  "$OUT/tip-gpu-nsys.nsys-rep" | tee "$OUT/tip-gpu-nsys-stats.txt"

if [[ ! -s $OUT/tip-gpu-nsys-stats.txt ]]; then
  echo "nsys stats file is empty: $OUT/tip-gpu-nsys-stats.txt" >&2
  exit 1
fi
echo DONE
