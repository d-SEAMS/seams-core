#!/usr/bin/env bash
#SBATCH --job-name=seams-offload
#SBATCH --partition=gpu-2xA100
#SBATCH --account=chem-ui
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --time=02:00:00
#SBATCH --output=offload-elja-%j.out
#
# OpenMP target offload build, Catch2 suite, and nsys profile on an
# Elja A100 node. Not the host-only pixi environment. Slurm copies
# the batch script into the spool, so BASH_SOURCE is not the repo
# path; SLURM_SUBMIT_DIR is.
#
# GPU batch shells do not have Lmod on PATH. Prefixes only: nvc++
# from the OpenHPC NVIDIA HPC SDK, Meson/Ninja/Eigen/FlexiBLAS from
# the GCCcore-13.3.0 EasyBuild tree, nsys from CUDA 12.4. Cluster
# Catch2 is 2.x; the tree carries a Catch2 3 wrap.
set -euo pipefail

EB=/hpcapps/lib-edda/easybuild/software
NVHPC_ROOT=/opt/ohpc/pub/compiler/nvhpc/22.3/Linux_x86_64/22.3
MESON_PRE=$EB/Meson/1.4.0-GCCcore-13.3.0
NINJA_PRE=$EB/Ninja/1.12.1-GCCcore-13.3.0
EIGEN_PRE=$EB/Eigen/3.4.0-GCCcore-13.3.0
FLEXI_PRE=$EB/FlexiBLAS/3.4.4-GCC-13.3.0
GCCCORE=$EB/GCCcore/13.3.0
PY312=$EB/Python/3.12.3-GCCcore-13.3.0
CUDA124=$EB/CUDA/12.4.0
CLANG17=$EB/Clang/17.0.6-GCCcore-13.2.0
Z3=$EB/Z3/4.13.0-GCCcore-13.2.0

export PATH=$CUDA124/bin:$NVHPC_ROOT/compilers/bin:$MESON_PRE/bin:$NINJA_PRE/bin:$PY312/bin:$PATH
export PKG_CONFIG_PATH=$EIGEN_PRE/share/pkgconfig:$FLEXI_PRE/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
export LD_LIBRARY_PATH=$PY312/lib:$NVHPC_ROOT/compilers/lib:$NVHPC_ROOT/cuda/lib64:$NVHPC_ROOT/math_libs/lib64:$GCCCORE/lib64:$FLEXI_PRE/lib:$Z3/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH=$MESON_PRE/lib/python3.12/site-packages${PYTHONPATH:+:$PYTHONPATH}
export NVHPC=/opt/ohpc/pub/compiler/nvhpc/22.3

ROOT=${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}
JOB_ID=${SLURM_JOB_ID:-manual}
OUT=${OUT_DIR:-$ROOT/offload-elja-$JOB_ID}
BUILD=$OUT/build
mkdir -p "$OUT"
cd "$ROOT"

export SEAMS_OFFLOAD=${SEAMS_OFFLOAD:-1}

# EasyBuild meson is a Python entry point with a `python` shebang.
# GPU images have python3 from the prefix, not /usr/bin/python.
meson() {
  "$PY312/bin/python3" "$MESON_PRE/bin/meson" "$@"
}

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
    echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
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
  if [[ -f $BUILD/meson-logs/meson-log.txt ]]; then
    cp "$BUILD/meson-logs/meson-log.txt" "$OUT/meson-log-$compiler.txt" || true
  fi
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
# nvc as CC fails Meson's C sanity check. Host C is GCCcore; nvc++
# is the C++ compiler that owns -mp=gpu.
export PATH=$GCCCORE/bin:$PATH
if configure_offload nvc++ nvc++ "$GCCCORE/bin/gcc"; then
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
