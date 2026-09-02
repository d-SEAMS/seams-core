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
# from EasyBuild NVHPC 23.7 (CUDA 12.2, matches the A100 driver),
# Meson/Ninja/Eigen/FlexiBLAS from the GCCcore-13.3.0 EasyBuild
# tree, nsys from CUDA 12.4. OpenHPC nvc++ 22.3 is the fallback
# (force -gpu=cuda11.6). Cluster Catch2 is 2.x; the tree carries
# a Catch2 3 wrap.
set -euo pipefail

EB=/hpcapps/lib-edda/easybuild/software
NVHPC_23=$EB/NVHPC/23.7-CUDA-12.2.0/Linux_x86_64/23.7
NVHPC_22=/opt/ohpc/pub/compiler/nvhpc/22.3/Linux_x86_64/22.3
NVHPC_ROOT=$NVHPC_23
MESON_PRE=$EB/Meson/1.4.0-GCCcore-13.3.0
NINJA_PRE=$EB/Ninja/1.12.1-GCCcore-13.3.0
EIGEN_PRE=$EB/Eigen/3.4.0-GCCcore-13.3.0
FLEXI_PRE=$EB/FlexiBLAS/3.4.4-GCC-13.3.0
GCCCORE=$EB/GCCcore/13.3.0
PY312=$EB/Python/3.12.3-GCCcore-13.3.0
CUDA124=$EB/CUDA/12.4.0
CLANG17=$EB/Clang/17.0.6-GCCcore-13.2.0
Z3=$EB/Z3/4.13.0-GCCcore-13.2.0
HWLOC=$EB/hwloc/2.9.2-GCCcore-13.2.0

export PATH=$CUDA124/bin:$NVHPC_ROOT/compilers/bin:$MESON_PRE/bin:$NINJA_PRE/bin:$PY312/bin:$PATH
export PKG_CONFIG_PATH=$EIGEN_PRE/share/pkgconfig:$FLEXI_PRE/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}
export PYTHONPATH=$MESON_PRE/lib/python3.12/site-packages${PYTHONPATH:+:$PYTHONPATH}
export NVHPC=$EB/NVHPC/23.7-CUDA-12.2.0
export NVHPC_CUDA_HOME=$NVHPC_ROOT/cuda/12.2

ROOT=${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}
JOB_ID=${SLURM_JOB_ID:-manual}
OUT=${OUT_DIR:-$ROOT/offload-elja-$JOB_ID}
BUILD=$OUT/build
mkdir -p "$OUT"
cd "$ROOT"
INC=${SEAMS_INC:-$ROOT/elja-sysroot/usr/include}
GCC_LIB=$GCCCORE/lib/gcc/x86_64-pc-linux-gnu/13.3.0

# Keep the vendor CUDA include list. OpenHPC 22.3 still names the
# system GCC 8 tree; rewrite only that. Always append login glibc
# headers.
write_nvc_localrc() {
  local sys=$1/compilers/bin/localrc
  if [[ ! -f $sys ]]; then
    unset NVLOCALRC
    return 0
  fi
  sed \
    -e "s|set GCCDIR=/usr/lib/gcc/x86_64-redhat-linux/8;|set GCCDIR=$GCC_LIB;|" \
    -e "s|set G77DIR=/usr/lib/gcc/x86_64-redhat-linux/8/;|set G77DIR=$GCC_LIB/;|" \
    -e "s|set GCCINC=\(.*\);|set GCCINC=\1 $INC;|" \
    -e "s|set GPPDIR=\(.*\);|set GPPDIR=\1 $INC;|" \
    "$sys" > "$OUT/nvc.localrc"
  export NVLOCALRC=$OUT/nvc.localrc
}

nvc_lib_path() {
  local root=$1
  local cuda_ver=$2
  local p=$root/compilers/lib
  if [[ -d $root/cuda/$cuda_ver/lib64 ]]; then
    p=$p:$root/cuda/$cuda_ver/lib64
  elif [[ -d $root/cuda/lib64 ]]; then
    p=$p:$root/cuda/lib64
  fi
  if [[ -d $root/math_libs/$cuda_ver/lib64 ]]; then
    p=$p:$root/math_libs/$cuda_ver/lib64
  elif [[ -d $root/math_libs/lib64 ]]; then
    p=$p:$root/math_libs/lib64
  fi
  echo "$p"
}

apply_nvhpc() {
  local root=$1
  local cuda_ver=$2
  NVHPC_ROOT=$root
  export PATH=$CUDA124/bin:$root/compilers/bin:$MESON_PRE/bin:$NINJA_PRE/bin:$PY312/bin:$PATH
  export LD_LIBRARY_PATH=$PY312/lib:$(nvc_lib_path "$root" "$cuda_ver"):$GCCCORE/lib64:$FLEXI_PRE/lib:$Z3/lib:$HWLOC/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
  if [[ -d $root/cuda/$cuda_ver ]]; then
    export NVHPC_CUDA_HOME=$root/cuda/$cuda_ver
  elif [[ -d $root/cuda ]]; then
    export NVHPC_CUDA_HOME=$root/cuda
  fi
  write_nvc_localrc "$root"
  # OpenHPC 22.3 needs GCCcore 13.3 stddef.h. nvc++ 23.7 already
  # uses GCCcore 12.3; putting 13.3 xmmintrin.h first ICEs cpp1.
  if [[ $cuda_ver == 11.6 ]]; then
    export CFLAGS="${CFLAGS:-} -I${GCC_LIB}/include"
    export CXXFLAGS="${CXXFLAGS:-} -I${GCC_LIB}/include"
  fi
}

apply_nvhpc "$NVHPC_23" 12.2

# GPU nodes have no glibc-devel. Login-node C runtime objects live
# in elja-crt/ (copied before submit). gcc and clang find crt1.o
# via LIBRARY_PATH.
CRT=${SEAMS_CRT:-$ROOT/elja-crt}
if [[ ! -f $CRT/crt1.o || ! -f $CRT/crti.o || ! -f $CRT/crtn.o ]]; then
  echo "missing C runtime objects in $CRT (copy crt1.o crti.o crtn.o Scrt1.o from a login /usr/lib64)" >&2
  exit 1
fi
# GPU images ship libc.so.6 but not the glibc-devel linker scripts.
# Write those scripts here so -lc resolves to /lib64/libc.so.6.
if [[ ! -f $CRT/libc.so ]]; then
  cat > "$CRT/libc.so" <<EOF
OUTPUT_FORMAT(elf64-x86-64)
GROUP ( /lib64/libc.so.6 $CRT/libc_nonshared.a AS_NEEDED ( /lib64/ld-linux-x86-64.so.2 ) )
EOF
fi
if [[ ! -f $CRT/libm.so ]]; then
  cat > "$CRT/libm.so" <<'EOF'
OUTPUT_FORMAT(elf64-x86-64)
GROUP ( /lib64/libm.so.6 )
EOF
fi
for pair in libpthread.so:libpthread.so.0 libdl.so:libdl.so.2 librt.so:librt.so.1; do
  name=${pair%%:*}
  target=${pair##*:}
  if [[ ! -e $CRT/$name && -e /lib64/$target ]]; then
    ln -sfn "/lib64/$target" "$CRT/$name"
  fi
done
export LIBRARY_PATH=$CRT:/lib64${LIBRARY_PATH:+:$LIBRARY_PATH}
# GPU nodes have no glibc-devel headers. Login /usr/include is
# copied to elja-sysroot/usr/include before submit.
INC=${SEAMS_INC:-$ROOT/elja-sysroot/usr/include}
if [[ ! -f $INC/features.h ]]; then
  echo "missing $INC/features.h (rsync /usr/include from a login node)" >&2
  exit 1
fi
# -idirafter keeps glibc after libstdc++ so include_next <stdlib.h>
# in cstdlib still finds the C header. CPATH would be searched first
# and include_next would then skip it.
export CFLAGS="${CFLAGS:-} -idirafter ${INC}"
export CXXFLAGS="${CXXFLAGS:-} -idirafter ${INC}"
# clang's driver does not search LIBRARY_PATH for startfiles.
# gcc accepts -B; nvc++ rejects it. Keep -B on CFLAGS only.
export CFLAGS="${CFLAGS:-} -B${CRT}"
# --start-group must appear before the Catch2 static archives so ld
# can rescan them for Catch::Session. tests/meson.build closes the
# group for nvc++.
export LDFLAGS="-Wl,--start-group ${LDFLAGS:-} -L${CRT} -L/lib64 -L${HWLOC}/lib -lhwloc"

export SEAMS_OFFLOAD=${SEAMS_OFFLOAD:-1}

# nvc++ hardcodes /usr/lib64/crt1.o. GPU nodes have no glibc-devel.
# Re-enter under a user-namespace overlay that adds the login crt
# objects and headers without hiding the runtime libc.
if [[ ${SEAMS_IN_SYSROOT:-0} != 1 ]]; then
  export SEAMS_IN_SYSROOT=1
  OVL=/tmp/seams-ovl-${SLURM_JOB_ID:-$$}
  mkdir -p "$OVL/lib64-upper" "$OVL/lib64-work" "$OVL/lib64-merged"
  export SEAMS_OVL=$OVL
  exec unshare --user --map-root-user --mount bash "$0" "$@"
fi

OVL=${SEAMS_OVL:-/tmp/seams-ovl-${SLURM_JOB_ID:-$$}}
if ! grep -q ' /usr/lib64 ' /proc/mounts; then
  mount -t overlay overlay \
    -o "lowerdir=/usr/lib64,upperdir=$OVL/lib64-upper,workdir=$OVL/lib64-work" \
    "$OVL/lib64-merged"
  cp "$CRT"/crt1.o "$CRT"/crti.o "$CRT"/crtn.o "$CRT"/Scrt1.o \
    "$OVL/lib64-merged/"
  if [[ -f $CRT/libc.so ]]; then
    cp "$CRT"/libc.so "$CRT"/libm.so "$OVL/lib64-merged/" 2>/dev/null || true
  fi
  mount --bind "$OVL/lib64-merged" /usr/lib64
fi

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
  local extra=()
  if [[ $compiler == nvc++ ]]; then
    # Meson 1.4 does not list C++20 for these nvc++ releases.
    # The tree uses concepts (CoordinateScalar); nvc++ 23.7 accepts
    # -std=c++20.
    extra+=(-Dcpp_std=none)
    export CXXFLAGS="${CXXFLAGS:-} -std=c++20"
  fi
  if meson setup "$BUILD" --buildtype=debugoptimized \
    -Dwith_openmp_offload=enabled \
    -Dwith_python=false \
    -Dwith_lua=disabled \
    -Dwith_mpi=disabled \
    -Dwith_gpulite=disabled \
    -Dcatch2:tests=false \
    "${extra[@]}" \
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
# nvc++ 23.7 ships CUDA 12.2, which matches the A100 driver.
# OpenHPC 22.3 ships only CUDA 11.6; the meson probe then needs
# -gpu=cuda11.6. Meson does not list C++20 for these nvc++
# releases, so setup uses -Dcpp_std=none -std=c++20. Clang 17
# has libomptarget bitcode but no NVPTX codegen target.
export PATH=$GCCCORE/bin:$PATH
if configure_offload nvc++ nvc++ "$GCCCORE/bin/gcc"; then
  USED_COMPILER=nvc++-23.7
else
  echo "nvc++ 23.7 offload configure failed; trying nvc++ 22.3" \
    | tee -a "$OUT/setup.log"
  apply_nvhpc "$NVHPC_22" 11.6
  export NVHPC=/opt/ohpc/pub/compiler/nvhpc/22.3
  if configure_offload nvc++ nvc++ "$GCCCORE/bin/gcc"; then
    USED_COMPILER=nvc++-22.3
  else
    echo "nvc++ offload configure failed; trying clang++ 17" \
      | tee -a "$OUT/setup.log"
    export PATH=$CLANG17/bin:$PATH
    export LD_LIBRARY_PATH=$CLANG17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
    if [[ -x $CLANG17/bin/clang++ ]] && \
       configure_offload clang++ "$CLANG17/bin/clang++" "$CLANG17/bin/clang"; then
      USED_COMPILER=clang++
    else
      echo "offload probe failed for nvc++ 23.7, nvc++ 22.3 and clang++; see $OUT/setup.log" >&2
      tail -80 "$OUT/setup.log" >&2 || true
      exit 1
    fi
  fi
fi
echo "used_compiler: $USED_COMPILER" | tee -a "$OUT/gpu-conditions.txt"
echo "NVHPC_ROOT=$NVHPC_ROOT NVHPC_CUDA_HOME=${NVHPC_CUDA_HOME:-}" \
  | tee -a "$OUT/gpu-conditions.txt"

if [[ ${PROBE_ONLY:-0} == 1 ]]; then
  echo PROBE_OK compiler=$USED_COMPILER
  exit 0
fi

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
