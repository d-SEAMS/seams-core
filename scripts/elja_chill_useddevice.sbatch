#!/usr/bin/env bash
#SBATCH --job-name=z4io-chill
#SBATCH --partition=gpu-1xA100
#SBATCH --account=chem-ui
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
# gpu-1xA100 advertises no GRES; exclusive node is the A100.
#SBATCH --time=00:25:00
#SBATCH --output=/users/home/rog32/dseams2_repro/results/chill_useddevice_%j.out
#SBATCH --error=/users/home/rog32/dseams2_repro/results/chill_useddevice_%j.err
#
# Prior 1813073: declare target y3+mic1 still
# cuMemcpyDtoHAsync 700 (std::vector[] on target). This job: raw
# pointers only in target loops. Hypothesis: usedDevice=1, labels match.

set -euo pipefail

EB=/hpcapps/lib-edda/easybuild/software
NVHPC_23=$EB/NVHPC/23.7-CUDA-12.2.0/Linux_x86_64/23.7
GCCCORE=$EB/GCCcore/13.3.0
PY312=$EB/Python/3.12.3-GCCcore-13.3.0
CUDA124=$EB/CUDA/12.4.0
FLEXI_PRE=$EB/FlexiBLAS/3.4.4-GCC-13.3.0
HWLOC=$EB/hwloc/2.9.2-GCCcore-13.2.0
EIGEN_PRE=$EB/Eigen/3.4.0-GCCcore-13.3.0

ROOT=${SEAMS_6RAT_ROOT:-/users/home/rog32/seams-core-6rat}
BUILD=${SEAMS_6RAT_BUILD:-$ROOT/offload-elja-1803067/build}
CRT=${SEAMS_CRT:-$ROOT/elja-crt}
INC=${SEAMS_INC:-$ROOT/elja-sysroot/usr/include}
IMPLA=/users/home/rog32/seams-core-jya1-implA
OUTDIR=/users/home/rog32/dseams2_repro/results/chill_useddevice_${SLURM_JOB_ID:-manual}
mkdir -p "$OUTDIR"

export PATH=$GCCCORE/bin:$CUDA124/bin:$NVHPC_23/compilers/bin:$PY312/bin:$PATH
export CC=$GCCCORE/bin/gcc
export NVHPC=$EB/NVHPC/23.7-CUDA-12.2.0
export NVHPC_CUDA_HOME=$NVHPC_23/cuda/12.2
export LD_LIBRARY_PATH=$PY312/lib:$NVHPC_23/compilers/lib:$NVHPC_23/cuda/12.2/lib64:$GCCCORE/lib64:$FLEXI_PRE/lib:$HWLOC/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export LIBRARY_PATH=$CRT:/lib64${LIBRARY_PATH:+:$LIBRARY_PATH}
export SEAMS_OFFLOAD=1
export OMP_NUM_THREADS=1

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
  cp "$CRT"/crt1.o "$CRT"/crti.o "$CRT"/crtn.o "$CRT"/Scrt1.o "$OVL/lib64-merged/"
  mount --bind "$OVL/lib64-merged" /usr/lib64
fi

{
  echo "HYPOTHESIS: nvc++ 23.7 -mp=gpu specializedChillPlus usedDevice=1"
  echo "PRIOR: terra host 3/3 usedDevice=0"
  echo "CHANGE: compile chill_offload.cpp with SEAMS_HAS_OFFLOAD -mp=gpu"
  echo "host=$(hostname) job=${SLURM_JOB_ID:-none} date=$(date -u +%Y%m%dT%H%M%SZ)"
  nvc++ --version | head -2
  nvidia-smi -L || true
} | tee "$OUTDIR/identity.txt"

FLAGS="-I$BUILD/src -I$ROOT/src -I$ROOT/src/include/internal -I$ROOT/src/include/external \
  -I$EIGEN_PRE/include -idirafter $INC -std=c++20 -fPIC -mp -Kieee -Mnofma \
  -DSEAMS_HAS_OPENMP=1 -DSEAMS_HAS_OFFLOAD=1 -mp=gpu -gpu=nofma -gpu=nomanaged"

nvc++ $FLAGS -c "$ROOT/src/chill_offload.cpp" -o "$OUTDIR/chill_offload.o"
nvc++ $FLAGS -c "$ROOT/tests/remasure_chill_useddevice.cpp" -o "$OUTDIR/remasure.o"
nvc++ -o "$OUTDIR/remasure_chill_useddevice" "$OUTDIR/remasure.o" "$OUTDIR/chill_offload.o" \
  -L$CRT -L/lib64 -L$HWLOC/lib -L$BUILD/src \
  $BUILD/src/libyodaLib.so \
  -mp=gpu -gpu=nofma -gpu=nomanaged -lhwloc -lm

export LD_LIBRARY_PATH=$BUILD/src:$LD_LIBRARY_PATH
cd "$ROOT/input"
set +e
"$OUTDIR/remasure_chill_useddevice" | tee "$OUTDIR/chill.txt"
rc=${PIPESTATUS[0]}
set -e
echo "CHILL_EXIT=$rc" | tee -a "$OUTDIR/identity.txt"
if grep -q "usedDevice=1" "$OUTDIR/chill.txt" && [[ $rc -eq 0 ]]; then
  echo "VERIFY_z4io=PASS" | tee -a "$OUTDIR/identity.txt"
else
  echo "VERIFY_z4io=FAIL" | tee -a "$OUTDIR/identity.txt"
fi
echo DONE
exit "$rc"
