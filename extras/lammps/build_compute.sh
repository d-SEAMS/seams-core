#!/usr/bin/env bash
# Build a LAMMPS that contains compute dseams and run extras/lammps/in.dseams.
# Usage: extras/lammps/build_compute.sh /path/to/seams-core-prefix [lammps-src]
set -euo pipefail
PREFIX="${1:?prefix that contains include/internal/seams_c_api.h and lib/libyodaLib*}"
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="${2:-$HERE/../../.lammps-src}"
if [[ ! -d "$SRC/.git" ]]; then
  git clone --depth 1 --branch stable_22Jul2025_update5 \
    https://github.com/lammps/lammps.git "$SRC"
fi
cp -f "$HERE/compute_dseams.cpp" "$HERE/compute_dseams.h" "$SRC/src/"
cmake -S "$SRC/cmake" -B "$SRC/build-dseams" \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=off \
  -D BUILD_SHARED_LIBS=off \
  -D PKG_MOLECULE=off \
  -D CMAKE_CXX_FLAGS="-I${PREFIX}/include/internal -I${PREFIX}/include" \
  -D CMAKE_EXE_LINKER_FLAGS="-L${PREFIX}/lib -L${PREFIX}/lib64 -Wl,-rpath,${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib64 -lyodaLib"
cmake --build "$SRC/build-dseams" -j"$(nproc)" --target lmp
LMP="$SRC/build-dseams/lmp"
test -x "$LMP"
cd "$HERE"
"$LMP" -in in.dseams
test -s dump.chill
# four atoms, last column is the CHILL+ integer
awk 'BEGIN{n=0} $1+0==$1 && NF>=6 {n++; if($NF+0<0 || $NF+0>8) bad=1} END{if(n<4 || bad) exit 1}' dump.chill
echo "compute dseams wrote dump.chill"
