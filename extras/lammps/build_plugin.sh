#!/usr/bin/env bash
# Build compute dseams as a LAMMPS plugin and run in.dseams.
# Usage: extras/lammps/build_plugin.sh LAMMPS_SRC LAMMPS_LIBDIR SEAMS_PREFIX
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
LMPSRC="${1:?lammps src dir with lammpsplugin.h}"
LMPLIB="${2:?dir containing liblammps.so}"
PREFIX="${3:?seams prefix with include/internal and lib/libyodaLib.so}"
OUT="$HERE/compute_dseams.so"
SPH=${SPHERICART_LIB:-$HOME/.local/lib/python3.12/site-packages/sphericart/lib}
CXX=${CXX:-c++}
"$CXX" -shared -fPIC -std=c++17 \
  -I"$LMPSRC" -I"$PREFIX/include/internal" -I"$PREFIX/include" \
  "$HERE/compute_dseams.cpp" "$HERE/plugin_dseams.cpp" \
  -L"$PREFIX/lib" -L"$LMPLIB" -L"$SPH" \
  -Wl,-rpath,"$PREFIX/lib" -Wl,-rpath,"$LMPLIB" -Wl,-rpath,"$SPH" \
  -lyodaLib -llammps -lsphericart -o "$OUT"
test -f "$OUT"
cd "$HERE"
LMP=${LMP:-lmp}
"$LMP" -in in.dseams
test -s dump.chill
echo "plugin wrote dump.chill"
