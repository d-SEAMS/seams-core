#!/bin/sh
# Build the linux-64 pixi stack and run the full Meson suite.
# A failure or a skip exits non-zero. Expected failures stay allowed.
set -eu

prefix=${CONDA_PREFIX:?pixi must export CONDA_PREFIX}
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$prefix/lib.AVX2_256/pkgconfig:$prefix/share/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LIBRARY_PATH="$prefix/lib:$prefix/lib.AVX2_256${LIBRARY_PATH:+:$LIBRARY_PATH}"
export CPATH="$prefix/include${CPATH:+:$CPATH}"
export LD_LIBRARY_PATH="$prefix/lib:$prefix/lib.AVX2_256${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

meson setup build --buildtype=release \
  -Dwith_tests=true -Dwith_python=true -Dwith_cli=true \
  -Dwith_gromacs=enabled -Dwith_lammps=enabled \
  -Dwith_openmm=enabled -Dwith_plumed=enabled
ninja -C build -j2

set +e
meson test -C build -j2 --print-errorlogs | tee build/suite.txt
status=$?
set -e

python3 - <<'PY'
import re
import sys

text = open("build/suite.txt", encoding="utf-8").read()

def count(name):
    match = re.search(rf"^{name}:\s+(\d+)\s*$", text, re.M)
    return int(match.group(1)) if match else 0

failed, skipped = count("Fail"), count("Skipped")
if failed or skipped:
    sys.exit(f"suite failed={failed} skipped={skipped}")
print("suite failed=0 skipped=0")
PY

exit "$status"
