#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build the linux-64 pixi stack and run the full Meson suite.
# A failure exits non-zero. Any skip is a failure.
set -eu

prefix=${CONDA_PREFIX:?pixi must export CONDA_PREFIX}
export PKG_CONFIG_PATH="$prefix/lib/pkgconfig:$prefix/lib.AVX2_256/pkgconfig:$prefix/share/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LIBRARY_PATH="$prefix/lib:$prefix/lib.AVX2_256${LIBRARY_PATH:+:$LIBRARY_PATH}"
export CPATH="$prefix/include${CPATH:+:$CPATH}"
export LD_LIBRARY_PATH="$prefix/lib:$prefix/lib.AVX2_256${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# conda-forge LAMMPS installs library.h only; compute dseams needs the class
# headers of the same release.
scratch=${RUNNER_TEMP:-$(mktemp -d)}
lammps_src=$(sh scripts/fetch_lammps_headers.sh "$scratch/lammps-src")
# docs_render_molecules executes from the default environment and needs a
# Python that imports bpy. Resolve it from the cycles environment explicitly.
pixi install -e cycles
cycles_python=$(
  pixi run -e cycles python - <<'PY'
import sys
import bpy
print(sys.executable)
PY
)
export MOLECULES_BPY="$cycles_python"

meson setup build --buildtype=release \
  -Dwith_tests=true -Dwith_cli=true \
  -Dwith_gromacs=enabled -Dwith_lammps=enabled \
  -Dwith_openmm=enabled -Dwith_plumed=enabled \
  -Dpaper_full=true \
  -Dlammps_src_include="$lammps_src" \
  ${SEAMS_MESON_ARGS:-}
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

skipped = []
for line in text.splitlines():
    match = re.match(r"^\s*\d+/\d+\s+(.*?)\s+SKIP\b", line)
    if match:
        skipped.append(match.group(1).rsplit(":", 1)[-1])
failed = count("Fail")
if failed or skipped:
    sys.exit(
        "suite failed={} skipped={}".format(failed, ", ".join(skipped))
    )
print("suite failed=0 skipped=0")
PY
# Skips alone leave meson at status 0. A failed test or a crash does not.
if [ "$status" -ne 0 ]; then
  exit "$status"
fi
