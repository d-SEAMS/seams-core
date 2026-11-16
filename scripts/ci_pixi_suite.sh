#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build the linux-64 pixi stack and run the full Meson suite.
# A failure exits non-zero. These optional cases may skip. Any other
# skip is a failure. Expected failures stay allowed.
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
optional = {
    "sqq_sh_triclinic",
    "docs_paper_interface_full",
    "docs_paper_cages_full",
    "docs_paper_rings_full",
    "docs_paper_brine_full",
    "docs_paper_nucleation_full",
    "docs_paper_reproduce_full",
    "docs_render_molecules",
    # libabigail's abidiff is in the abi pixi environment only. The C API
    # workflow's capi-abi job runs both there.
    "capi_abi_2_13",
    "capi_abi_2_13_selftest",
}

def count(name):
    match = re.search(rf"^{name}:\s+(\d+)\s*$", text, re.M)
    return int(match.group(1)) if match else 0

skipped = []
for line in text.splitlines():
    match = re.match(r"^\s*\d+/\d+\s+(.*?)\s+SKIP\b", line)
    if match:
        skipped.append(match.group(1).rsplit(":", 1)[-1])
unexpected = [name for name in skipped if name not in optional]
failed = count("Fail")
if failed or unexpected:
    sys.exit(
        "suite failed={} unexpected_skips={}".format(failed, ", ".join(unexpected))
    )
print("suite failed=0 skipped={} optional".format(len(skipped)))
PY
# Skips alone leave meson at status 0. A failed test or a crash does not.
if [ "$status" -ne 0 ]; then
  exit "$status"
fi
