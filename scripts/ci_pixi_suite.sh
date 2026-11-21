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

# The docs snippet and doc-block tests import dseams and pydseams from pydseamslib.
sh scripts/install_pydseamslib.sh --if-missing

# abi is solved apart from default and carries libabigail's abidiff for
# capi_abi_2_13. Supply it via a Meson native file when the caller did not
# already pass one and abidiff is not on PATH, matching CI.
abidiff_ini=""
case " ${SEAMS_MESON_ARGS:-} " in
  *" --native-file "*|*" --native-file="*) ;;
  *)
    if ! command -v abidiff >/dev/null 2>&1; then
      abidiff_bin=""
      if [ -x "$prefix/../abi/bin/abidiff" ]; then
        abidiff_bin=$(CDPATH= cd -- "$prefix/../abi/bin" && pwd -P)/abidiff
      elif command -v pixi >/dev/null 2>&1; then
        abidiff_bin=$(pixi run --frozen -e abi sh -c 'command -v abidiff' 2>/dev/null || true)
      fi
      if [ -n "$abidiff_bin" ] && [ -x "$abidiff_bin" ]; then
        abidiff_ini="$scratch/abidiff.ini"
        printf "[binaries]\nabidiff = '%s'\n" "$abidiff_bin" > "$abidiff_ini"
      fi
    fi
    ;;
esac

if [ -n "$abidiff_ini" ]; then
  set -- --native-file "$abidiff_ini"
else
  set --
fi

meson setup build --buildtype=release \
  -Dwith_tests=true -Dwith_cli=true \
  -Dwith_gromacs=enabled -Dwith_lammps=enabled \
  -Dwith_openmm=enabled -Dwith_plumed=enabled \
  -Dpaper_full=true \
  -Dlammps_src_include="$lammps_src" \
  "$@" \
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

# pin_host runs its placement checks only with eight cpus or more
# (tests/pin_host.py says why); on the four-cpu runners it reports a skip.
ALLOWED_SKIPS = {"pin_host"}

skipped = []
allowed = []
for line in text.splitlines():
    match = re.match(r"^\s*\d+/\d+\s+(.*?)\s+SKIP\b", line)
    if match:
        name = match.group(1).rsplit(":", 1)[-1]
        (allowed if name in ALLOWED_SKIPS else skipped).append(name)
failed = count("Fail")
if failed or skipped:
    sys.exit(
        "suite failed={} skipped={}".format(failed, ", ".join(skipped))
    )
print("suite failed=0 skipped=0 allowed={}".format(", ".join(allowed) or "none"))
PY
# Skips alone leave meson at status 0. A failed test or a crash does not.
if [ "$status" -ne 0 ]; then
  exit "$status"
fi
