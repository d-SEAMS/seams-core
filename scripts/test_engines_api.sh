#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build libseams and run the C API checks and the OpenMM reporter.
# LAMMPS, GROMACS, OpenMM, and PLUMED are the meson engine tests.
set -e
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
build=$root/build
lib=$build/src
eessi_lib=/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/sapphirerapids/software/GCCcore/15.2.0/lib64

if [ ! -x "$lib/test_seams_api" ]; then
  meson setup "$build" --wrap-mode=nofallback --force-fallback-for=linkcell \
    -Dwith_tests=false \
    -Dwith_ira=disabled -Dwith_sphericart=disabled -Dwith_nauty=disabled \
    -Dwith_gpulite=disabled -Dwith_openmp_offload=disabled
  meson compile -C "$build" test_seams_api
fi
"$lib/test_seams_api"
python3 "$root/docs/build_docs.py"
python3 "$root/engines/make_ice.py"

# OpenMM reporter. SEAMS_PY is a virtualenv whose python imports openmm
# and numpy. The language bindings test themselves in their own
# repositories against a seams-core tag.
if [ -n "$SEAMS_PY" ]; then
  (cd "$root/engines/openmm" && \
    PYTHONPATH="$root/adapters/python" \
    SEAMS_LANG_LIB="$root/build/bindings/libseams_lang.so" \
    LD_LIBRARY_PATH="$eessi_lib:$lib:$root/build/bindings" \
    "$SEAMS_PY" test_reporter.py)
fi

echo "engines-api ok"
