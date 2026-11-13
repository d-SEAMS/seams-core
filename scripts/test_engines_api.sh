#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build libseams and run the C API and the language fronts.
# LAMMPS, GROMACS, OpenMM, and PLUMED are the meson engine tests.
set -e
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
build=$root/build
lib=$build/src
eessi_lib=/cvmfs/software.eessi.io/versions/2026.06/software/linux/x86_64/intel/sapphirerapids/software/GCCcore/15.2.0/lib64

if [ ! -x "$lib/test_seams_api" ]; then
  meson setup "$build" --wrap-mode=nofallback -Dwith_tests=false \
    -Dwith_ira=disabled -Dwith_sphericart=disabled -Dwith_nauty=disabled \
    -Dwith_gpulite=disabled -Dwith_openmp_offload=disabled
  meson compile -C "$build" test_seams_api
fi
"$lib/test_seams_api"
python3 "$root/docs/build_docs.py"
python3 "$root/engines/make_ice.py"

# Fortran
mod=$(mktemp -d)
gfortran -J"$mod" -c "$root/bindings/fortran/src/seams_kinds.f90" -o "$mod/seams_kinds.o"
gfortran -J"$mod" -I"$mod" -c "$root/bindings/fortran/src/seams_c.f90" -o "$mod/seams_c.o"
gfortran -J"$mod" -I"$mod" -c "$root/bindings/fortran/src/seams_support.f90" -o "$mod/seams_support.o"
gfortran -J"$mod" -I"$mod" -c "$root/bindings/fortran/src/seams.f90" -o "$mod/seams.o"
gfortran -J"$mod" -I"$mod" -c "$root/bindings/fortran/src/seams_methods.f90" -o "$mod/seams_methods.o"
gfortran -I"$mod" "$root/bindings/fortran/test_seams.f90" "$mod/seams_c.o" \
  -o "$mod/test_seams" -L"$lib" -lseams -Wl,-rpath,"$lib"
gfortran -I"$mod" "$root/bindings/fortran/test/test_high.f90" \
  "$mod/seams_kinds.o" "$mod/seams_c.o" "$mod/seams_support.o" "$mod/seams.o" "$mod/seams_methods.o" \
  -o "$mod/test_high" -L"$lib" -lseams -Wl,-rpath,"$lib"
"$mod/test_seams" "$root"
"$mod/test_high" "$root"

# Julia 1.11 ccall suite, when julia is on PATH.
if command -v julia >/dev/null 2>&1; then
  (cd "$root/bindings/julia/Dseams" && \
    SEAMS_LANG_LIB="$root/build/bindings/libseams_lang.so" \
    DSEAMS_ROOT="$root" \
    LD_LIBRARY_PATH="$root/build/bindings:$lib" \
    julia --project=. -e 'using Pkg; Pkg.instantiate(); include("test/runtests.jl")')
fi

# Rust
SEAMS_LIB_DIR=$lib LD_LIBRARY_PATH="$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  cargo test --manifest-path "$root/bindings/rust/Cargo.toml" --offline -- --test-threads=1

# Python (nanobind). SEAMS_PY is a virtualenv whose python can import numpy.
if [ -n "$SEAMS_PY" ]; then
  nb=$("$SEAMS_PY" -c 'import nanobind, pathlib; print(pathlib.Path(nanobind.__file__).parent)')
  inc=$("$SEAMS_PY" -c 'import sysconfig; print(sysconfig.get_path("include"))')
  suf=$("$SEAMS_PY" -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
  c++ -shared -fPIC -O2 -std=c++17 -fvisibility=hidden \
    -I"$root/include" -I"$build/src" -I"$nb/include" -I"$nb/ext/robin_map/include" -I"$inc" \
    "$root/bindings/python/dseams/_core.cpp" "$nb/src/nb_combined.cpp" \
    -L"$lib" -lseams -Wl,-rpath,"$lib" \
    -o "$root/bindings/python/dseams/_core$suf"
  (cd "$root/bindings/python" && PYTHONPATH="$PWD" LD_LIBRARY_PATH="$eessi_lib:$lib" \
    "$SEAMS_PY" -m unittest -v tests/test_dseams.py)
  (cd "$root/engines/openmm" && \
    PYTHONPATH="$root/bindings/python:$root/adapters/python" \
    SEAMS_LANG_LIB="$root/build/bindings/libseams_lang.so" \
    LD_LIBRARY_PATH="$eessi_lib:$lib:$root/build/bindings" \
    "$SEAMS_PY" test_reporter.py)
fi

rm -rf "$mod"
echo "engines-api ok"
