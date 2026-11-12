#!/bin/sh
# Run the language fronts against the build tree.
# Exit 77 when a toolchain is absent, unless SEAMS_REQUIRE_LANG=1.
set -e
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
build=${1:-$root/build}
build=$(CDPATH= cd -- "$build" && pwd)
lib=$build/src
lang=$build/bindings
missing=0

note_missing() {
  echo "language front skipped: $1" >&2
  missing=1
}

if command -v python3 >/dev/null 2>&1; then
  SEAMS_LANG_LIB=$lang/libseams_lang.so \
    LD_LIBRARY_PATH=$lang:$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
    python3 "$root/tests/test_frames_abi.py"
else
  note_missing python3
fi

if command -v gfortran >/dev/null 2>&1; then
  mod=$(mktemp -d)
  gfortran -J"$mod" -c "$root/bindings/fortran/seams.f90" -o "$mod/seams.o"
  gfortran -I"$mod" "$root/bindings/fortran/test_seams.f90" "$mod/seams.o" \
    -o "$mod/test_seams" -L"$lib" -lseams -Wl,-rpath,"$lib"
  LD_LIBRARY_PATH=$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} "$mod/test_seams"
  rm -rf "$mod"
else
  note_missing gfortran
fi

if command -v cargo >/dev/null 2>&1; then
  SEAMS_LIB_DIR=$lib LD_LIBRARY_PATH=$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
    cargo test --manifest-path "$root/bindings/rust/Cargo.toml" -- --test-threads=1
else
  note_missing cargo
fi

if command -v julia >/dev/null 2>&1; then
  (cd "$root/bindings/julia/Dseams" && \
    SEAMS_LANG_LIB=$lang/libseams_lang.so \
    DSEAMS_ROOT=$root \
    LD_LIBRARY_PATH=$lang:$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
    julia --project=. test/runtests.jl)
else
  note_missing julia
fi

if [ "$missing" -ne 0 ]; then
  if [ "${SEAMS_REQUIRE_LANG:-0}" = 1 ]; then
    exit 1
  fi
  exit 77
fi
