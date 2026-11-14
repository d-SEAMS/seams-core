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

# Any conda interpreter: the active prefix, then python3, then python.
# A missing numpy, or a libstdc++/libpython that cannot load the built
# library, is a skip of this front. An assertion failure is not.
pick_python() {
  # An active conda env wins over the interpreter recorded at configure
  # time, so the front follows whichever conda python is on the prefix.
  if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/python" ]; then
    printf '%s\n' "${CONDA_PREFIX}/bin/python"
    return
  fi
  if [ -n "${DSEAMS_PYTHON:-}" ] && [ -x "${DSEAMS_PYTHON}" ]; then
    printf '%s\n' "${DSEAMS_PYTHON}"
    return
  fi
  command -v python3 || command -v python || true
}

# A conda interpreter is started by the host loader. A library built with
# a newer glibc then fails to open. Retry with the program interpreter
# recorded in the seams binary, which is the loader that already runs it.
build_interp() {
  seams=$build/seams
  if [ ! -x "$seams" ] || ! command -v readelf >/dev/null 2>&1; then
    return 1
  fi
  readelf -l "$seams" | awk '/Requesting program interpreter/ {print $NF}' | tr -d '[]'
}

run_python() {
  runtime=""
  if command -v g++ >/dev/null 2>&1; then
    stdcxx=$(g++ -print-file-name=libstdc++.so.6)
    case "$stdcxx" in
      /*) runtime=$(dirname "$stdcxx") ;;
    esac
  fi
  front_ld=${runtime:+$runtime:}$lang:$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
  SEAMS_LANG_LIB=$lang/libseams_lang.so \
    LD_LIBRARY_PATH=$front_ld \
    "$py" "$root/tests/test_frames_abi.py"
  status=$?
  if [ "$status" -eq 0 ] || [ "$status" -ne 77 ]; then
    return "$status"
  fi
  interp=$(build_interp) || return 77
  [ -n "$interp" ] && [ -x "$interp" ] || return 77
  libm=$(LD_LIBRARY_PATH=$front_ld ldd "$lang/libseams_lang.so" 2>/dev/null | awk '/libm\.so/{print $3; exit}')
  extra=""
  if [ -n "$libm" ]; then
    extra=$(dirname "$libm"):$(dirname "$interp"):
  fi
  SEAMS_LANG_LIB=$lang/libseams_lang.so \
    "$interp" --library-path "${extra}${front_ld}" \
    "$py" "$root/tests/test_frames_abi.py"
}

py=$(pick_python)
if [ -n "$py" ]; then
  set +e
  run_python
  py_status=$?
  set -e
  if [ "$py_status" -eq 77 ]; then
    note_missing "python ($py)"
  elif [ "$py_status" -ne 0 ]; then
    exit "$py_status"
  fi
else
  note_missing python
fi

if command -v gfortran >/dev/null 2>&1; then
  mod=$(mktemp -d)
  gfortran -J"$mod" -c "$root/bindings/fortran/seams.f90" -o "$mod/seams.o"
  gfortran -I"$mod" "$root/bindings/fortran/test_seams.f90" "$mod/seams.o" \
    -o "$mod/test_seams" -L"$lib" -lseams -Wl,-rpath,"$lib"
  LD_LIBRARY_PATH=$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} "$mod/test_seams" "$root"
  rm -rf "$mod"
else
  note_missing gfortran
fi

if command -v cargo >/dev/null 2>&1; then
  # The docs jobs also run cargo test in this tree. One lock keeps the
  # target directory from being rebuilt under a live test.
  flock "$root/bindings/rust/.cargo-test.lock" \
    env SEAMS_LIB_DIR=$lib LD_LIBRARY_PATH=$lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} \
    cargo test --manifest-path "$root/bindings/rust/Cargo.toml" -- --test-threads=1
else
  note_missing cargo
fi

if command -v julia >/dev/null 2>&1; then
  (cd "$root/bindings/julia/Dseams" && \
    julia --project=. dev/develop_jll.jl >/dev/null 2>&1 && \
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
