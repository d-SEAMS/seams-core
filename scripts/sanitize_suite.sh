#!/bin/sh
# Copyright (c) d-SEAMS developers
# Re-run the suite under AddressSanitizer and UndefinedBehaviorSanitizer.
#
#   scripts/sanitize_suite.sh
#   scripts/sanitize_suite.sh --reconfigure
#
# Extra arguments are passed to `meson setup`. The build directory is
# build-san next to this repository. ThreadSanitizer is a separate
# configure: -Db_sanitize=thread.
#
# R and Julia are not part of this run. Preloading the sanitizer runtime
# makes R abort inside its own startup (readRDS of the base package), and
# Julia's JIT frees a pointer it did not obtain from the sanitizer
# allocator before any seams call. bindings/meson.build therefore does
# not register those fronts when b_sanitize contains address. The Python
# ctypes front does run, because the interpreter can host the runtime
# when it is loaded first.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
build=$root/build-san
# Resolved before meson setup so find_program('cargo') records the wrapper.
export PATH="$root/scripts/asan-path:$PATH"
reconfigure=
if [ "${1:-}" = "--reconfigure" ]; then
  reconfigure=--reconfigure
  shift
fi
if [ -d "$build" ] && [ -z "$reconfigure" ]; then
  reconfigure=--reconfigure
fi
meson setup "$build" $reconfigure \
  --buildtype=debugoptimized \
  -Db_sanitize=address,undefined \
  -Dwith_tests=true \
  -Dwith_cli=true \
  "$@"
ninja -C "$build"
# Python, Lua, and the engine hosts are not linked with the sanitizer.
# A sanitized libseams loaded into them aborts with "ASan runtime does
# not come first" unless this runtime is already mapped.
# The runtime is preloaded for meson test only. A host tool this script
# runs itself (mkdir, cat) is not built with the sanitizer, and a preloaded
# LeakSanitizer reports that tool's own allocations at its exit.
_preload=${LD_PRELOAD:-}
if [ -z "$_preload" ]; then
  asan=$(${CC:-cc} -print-file-name=libasan.so 2>/dev/null || true)
  if [ -n "$asan" ] && [ "$asan" != "libasan.so" ] && [ -f "$asan" ]; then
    _preload=$asan
  fi
fi
# The cargo wrapper reads this and hands it to the test runner only.
if [ -n "$_preload" ]; then
  export ASAN_TEST_PRELOAD=${_preload%%:*}
fi
# LeakSanitizer stays on and honors one suppressions= file. A second
# suppressions= key replaces the first. scripts/lsan.supp names unsanitized
# hosts (CPython, NumPy, Qt, OVITO, OpenMM). ci/lsan-python.supp names the
# Python and matplotlib frames capi_crosscheck leaves behind. An incoming
# suppressions= file is appended, and any other LSAN_OPTIONS keys are kept.
# A seams frame is not listed.
mkdir -p "$build"
_lsan_combined=$build/lsan-combined.supp
{
  cat "$root/scripts/lsan.supp"
  printf '\n'
  cat "$root/ci/lsan-python.supp"
  printf '\n'
} > "$_lsan_combined"
_lsan_rest=
if [ -n "${LSAN_OPTIONS:-}" ]; then
  _lsan_ifs=$IFS
  IFS=:
  for _lsan_item in $LSAN_OPTIONS; do
    case $_lsan_item in
      suppressions=*)
        _lsan_file=${_lsan_item#suppressions=}
        if [ -f "$root/$_lsan_file" ]; then
          _lsan_file=$root/$_lsan_file
        fi
        if [ ! -f "$_lsan_file" ]; then
          echo "lsan suppressions file not found: $_lsan_file" >&2
          exit 1
        fi
        # The two files above are already in the combined list.
        if [ "$_lsan_file" = "$root/scripts/lsan.supp" ] || [ "$_lsan_file" = "$root/ci/lsan-python.supp" ]; then
          continue
        fi
        cat "$_lsan_file" >> "$_lsan_combined"
        printf '\n' >> "$_lsan_combined"
        ;;
      '')
        ;;
      *)
        _lsan_rest="${_lsan_rest:+$_lsan_rest:}$_lsan_item"
        ;;
    esac
  done
  IFS=$_lsan_ifs
fi
export LSAN_OPTIONS="suppressions=$_lsan_combined${_lsan_rest:+:$_lsan_rest}"
# The dev profile is unoptimized. opt-level 1 keeps the tests inside
# the cargo-test limit once the compiler itself is not sanitized.
export CARGO_PROFILE_DEV_OPT_LEVEL=1
export CARGO_PROFILE_TEST_OPT_LEVEL=1
# Tests compile snippets with $CC. That driver must not be preloaded.
if [ -n "${CC:-}" ]; then
  export REAL_CC=$CC
  export CC=$root/scripts/asan-path/cc
fi
if [ -n "${CXX:-}" ]; then
  export REAL_CXX=$CXX
  export CXX=$root/scripts/asan-path/c++
fi
if [ -n "${FC:-}" ]; then
  export REAL_FC=$FC
  export FC=$root/scripts/asan-path/gfortran
fi
if [ -n "$_preload" ]; then
  LD_PRELOAD=$_preload meson test -C "$build" --print-errorlogs
else
  meson test -C "$build" --print-errorlogs
fi
