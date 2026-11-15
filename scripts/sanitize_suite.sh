#!/bin/sh
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
  -Dwith_python=true \
  -Dwith_cli=true \
  "$@"
ninja -C "$build"
# Python, Lua, and the engine hosts are not linked with the sanitizer.
# A sanitized libseams loaded into them aborts with "ASan runtime does
# not come first" unless this runtime is already mapped.
if [ -z "${LD_PRELOAD:-}" ]; then
  asan=$(${CC:-cc} -print-file-name=libasan.so 2>/dev/null || true)
  if [ -n "$asan" ] && [ "$asan" != "libasan.so" ] && [ -f "$asan" ]; then
    export LD_PRELOAD="$asan"
  fi
fi
# LeakSanitizer stays on. scripts/lsan.supp names only unsanitized
# hosts (CPython, NumPy, Qt, OVITO, OpenMM). A seams frame is not listed.
export LSAN_OPTIONS="suppressions=$root/scripts/lsan.supp${LSAN_OPTIONS:+:$LSAN_OPTIONS}"
# cargo test builds the dev profile, which is unoptimized. Linking that
# profile into this suite does not finish inside the cargo-test limit.
export CARGO_PROFILE_DEV_OPT_LEVEL=1
export CARGO_PROFILE_TEST_OPT_LEVEL=1
meson test -C "$build" --print-errorlogs
