#!/bin/sh
# ThreadSanitizer suite for the OpenMP tests.
#
#   scripts/sanitize_thread.sh
#
# C++ is clang -fsanitize=thread linked to an OpenMP runtime whose
# barriers Archer reports to ThreadSanitizer. Linkcell is rebuilt with
# nightly -Zsanitizer=thread -Zbuild-std so Rayon's join is visible and
# the pool stays multithreaded. Archer asks for
# ignore_noninstrumented_modules=1 so the runtime's own atomics are not
# reported; an application race still fails the run.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
prefix=${TSAN_OPENMP_PREFIX:-$root/.tsan-openmp}
sh "$root/scripts/build_tsan_openmp.sh" "$prefix"

if ! rustup run nightly rustc --version >/dev/null 2>&1; then
  echo "nightly rustc with rust-src is required for the linkcell ThreadSanitizer build" >&2
  exit 1
fi

# Pixi exports CC as the conda gcc. The Rust sanitizer runtime is
# clang's, so this job has to compile and link with clang.
export CC=clang
export CXX=clang++
# Conda exports LDFLAGS with --as-needed. Meson places -fsanitize=thread
# before the Rust staticlib, so that flag drops the sanitizer runtime.
# The runtime is linked into the test executables, not each shared
# library, so --no-undefined cannot be required of those libraries.
if [ -n "${LDFLAGS:-}" ]; then
  LDFLAGS=$(printf '%s\n' "$LDFLAGS" | sed 's/-Wl,--as-needed//g; s/--as-needed//g')
  export LDFLAGS
fi
unset LDFLAGS_LD

# Clang does not search the conda prefix the way the conda gcc wrapper
# does. mpfr, blas and the other pixi libraries live there. The OpenMP
# runtime stays first so -fopenmp links Archer's libomp.
conda_lib=""
conda_inc=""
if [ -n "${CONDA_PREFIX:-}" ]; then
  conda_lib="$CONDA_PREFIX/lib"
  conda_inc="$CONDA_PREFIX/include"
fi
export LIBRARY_PATH="$prefix/lib${conda_lib:+:$conda_lib}${LIBRARY_PATH:+:$LIBRARY_PATH}"
export CPATH="$prefix/include${conda_inc:+:$conda_inc}${CPATH:+:$CPATH}"
export LD_LIBRARY_PATH="$prefix/lib${conda_lib:+:$conda_lib}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_TOOL_LIBRARIES="$prefix/lib/libarcher.so"
# Only the linkcell cargo invocation may see the nightly sanitizer
# flags. A process-wide RUSTFLAGS would break the other crates, which
# are built with stable rustc.
unset RUSTFLAGS
export LINKCELL_SANITIZE=thread
# Archer's own requirement. It does not hide races in instrumented code.
case ${TSAN_OPTIONS:-} in
  *ignore_noninstrumented_modules*) ;;
  *) TSAN_OPTIONS="halt_on_error=1:ignore_noninstrumented_modules=1${TSAN_OPTIONS:+:$TSAN_OPTIONS}" ;;
esac
export TSAN_OPTIONS
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"

build=$root/build-tsan
reconfigure=
if [ "${1:-}" = "--reconfigure" ]; then
  reconfigure=--reconfigure
  shift
fi
if [ -d "$build" ] && [ -z "$reconfigure" ]; then
  reconfigure=--reconfigure
fi
# Linkcell's staticlib contains a ThreadSanitizer build of std.
# Minimage's staticlib contains the ordinary std. Both define
# rust_eh_personality, and the two copies are not the same COMDAT, so
# the link needs to keep the first definition.
meson setup "$build" $reconfigure \
  --buildtype=debugoptimized \
  -Dwith_tests=true \
  -Dwith_python=false \
  -Dwith_cli=true \
  -Db_sanitize=thread \
  -Db_asneeded=false \
  -Db_lundef=false \
  -Dc_link_args='-Wl,--allow-multiple-definition' \
  -Dcpp_link_args='-Wl,--allow-multiple-definition' \
  "$@"
ninja -C "$build" \
  tests/test_neighbours \
  tests/test_bop \
  tests/test_chill_plus_rule \
  tests/test_cage_enum \
  tests/test_rings_ih
if command -v sudo >/dev/null 2>&1 && sudo -n true 2>/dev/null; then
  sudo sysctl vm.mmap_rnd_bits=28 || true
fi
meson test -C "$build" --print-errorlogs \
  neighbours bop chill_plus_rule cage_enum rings_ih
