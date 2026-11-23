#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build LLVM's OpenMP runtime with OMPT and Archer.
#
# Archer publishes the runtime's barriers to ThreadSanitizer. The
# install prefix is the first argument, or .tsan-openmp next to this
# repository. Nothing is written under /lib, /lib64, or /usr/lib.
#
# The runtime is built from the LLVM release that matches the compiler.
# Each release adds entry points the compiler emits, such as
# __kmpc_dispatch_deinit in 20, so an older libomp leaves them undefined
# at load time. The tag is the version of the clang that compiles the
# sanitized code, whichever compiler builds the runtime itself.
# TSAN_OPENMP_REF overrides the tag.
set -eu
root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
prefix=${1:-${TSAN_OPENMP_PREFIX:-$root/.tsan-openmp}}
ref=${TSAN_OPENMP_REF:-llvmorg-$(clang -dumpversion)}
if [ -f "$prefix/lib/libomp.so" ] && [ -f "$prefix/lib/libarcher.so" ] &&
  [ "$(cat "$prefix/ref" 2>/dev/null)" = "$ref" ]; then
  echo "OpenMP ThreadSanitizer runtime $ref already at $prefix"
  exit 0
fi
src=${TSAN_OPENMP_SRC:-$root/.tsan-openmp-src}
if [ ! -f "$src/openmp/CMakeLists.txt" ] || [ "$(cat "$src/ref" 2>/dev/null)" != "$ref" ]; then
  rm -rf "$src"
  git clone --depth 1 --filter=blob:none --sparse --branch "$ref" \
    https://github.com/llvm/llvm-project.git "$src"
  printf '%s\n' "$ref" > "$src/ref"
  # Standalone openmp looks for LLVMCheckCompilerLinkerFlag in
  # ../cmake/Modules. A checkout of only openmp leaves that include
  # missing and cmake stops with an unknown command.
  git -C "$src" sparse-checkout set openmp cmake/Modules
fi
build=$src/build
rm -rf "$build"
cmake -S "$src/openmp" -B "$build" \
  -DCMAKE_C_COMPILER="${CC:-clang}" \
  -DCMAKE_CXX_COMPILER="${CXX:-clang++}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DOPENMP_ENABLE_LIBOMPTARGET=OFF \
  -DLIBOMP_OMPT_SUPPORT=ON \
  -DOPENMP_ENABLE_OMPT_TOOLS=ON \
  -DLIBOMP_ARCHER_SUPPORT=ON \
  -DLIBOMP_OMPD_SUPPORT=OFF \
  -DLIBOMP_OMPD_GDB_SUPPORT=OFF
cmake --build "$build" -j"$(nproc)"
cmake --install "$build"
test -f "$prefix/lib/libomp.so"
test -f "$prefix/lib/libarcher.so"
printf '%s\n' "$ref" > "$prefix/ref"
