#!/usr/bin/env bash
# Run the checked-in clang-tidy configuration on the C++ sources.
# Usage: scripts/clang_tidy_src.sh <build-dir>
# The build directory must contain compile_commands.json.
set -euo pipefail

build=${1:?build directory with compile_commands.json}
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ -n ${CLANG_TIDY:-} ]]; then
  tidy=$CLANG_TIDY
elif command -v clang-tidy >/dev/null 2>&1; then
  tidy=clang-tidy
elif command -v clang-tidy-18 >/dev/null 2>&1; then
  tidy=clang-tidy-18
else
  echo "clang-tidy is not on PATH" >&2
  exit 1
fi

if [[ -n ${CLANG:-} ]]; then
  clang=$CLANG
elif command -v clang-18 >/dev/null 2>&1; then
  clang=clang-18
elif command -v clang >/dev/null 2>&1; then
  clang=clang
else
  echo "clang is not on PATH; install clang and libomp-dev" >&2
  exit 1
fi

if [[ ! -f "$build/compile_commands.json" ]]; then
  echo "missing $build/compile_commands.json" >&2
  exit 1
fi

build=$(cd "$build" && pwd)
resource=$("$clang" -print-resource-dir)
omp="$resource/include/omp.h"
if [[ ! -f "$omp" ]]; then
  echo "missing $omp; install libomp-dev (Clang's OpenMP, not GCC's)" >&2
  exit 1
fi

mapfile -t files < <(find "$root/src" -name '*.cpp' ! -path '*/include/*' | sort)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "no C++ sources under src/" >&2
  exit 1
fi

# -I the Clang resource directory so <omp.h> is libomp's header. GCC's
# omp.h uses attributes Clang rejects, and a stub hides a missing install.
"$tidy" -p "$build" --extra-arg="-I$resource/include" --quiet "${files[@]}"
