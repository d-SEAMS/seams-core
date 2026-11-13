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

if [[ ! -f "$build/compile_commands.json" ]]; then
  echo "missing $build/compile_commands.json" >&2
  exit 1
fi

build=$(cd "$build" && pwd)
omp="$build/tidy-omp"
mkdir -p "$omp"
rm -f "$omp/omp.h"
# A short stand-in for omp.h. Clang's parser rejects the malloc attribute
# in the GCC 13+ header, and the real declarations are not what the checks
# read. The pragmas do not need this file.
cat > "$omp/omp.h" << 'EOF'
#pragma once
#ifdef __cplusplus
extern "C" {
#endif
int omp_get_max_threads(void);
int omp_get_num_threads(void);
int omp_get_thread_num(void);
void omp_set_num_threads(int);
int omp_in_parallel(void);
double omp_get_wtime(void);
void omp_set_dynamic(int);
int omp_get_num_devices(void);
void omp_pause_resource_all(int);
#ifdef __cplusplus
}
#endif
#define omp_pause_soft 1
#define omp_pause_hard 2
EOF

mapfile -t files < <(find "$root/src" -name '*.cpp' ! -path '*/include/*' | sort)
if [[ ${#files[@]} -eq 0 ]]; then
  echo "no C++ sources under src/" >&2
  exit 1
fi

"$tidy" -p "$build" --extra-arg="-I$omp" --quiet "${files[@]}"
