#!/bin/sh
# Compare the exported C ABI with the committed baseline.
# The baseline was taken against main 2dbefa97. A difference is an
# unannounced break and fails the job. Pass "write" to refresh it.
set -eu
build=${1:?build directory}
baseline=${2:?baseline file}
mode=${3:-check}

dump_lib() {
  lib=$1
  title=$2
  echo "### $title"
  # Address-free, stable text. C linkage only: C++ symbols start with _Z.
  nm -D --defined-only "$lib" | awk '{print $2, $3}' | awk '$2 !~ /^_Z/ && $2 ~ /^seams_/ {print}' | sort
}

{
  echo "# d-SEAMS exported C ABI"
  echo "# reference: 2dbefa97d10db3fbacea04bca935b65d21ce8a6f"
  echo "# libseams.so and libseams_lang.so are new since that commit."
  echo "# libyodaLib.so lists the seams_* C symbols still exported from the core library."
  dump_lib "$build/src/libseams.so" libseams.so
  dump_lib "$build/bindings/libseams_lang.so" libseams_lang.so
  dump_lib "$build/src/libyodaLib.so" libyodaLib.so
  if [ -n "${SEAMS_ABI_OLD:-}" ] && [ -f "$SEAMS_ABI_OLD" ]; then
    echo "### abidiff libyodaLib against ${SEAMS_ABI_REF:-2dbefa97}"
    abidiff --no-added-syms "$SEAMS_ABI_OLD" "$build/src/libyodaLib.so" || true
  fi
} > "$build/capi-abi.txt"

if [ "$mode" = "write" ]; then
  cp "$build/capi-abi.txt" "$baseline"
  exit 0
fi
if [ ! -f "$baseline" ]; then
  echo "missing ABI baseline $baseline" >&2
  exit 1
fi
if ! diff -u "$baseline" "$build/capi-abi.txt"; then
  echo "C ABI differs from the committed baseline (unannounced break)" >&2
  exit 1
fi
