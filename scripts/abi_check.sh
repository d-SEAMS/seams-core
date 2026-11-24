#!/usr/bin/env -S -u LD_PRELOAD /bin/sh
# Copyright (c) d-SEAMS developers
# Compare the exported C ABI with the committed baseline.
# The baseline was taken against main 2dbefa97. A difference is an
# unannounced break and fails the job. Pass "write" to refresh it.
#
# When SEAMS_ABIDIFF names abidiff and SEAMS_ABI_DEBUG_LIB a libseams
# built with -g, libseams is also compared type by type with the ABI 1
# corpus in abi_baseline.abi. Added functions pass; any other change
# fails. Without either, or with an abidiff too old to read the corpus,
# that step is skipped and says so. `pixi run -e abi` provides abidiff.
# The corpus is libseams at 0127f41e (2.14.0, ABI 1), written from a -g
# build by libabigail 2.9 (corpus format 2.4, which 2.9 and 3.0 read)
# with
#   abidw --headers-dir include --drop-private-types --drop-undefined-syms \
#     --no-corpus-path --no-comp-dir-path --no-show-locs \
#     --no-load-undefined-interfaces -o scripts/abi_baseline.abi libseams.so
# and the abi-instr path cut to src/seams_api.cpp.
set -eu
# The shebang drops LD_PRELOAD before this shell is mapped. Unset inside
# the script is too late: LeakSanitizer already owns the shell and reports
# its allocations at exit. nm and abidiff abort the same way, so children
# must not inherit a preload either.
unset LD_PRELOAD
build=${1:?build directory}
baseline=${2:?baseline file}
mode=${3:-check}
corpus=${baseline%/*}/abi_baseline.abi
case $0 in
  */*) _abi_dir=${0%/*} ;;
  *) _abi_dir=. ;;
esac
include=$(cd "$_abi_dir/../include" && pwd)

# The file Meson wrote. The baseline headings stay lib*.so so the text
# matches the ELF dump; only the path changes.
lib_file() {
  stem=$1
  if [ -f "${stem}.so" ]; then
    printf '%s\n' "${stem}.so"
  elif [ -f "${stem}.dylib" ]; then
    printf '%s\n' "${stem}.dylib"
  elif [ -f "${stem}.dll" ]; then
    printf '%s\n' "${stem}.dll"
  else
    printf '%s\n' "${stem}.so"
  fi
}

# Address-free, stable text. C linkage only: C++ symbols start with _Z.
# GNU nm -D is the dynamic symbol table. On Darwin, llvm-nm -D reads the
# export trie the same way; Apple nm -gU is the defined global symbols.
nm_defined() {
  lib=$1
  if [ "$(uname -s)" != Darwin ]; then
    nm -D --defined-only "$lib" | awk '{print $2, $3}' | awk '$2 !~ /^_Z/ && $2 ~ /^seams_/ {print}' | LC_ALL=C sort
    return
  fi
  nm_bin=$(command -v nm || true)
  llvm=""
  if [ -n "$nm_bin" ] && [ -x "${nm_bin%/*}/llvm-nm" ]; then
    llvm=${nm_bin%/*}/llvm-nm
  elif command -v llvm-nm >/dev/null 2>&1; then
    llvm=$(command -v llvm-nm)
  fi
  raw=$build/nm-defined.txt
  if [ -n "$llvm" ] && "$llvm" -D --defined-only "$lib" >"$raw" 2>/dev/null; then
    :
  else
    nm -gU "$lib" >"$raw"
  fi
  # Mach-O prefixes every C symbol with one underscore. llvm-nm
  # sometimes prints the C spelling already; strip only when the
  # extra underscore is in the listing.
  awk '{
    name[NR] = $3
    type[NR] = $2
    n = NR
    if ($3 ~ /^_seams_/ || $3 ~ /^__Z/) macho = 1
  }
  END {
    for (i = 1; i <= n; i++) {
      sym = name[i]
      if (macho && sym ~ /^_/) sub(/^_/, "", sym)
      if (sym ~ /^seams_/ && sym !~ /^_Z/) print type[i], sym
    }
  }' "$raw" | LC_ALL=C sort
}

dump_lib() {
  lib=$1
  title=$2
  echo "### $title"
  nm_defined "$lib"
}

{
  echo "# d-SEAMS exported C ABI"
  echo "# reference: 2dbefa97d10db3fbacea04bca935b65d21ce8a6f"
  echo "# libseams.so and libseams_lang.so are new since that commit."
  echo "# libyodaLib.so lists the seams_* C symbols still exported from the core library."
  dump_lib "$(lib_file "$build/src/libseams")" libseams.so
  dump_lib "$(lib_file "$build/bindings/libseams_lang")" libseams_lang.so
  dump_lib "$(lib_file "$build/src/libyodaLib")" libyodaLib.so
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
echo "exported symbols match $baseline"

abidiff=${SEAMS_ABIDIFF:-}
debug_lib=${SEAMS_ABI_DEBUG_LIB:-}
if [ -z "$abidiff" ] || [ ! -x "$abidiff" ]; then
  echo "abidiff not found; type-level check against ABI 1 skipped"
  exit 0
fi
if [ -z "$debug_lib" ] || [ ! -f "$debug_lib" ]; then
  echo "no libseams with debug info; type-level check against ABI 1 skipped"
  exit 0
fi
if [ ! -f "$corpus" ]; then
  echo "missing ABI 1 corpus $corpus" >&2
  exit 1
fi

# abidiff exit status is a bit set: 1 error, 2 usage, 4 change,
# 8 incompatible change. With added symbols left out, any bit fails.
status=0
"$abidiff" --no-added-syms --drop-private-types --headers-dir2 "$include" \
  "$corpus" "$debug_lib" > "$build/capi-abidiff.txt" 2>&1 || status=$?
if [ "$status" -eq 1 ] && grep -q "incompatible format version" "$build/capi-abidiff.txt"; then
  cat "$build/capi-abidiff.txt"
  echo "$abidiff cannot read the ABI 1 corpus format; type-level check skipped"
  exit 0
fi
if [ "$status" -ne 0 ]; then
  cat "$build/capi-abidiff.txt"
  echo "libseams differs from the ABI 1 corpus beyond added functions (abidiff status $status)" >&2
  exit 1
fi
status=0
"$abidiff" --drop-private-types --headers-dir2 "$include" \
  "$corpus" "$debug_lib" > "$build/capi-abidiff-added.txt" 2>&1 || status=$?
if [ $((status & 11)) -ne 0 ]; then
  cat "$build/capi-abidiff-added.txt"
  echo "abidiff failed (status $status)" >&2
  exit 1
fi
cat "$build/capi-abidiff-added.txt"
echo "libseams keeps the ABI 1 corpus; additions only"
