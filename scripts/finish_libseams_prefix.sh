#!/bin/bash
# Copyright (c) d-SEAMS developers
# Make an installed libseams prefix safe to relocate.
set -euo pipefail

if [[ $# -ne 1 ]]; then
  printf '%s\n' "usage: finish_libseams_prefix.sh PREFIX" >&2
  exit 1
fi

prefix=$1
shopt -s nullglob
pc_files=("$prefix"/lib/pkgconfig/*.pc)
if [[ ${#pc_files[@]} -ne 3 ]]; then
  printf '%s\n' "expected seams.pc, seams-core.pc and seams-lang.pc" >&2
  exit 1
fi

for pc in "${pc_files[@]}"; do
  base=$(basename "$pc")
  case $base in
    seams.pc|seams-core.pc|seams-lang.pc) ;;
    *)
      printf '%s\n' "unexpected pkg-config file $base" >&2
      exit 1
      ;;
  esac
  sed -i 's|^prefix=.*|prefix=${pcfiledir}/../..|' "$pc"
  line=$(grep '^prefix=' "$pc")
  if [[ $line != 'prefix=${pcfiledir}/../..' ]]; then
    printf '%s\n' "$base still has prefix line: $line" >&2
    exit 1
  fi
  if grep -F "$prefix" "$pc"; then
    printf '%s\n' "$base still names the build prefix" >&2
    exit 1
  fi
done
