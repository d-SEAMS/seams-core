#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
  echo "usage: $0 BUILD-DIR INPUT-PREFIX OUTPUT-PREFIX" >&2
  exit 2
fi

build=$1
input=$2
output=$3

meson compile -C "$build" derive_confined_water_fixture
"$build/tests/derive_confined_water_fixture" \
  "$input.gro" "$input.xtc" "$output.gro" "$output.xtc"
