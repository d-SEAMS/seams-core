#!/bin/sh
# Copyright (c) d-SEAMS developers
# SPDX-License-Identifier: MIT
#
# Build the GPL-2 shape program as its own executable. It is not linked
# into libyodaLib. The C sources are K&R, so they stay on the C compiler.

set -eu
out=$1
root=$2
src=$root/third_party/csm/src
CXX=/usr/bin/g++
CC=/usr/bin/gcc
if [ ! -x "$CXX" ]; then CXX=g++; fi
if [ ! -x "$CC" ]; then CC=gcc; fi
if [ -x /usr/bin/pkg-config ]; then
  PKG=/usr/bin/pkg-config
else
  PKG=pkg-config
fi
if "$PKG" --exists openbabel-3; then
  cflags=$("$PKG" --cflags openbabel-3)
  libs=$("$PKG" --libs openbabel-3)
else
  cflags=$("$PKG" --cflags openbabel-2.0)
  libs=$("$PKG" --libs openbabel-2.0)
fi
obj=$(dirname "$out")/csm-obj
mkdir -p "$obj"
$CC -O2 -std=c99 -I"$src" -c "$src/matrix.c" -o "$obj/matrix.o"
$CC -O2 -std=c99 -I"$src" -c "$src/Molecule.c" -o "$obj/Molecule.o"
$CC -O2 -std=c99 -I"$src" -c "$src/svdcmp.c" -o "$obj/svdcmp.o"
$CC -O2 -std=c99 -I"$src" -c "$src/parseFunctions.c" -o "$obj/parseFunctions.o"
$CC -O2 -std=c99 -I"$src" -c "$src/permuter.c" -o "$obj/permuter.o"
$CXX -O2 $cflags -I"$src" -c "$src/babelAdapter.cpp" -o "$obj/babelAdapter.o"
$CXX -O2 $cflags -I"$src" -c "$src/shape.cpp" -o "$obj/shape.o"
$CXX -O2 -o "$out" "$obj/matrix.o" "$obj/Molecule.o" "$obj/svdcmp.o" \
  "$obj/parseFunctions.o" "$obj/permuter.o" "$obj/babelAdapter.o" "$obj/shape.o" \
  $libs -lm
