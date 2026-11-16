#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build LAMMPS against the PLUMED already on PATH and install it under DEST.
#
# conda-forge LAMMPS is linked to PLUMED 2.9 only. A 2.10 kernel has to be
# the one this LAMMPS is compiled with. The install prefix is DEST, never
# /lib, /lib64, or /usr/lib.
#
# usage: scripts/build_lammps_plumed.sh DEST
# prints the LAMMPS src directory meson should use as -Dlammps_src_include
set -eu

dest=${1:?usage: build_lammps_plumed.sh DEST}
case $dest in
  /lib|/lib/*|/lib64|/lib64/*|/usr/lib|/usr/lib/*)
    echo "refusing to install LAMMPS under $dest" >&2
    exit 1
    ;;
esac

ver=$(pkg-config --modversion plumed)
case $ver in
  2.10.*) ;;
  *)
    echo "this script builds LAMMPS against PLUMED 2.10, pkg-config reports $ver" >&2
    exit 1
    ;;
esac

src=${LAMMPS_SRC:-$dest/src}
headers="$src/src"
tag=${LAMMPS_TAG:-stable_22Jul2025_update1}
stamp="$dest/.lammps-tag"
# A restored CI cache already has this install. Rebuild when the tag,
# the binary or the installed headers are missing. Build logs stay on stderr so the only
# stdout is the header directory.
if [ -x "$dest/bin/lmp" ] && [ -f "$headers/library.h" ] && [ -f "$dest/include/lammps/library.h" ] && [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$tag" ]; then
  printf '%s\n' "$headers"
  exit 0
fi
if [ ! -f "$src/cmake/CMakeLists.txt" ]; then
  rm -rf "$src"
  git clone --depth 1 --branch "$tag" https://github.com/lammps/lammps.git "$src" >&2
fi

# runtime mode dlopens libplumedKernel.so from PLUMED_KERNEL. MANYBODY is
# the sw pair style the mW deck uses. MPI stays off: plumed_e2e is one rank.
cmake -S "$src/cmake" -B "$dest/build" \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$dest" \
  -D CMAKE_PREFIX_PATH="${CONDA_PREFIX:-}" \
  -D CMAKE_INSTALL_RPATH="$dest/lib" \
  -D CMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -D BUILD_MPI=no \
  -D BUILD_SHARED_LIBS=yes \
  -D PKG_MANYBODY=yes \
  -D PKG_PLUMED=yes \
  -D DOWNLOAD_PLUMED=no \
  -D PLUMED_MODE=runtime >&2
cmake --build "$dest/build" --parallel "${NPROC:-$(nproc)}" >&2
cmake --install "$dest/build" >&2
test -x "$dest/bin/lmp"
printf '%s\n' "$tag" > "$stamp"
printf '%s\n' "$headers"
