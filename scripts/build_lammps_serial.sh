#!/bin/sh
# Copyright (c) d-SEAMS developers
# Build a serial LAMMPS (the src/STUBS MPI) and install it under DEST.
#
# conda-forge LAMMPS on linux-64 is linked to Open MPI. Open MPI 5 still
# holds component and PMIx allocations after MPI_Finalize, and a
# sanitized run_hook that initializes it fails LeakSanitizer on those.
# This LAMMPS never loads an MPI library. The install prefix is DEST,
# never /lib, /lib64, or /usr/lib.
#
# usage: scripts/build_lammps_serial.sh DEST
# prints DEST/lib/pkgconfig, the directory of its liblammps.pc
set -eu

dest=${1:?usage: build_lammps_serial.sh DEST}
case $dest in
  /lib|/lib/*|/lib64|/lib64/*|/usr/lib|/usr/lib/*)
    echo "refusing to install LAMMPS under $dest" >&2
    exit 1
    ;;
esac

# The release and OpenMP setting are conda-forge's lammps 2024.08.29, so
# run_hook prints what docs/videos/lammps/recorded.txt holds. LAMMPS
# trims error paths at the first "src/", so the clone is not named src.
src=${LAMMPS_SRC:-$dest/lammps}
tag=${LAMMPS_TAG:-stable_29Aug2024}
stamp="$dest/.lammps-tag"
# A restored CI cache already has this install. Build logs stay on
# stderr so the only stdout is the pkg-config directory.
if [ -f "$dest/lib/pkgconfig/liblammps.pc" ] && [ -f "$dest/include/lammps/library.h" ] && [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$tag" ]; then
  printf '%s\n' "$dest/lib/pkgconfig"
  exit 0
fi
if [ -z "${LAMMPS_SRC:-}" ] && { [ ! -f "$src/cmake/CMakeLists.txt" ] || [ "$(git -C "$src" describe --tags 2>/dev/null)" != "$tag" ]; }; then
  rm -rf "$src" "$dest/build"
  git clone --depth 1 --branch "$tag" https://github.com/lammps/lammps.git "$src" >&2
fi

# The core styles are enough for engines/lammps/in.dseams (lj/cut and
# fix external).
cmake -S "$src/cmake" -B "$dest/build" \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$dest" \
  -D CMAKE_INSTALL_LIBDIR=lib \
  -D CMAKE_INSTALL_RPATH="$dest/lib" \
  -D CMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -D BUILD_MPI=no \
  -D BUILD_OMP=yes \
  -D BUILD_SHARED_LIBS=yes >&2
cmake --build "$dest/build" --parallel "${NPROC:-$(nproc)}" >&2
cmake --install "$dest/build" >&2
test -f "$dest/lib/pkgconfig/liblammps.pc"
printf '%s\n' "$tag" > "$stamp"
printf '%s\n' "$dest/lib/pkgconfig"
