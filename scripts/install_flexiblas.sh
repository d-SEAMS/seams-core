#!/bin/sh
# Ubuntu 24.04 has no libflexiblas-dev. Later Ubuntu releases do. Install
# the distro package when it exists, and otherwise build FlexiBLAS 3.4.4
# with the netlib backend it ships, into /usr/local.
set -eu
if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists flexiblas; then
  exit 0
fi
apt-get update
if apt-cache show libflexiblas-dev >/dev/null 2>&1; then
  apt-get install -y libflexiblas-dev
  exit 0
fi
apt-get install -y gcc g++ gfortran cmake ninja-build pkg-config libopenblas-dev
version=3.4.4
url="https://github.com/mpimd-csc/flexiblas/releases/download/v${version}/flexiblas-${version}.tar.gz"
want=05040ae092142dd0bf38d1bb9ce33f6b475d9f9bb455e33be997932ae855c22b
src=$(mktemp -d)
trap 'rm -rf "$src"' EXIT
curl -fsSL -o "$src/flexiblas.tar.gz" "$url"
got=$(sha256sum "$src/flexiblas.tar.gz" | awk '{print $1}')
if [ "$got" != "$want" ]; then
  echo "flexiblas ${version} sha256 ${got}" >&2
  echo "expected ${want}" >&2
  exit 1
fi
tar -xzf "$src/flexiblas.tar.gz" -C "$src"
# gcc provides omp.h. The reference LAPACK test trees are not the library.
CC=gcc CXX=g++ FC=gfortran cmake -S "$src/flexiblas-${version}" -B "$src/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr/local \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DTESTS=OFF \
  -DFLEXIBLAS_DEFAULT_BLAS=OPENBLAS
cmake --build "$src/build" -j2
cmake --install "$src/build"
ldconfig
