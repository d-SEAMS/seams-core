#!/bin/bash
# Copyright (c) d-SEAMS developers
set -euxo pipefail

# conda-forge's compiler activation sets MESON_ARGS (prefix, libdir,
# buildtype, and the cross file when cross-compiling).
if [[ -z "${MESON_ARGS:-}" ]]; then
  MESON_ARGS="-Dbuildtype=release --prefix=${PREFIX} -Dlibdir=lib"
fi

# Every optional backend that is not a conda-forge package is off, and
# nodownload makes a missing subproject an error instead of a fetch.
# chemfiles stays off until conda-forge ships one that passes the
# mmCIF probe (cmake/chemfiles_mmcif_probe.cpp); the build would
# otherwise download and patch 0.10.4.
# The linux archive script sets SEAMS_WITH_READCON=enabled after it
# has fetched the wraps. Every other build, including conda-forge and
# a checkout where the wrap directory is already present, leaves the
# reader off.
if [[ ${SEAMS_WITH_READCON:-} == enabled ]]; then
  readcon_opt=-Dwith_readcon=enabled
else
  readcon_opt=-Dwith_readcon=disabled
fi
meson setup _build_lib ${MESON_ARGS} \
  --wrap-mode=nodownload \
  -Dwith_tests=false \
  -Dwith_cli=true \
  -Dwith_blas=disabled \
  -Dwith_chemfiles=disabled \
  "$readcon_opt" \
  -Dwith_linkcell=disabled \
  -Dwith_minimage=disabled \
  -Dwith_gpulite=disabled \
  -Dwith_openmp_offload=disabled \
  -Dwith_mpi=disabled \
  -Dwith_lammps=disabled \
  -Dwith_gromacs=disabled \
  -Dwith_openmm=disabled \
  -Dwith_plumed=disabled \
  -Dwith_ira=disabled \
  -Dwith_sphericart=disabled \
  -Dwith_csm=disabled \
  -Dwith_spglib=disabled \
  -Dwith_nauty=disabled

meson compile -C _build_lib -j "${CPU_COUNT:-2}"
meson install -C _build_lib --no-rebuild
