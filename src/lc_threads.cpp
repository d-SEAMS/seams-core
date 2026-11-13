//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Linkcell publishes lc_set_threads. A build of that library from before
// the symbol existed still has to link. These definitions are weak, so a
// strong copy in liblinkcell replaces them.

#include <algorithm>
#include <cstddef>
#include <cstdint>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#ifdef SEAMS_HAS_LINKCELL
#include <linkcell.h>
#endif

// linkcell 0.3.11 does not export the mesh partition used by the domain
// split. A nonzero return is the failure the callers already handle.
// A newer liblinkcell that defines these symbols replaces the weak copies.
#ifdef SEAMS_HAS_LINKCELL
extern "C" __attribute__((weak)) int
lc_halo_width(double cutoff, int images, const double *width, std::uint32_t *halo) {
  (void)cutoff;
  (void)images;
  (void)width;
  (void)halo;
  return 1;
}

extern "C" __attribute__((weak)) int
lc_mesh_cut(const double *xyz, std::size_t n, const lc_cell *box, double cutoff,
            const void *mask, double hint, std::size_t parts, int order, double skin,
            int flags, int *nxyz, double *bin, int *assign, int *point,
            std::size_t assign_cap, std::size_t point_cap, std::size_t *ncell) {
  (void)xyz;
  (void)n;
  (void)box;
  (void)cutoff;
  (void)mask;
  (void)hint;
  (void)parts;
  (void)order;
  (void)skin;
  (void)flags;
  (void)nxyz;
  (void)bin;
  (void)assign;
  (void)point;
  (void)assign_cap;
  (void)point_cap;
  (void)ncell;
  return 1;
}

extern "C" __attribute__((weak)) int
lc_partition_weights(const int *nxyz, const std::uint64_t *weight, std::size_t ncell,
                     std::size_t parts, int balance, double slack, int *assign) {
  (void)nxyz;
  (void)weight;
  (void)ncell;
  (void)parts;
  (void)balance;
  (void)slack;
  (void)assign;
  return 1;
}
#endif

extern "C" __attribute__((weak)) void lc_set_threads(int n) {
#ifdef SEAMS_HAS_OPENMP
  if (n > 0) {
    omp_set_num_threads(n);
  }
#else
  (void)n;
#endif
}

extern "C" __attribute__((weak)) int lc_num_threads(void) {
#ifdef SEAMS_HAS_OPENMP
  return std::max(1, omp_get_max_threads());
#else
  return 1;
#endif
}
