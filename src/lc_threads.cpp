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

// linkcell 0.3.11 publishes pairs and k-nearest, not the mesh cut.
// A nonzero return sends the caller down the in-tree block path.
// A later library that defines these symbols replaces the weak copies.
struct lc_cell;

extern "C" __attribute__((weak)) int
lc_mesh_cut(const double *xyz, std::size_t n, const lc_cell *box, double cutoff,
            const void *mask, double hint, std::size_t parts, int method, double slack,
            int flags, int *nxyz, double *width, int *assign, int *point_cell,
            std::size_t assign_n, std::size_t point_n, std::size_t *ncell) {
  (void)xyz;
  (void)n;
  (void)box;
  (void)cutoff;
  (void)mask;
  (void)hint;
  (void)parts;
  (void)method;
  (void)slack;
  (void)flags;
  (void)nxyz;
  (void)width;
  (void)assign;
  (void)point_cell;
  (void)assign_n;
  (void)point_n;
  (void)ncell;
  return 1;
}

extern "C" __attribute__((weak)) int lc_halo_width(double cutoff, int images, const double *bin,
                                                   std::uint32_t *out) {
  (void)cutoff;
  (void)images;
  (void)bin;
  (void)out;
  return 1;
}

extern "C" __attribute__((weak)) int
lc_partition_weights(const int *nxyz, const std::uint64_t *weight, std::size_t ncell,
                     std::size_t parts, int method, double slack, int *assign) {
  (void)nxyz;
  (void)weight;
  (void)ncell;
  (void)parts;
  (void)method;
  (void)slack;
  (void)assign;
  return 1;
}
