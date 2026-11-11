#pragma once

#include <cstdlib>

// One width for OpenMP and Rayon. OpenMP already follows
// OMP_NUM_THREADS. Linkcell's pool does not, so a process whose
// hardware reports one cpu keeps the neighbour walk on one thread
// even when OpenMP was given several. lc_set_threads is read by the
// next linkcell search.

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

extern "C" void lc_set_threads(int n);
extern "C" int lc_num_threads(void);

namespace seams {

inline void adoptThreads() {
  // Idle OpenMP threads spin by default and take cores from Rayon.
  setenv("OMP_WAIT_POLICY", "PASSIVE", 0);
#ifdef SEAMS_HAS_OPENMP
  // A domain is searched by the thread that owns it. Rayon stays at one
  // so the OpenMP team is the only pool.
  if (omp_in_parallel()) {
#ifdef SEAMS_HAS_LINKCELL
    lc_set_threads(1);
#endif
    return;
  }
  int n = omp_get_max_threads();
#else
  int n = 1;
#endif
  if (n < 1) {
    n = 1;
  }
#ifdef SEAMS_HAS_LINKCELL
  lc_set_threads(n);
#else
  (void)n;
#endif
}

} // namespace seams
