#pragma once

#include <cstdlib>

// One width for OpenMP and Rayon. OpenMP already follows
// OMP_NUM_THREADS. Linkcell's pool does not, so a process whose
// hardware reports one cpu keeps the neighbour walk on one thread
// even when OpenMP was given several. lc_set_threads is read by the
// next linkcell search.
//
// OMP_PROC_BIND pins the calling thread before Rayon starts its pool,
// and the workers then inherit that single cpu. The mask is widened to
// every cpu the kernel still allows, once, before the pool is created.
// A live OpenMP team is asked to stop spinning so the search can use
// those cpus.

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#ifdef __linux__
#include <sched.h>
#endif

extern "C" void lc_set_threads(int n);
extern "C" int lc_num_threads(void);

namespace seams {

inline void widenAffinity() {
#ifdef __linux__
  static int done = 0;
  if (done) {
    return;
  }
  done = 1;
  cpu_set_t trial;
  CPU_ZERO(&trial);
  if (sched_getaffinity(0, sizeof(trial), &trial) != 0) {
    return;
  }
  for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
    if (CPU_ISSET(cpu, &trial)) {
      continue;
    }
    CPU_SET(cpu, &trial);
    if (sched_setaffinity(0, sizeof(trial), &trial) != 0) {
      CPU_CLR(cpu, &trial);
    }
  }
#else
#endif
}

inline void releaseTeam() {
#ifdef SEAMS_HAS_OPENMP
  omp_pause_resource_all(omp_pause_soft);
#endif
}

inline void adoptThreads() {
  // Idle OpenMP threads spin by default and take cores from Rayon.
  setenv("OMP_WAIT_POLICY", "PASSIVE", 0);
  widenAffinity();
  releaseTeam();
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
