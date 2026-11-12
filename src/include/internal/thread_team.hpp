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

// OpenMP work that does not hand the cores to Rayon. Pausing the team
// here parks every worker and the next parallel region pays a wake-up
// that moves by several milliseconds from call to call. Widen the mask
// and freeze the team size, and leave the workers asleep only under
// the passive wait policy.
inline void prepareCompute() {
  // Do not force a passive wait. A team that sleeps between calls wakes
  // up on a spread of several milliseconds; the benchmark keeps the
  // team with OMP_WAIT_POLICY=active.
  widenAffinity();
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
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
  // The next linkcell search uses this width. OpenMP is the team that
  // walks the cutoff rows; an idle Rayon pool waits passively.
  lc_set_threads(n);
#else
  (void)n;
#endif
}

// Thread count for a call that is too small for the fork. Measured on
// this host: below 2048 atoms an eight-thread region is slower than the
// serial walk, for both the cutoff list and the CHILL+ cell walk.
inline int teamFor(int n) {
#ifdef SEAMS_HAS_OPENMP
  const int team = std::max(1, omp_get_max_threads());
  if (n < 2048) {
    return 1;
  }
  return team;
#else
  (void)n;
  return 1;
#endif
}

} // namespace seams
