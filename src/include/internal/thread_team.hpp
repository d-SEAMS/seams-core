#pragma once

#include <algorithm>
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
// The pool is built once and its workers sleep between searches. The
// OpenMP team is not paused: a pause makes the next region pay a
// wake-up that moves by tens of milliseconds.

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
  // Linkcell's pool sleeps on a condvar between searches. Forcing OpenMP
  // passive as well made the next region wake on a spread of tens of
  // milliseconds, which is the slow 8-thread sample. Leave OpenMP on the
  // active wait the benchmark already uses.
  widenAffinity();
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
#endif
}

inline void adoptThreads() {
  // Record the Rayon width and leave the OpenMP team running. Pausing
  // it, or forcing a passive wait, made the next region late by tens of
  // milliseconds. The pool sleeps on its own condvar between searches.
  prepareCompute();
#ifdef SEAMS_HAS_OPENMP
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

// Thread count for a call that is too small for the fork. Below 2048
// atoms an eight-thread region is slower than the serial walk. Below
// 8192 atoms four threads beat eight: the extra workers do not have a
// full share of the walk.
inline int teamFor(int n) {
#ifdef SEAMS_HAS_OPENMP
  const int team = std::max(1, omp_get_max_threads());
  if (n < 2048) {
    return 1;
  }
  if (n < 8192 && team > 4) {
    return 4;
  }
  return team;
#else
  (void)n;
  return 1;
#endif
}

// Sets the OpenMP width from teamFor for the rest of the scope, then
// puts the caller's width back. An OpenMP-only kernel uses this instead
// of adoptThreads, which would also start the Rayon pool.
struct TeamWidth {
  int saved = 1;
  bool on = false;
#ifdef __linux__
  cpu_set_t oldMask;
#endif
  explicit TeamWidth(int n) {
#ifdef SEAMS_HAS_OPENMP
    if (!omp_in_parallel()) {
      saved = std::max(1, omp_get_max_threads());
      const int want = teamFor(n);
      if (want != saved) {
        omp_set_num_threads(want);
        on = true;
      }
#ifdef __linux__
      // A four-thread team on an eight-core mask still migrates across
      // all eight. Bind it to the first cores the process is allowed,
      // before the workers are created, so they inherit that mask.
      if (want > 0 && want < saved &&
          sched_getaffinity(0, sizeof(oldMask), &oldMask) == 0) {
        int cpus[256];
        int ncpu = 0;
        for (int cpu = 0; cpu < CPU_SETSIZE && ncpu < want && ncpu < 256; ++cpu) {
          if (CPU_ISSET(cpu, &oldMask)) {
            cpus[ncpu++] = cpu;
          }
        }
        if (ncpu > 0) {
          // Once per process. Repeating the bind on a call of a few
          // hundred microseconds was the whole gap between four threads
          // requested directly and eight narrowed to four. OpenMP's own
          // bind would put the workers back on the whole allowed set,
          // so it is turned off before the team is created.
          static int bound = 0;
          if (bound != want) {
            setenv("OMP_PROC_BIND", "false", 1);
#pragma omp parallel num_threads(want)
            {
              const int tid = omp_get_thread_num();
              cpu_set_t mask;
              CPU_ZERO(&mask);
              CPU_SET(cpus[tid < ncpu ? tid : 0], &mask);
              sched_setaffinity(0, sizeof(mask), &mask);
            }
            bound = want;
          }
        }
      }
#endif
    }
#else
    (void)n;
#endif
  }
  ~TeamWidth() {
#ifdef SEAMS_HAS_OPENMP
    if (on) {
      omp_set_num_threads(saved);
    }
#endif
  }
  TeamWidth(const TeamWidth &) = delete;
  TeamWidth &operator=(const TeamWidth &) = delete;
};

} // namespace seams
