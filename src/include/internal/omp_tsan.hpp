#pragma once

// libgomp's barrier does not publish a happens-before that ThreadSanitizer
// can see, so a write inside a parallel region and a read after the join
// are reported as a race. This counter is the same join, made of atomics
// the sanitizer instruments. It is empty on every other build.

#include <atomic>

#if defined(__SANITIZE_THREAD__) && defined(_OPENMP)
#define SEAMS_TSAN_JOIN 1
#include <omp.h>
#endif

namespace seams::tsan {

struct Slot {
  std::atomic<int> arrived{0};
  std::atomic<int> epoch{0};
};

inline void arrive_and_wait(Slot &slot, int threads) {
#ifdef SEAMS_TSAN_JOIN
  if (threads < 1) {
    return;
  }
  const int seen = slot.epoch.load(std::memory_order_relaxed);
  const int turn = slot.arrived.fetch_add(1, std::memory_order_acq_rel);
  if (turn + 1 == threads) {
    slot.arrived.store(0, std::memory_order_relaxed);
    slot.epoch.fetch_add(1, std::memory_order_release);
    return;
  }
  while (slot.epoch.load(std::memory_order_acquire) == seen) {
  }
#else
  (void)slot;
  (void)threads;
#endif
}

} // namespace seams::tsan
