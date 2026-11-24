#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "thread_team.hpp"

extern "C" int lc_team_run(std::size_t nstages, const std::size_t *chunks,
                           void (*fn)(void *, std::size_t, std::size_t, std::size_t),
                           void *ctx);

// Stages of chunks in one parallel region. A run claims each chunk once
// and every chunk of a stage returns before the next stage starts.
//
// While other processes keep half of the caller's cpus busy or more, the
// region runs on linkcell's pool, the caller one of the workers: a chunk
// is claimed from a counter and a stage closes when its last claimed
// chunk returns, so neither a stage boundary nor the return waits for a
// worker the scheduler has not run. An OpenMP join or barrier waits for
// every thread of the team: beside an eight-thread spinner on eight cpus
// one join of an empty loop took 3 to 11 ms, and a q6 call of 4 ms took
// 8 to 19 ms in the processes whose team met one. On a machine the caller
// has to itself, the region is one OpenMP team with a loop per stage,
// whose spinning workers start each stage at once; Rayon workers that
// have gone to sleep do not, and the q6 call took 3.1 ms on them instead
// of 2.3.
namespace seams::team {

// Workers a run started now may use. The `worker` index a run passes is
// below it, so it sizes per-worker scratch. Linkcell's team is the pool
// run() selects with lc_set_threads(poolFor(0)), so it is never wider
// than poolFor(0); asking lc_team_width would start that pool on a
// machine where the run stays on OpenMP.
inline int width() {
#ifdef SEAMS_HAS_OPENMP
  if (omp_in_parallel()) {
    return 1;
  }
#endif
  return std::max(teamWidth(), poolFor(0));
}

// Chunks per stage for `width` workers: one per OpenMP thread on a
// machine the caller has to itself, so each thread keeps one share as the
// original loops did, and two per worker while other processes crowd the
// cpus, so a worker the scheduler parks holds up a smaller share. The
// OpenMP team can be narrower than `width`, which also covers linkcell's
// pool: eight chunks on a team narrowed to seven gave one thread two, and
// the 110592-molecule cutoff list took 2.8 ms instead of 1.5.
inline int parts(int width) {
  if (width <= 1) {
    return 1;
  }
  return othersHoldHalf() ? 2 * width : std::max(1, std::min(width, teamWidth()));
}

namespace detail {

template <class F> struct Call {
  F *fn;
  std::atomic<bool> failed{false};
  std::exception_ptr failure;
  std::mutex lock;

  void one(std::size_t stage, std::size_t chunk, std::size_t worker) noexcept {
    if (failed.load(std::memory_order_relaxed)) {
      return;
    }
    try {
      (*fn)(stage, chunk, static_cast<int>(worker));
    } catch (...) {
      const std::lock_guard<std::mutex> hold(lock);
      if (!failure) {
        failure = std::current_exception();
      }
      failed.store(true, std::memory_order_relaxed);
    }
  }

  static void trampoline(void *ctx, std::size_t stage, std::size_t chunk,
                         std::size_t worker) noexcept {
    static_cast<Call *>(ctx)->one(stage, chunk, worker);
  }
};

template <class C> void ompRun(const std::vector<std::size_t> &sizes, C &call) {
  std::size_t most = 0;
  for (const std::size_t n : sizes) {
    most = std::max(most, n);
  }
#ifdef SEAMS_HAS_OPENMP
  const int threads =
      omp_in_parallel() ? 1 : static_cast<int>(std::min<std::size_t>(
                                      static_cast<std::size_t>(omp_get_max_threads()), most));
  if (threads > 1) {
#pragma omp parallel num_threads(threads)
    {
      const std::size_t worker = static_cast<std::size_t>(omp_get_thread_num());
      for (std::size_t s = 0; s < sizes.size(); ++s) {
        const long long count = static_cast<long long>(sizes[s]);
        // Each thread takes one contiguous run of chunks, as the loops
        // this replaced did.
#pragma omp for schedule(static)
        for (long long c = 0; c < count; ++c) {
          call.one(s, static_cast<std::size_t>(c), worker);
        }
      }
    }
    return;
  }
#endif
  for (std::size_t s = 0; s < sizes.size(); ++s) {
    for (std::size_t c = 0; c < sizes[s]; ++c) {
      call.one(s, c, 0);
    }
  }
}

} // namespace detail

// Calls fn(stage, chunk, worker) once for every chunk of every stage.
// Every chunk of a stage returns before any chunk of the next starts.
// The first exception is rethrown once the region ends; chunks not yet
// started when it was thrown are skipped.
template <class F> void run(const std::vector<std::size_t> &sizes, F &&fn) {
  if (sizes.empty()) {
    return;
  }
  detail::Call<std::remove_reference_t<F>> call;
  call.fn = &fn;
  bool nested = false;
#ifdef SEAMS_HAS_OPENMP
  nested = omp_in_parallel() != 0;
#endif
  if (!nested && othersHoldHalf()) {
    lc_set_threads(poolFor(0));
    if (lc_team_run(sizes.size(), sizes.data(), &decltype(call)::trampoline, &call) != 0) {
      throw std::runtime_error("linkcell team run failed");
    }
  } else {
    detail::ompRun(sizes, call);
  }
  if (call.failure) {
    std::rethrow_exception(call.failure);
  }
}

template <class F> void run(std::initializer_list<std::size_t> chunks, F &&fn) {
  run(std::vector<std::size_t>(chunks), std::forward<F>(fn));
}

// Chunks of `grain` items covering 0..n.
inline std::size_t chunksOf(std::size_t n, std::size_t grain) {
  return grain == 0 ? 0 : (n + grain - 1) / grain;
}

// One stage over 0..n in chunks of about `grain`: fn(begin, end,
// worker). The chunk count is rounded up to a multiple of the workers
// the run starts, so each takes the same number of equal chunks.
template <class F> void forRange(std::size_t n, std::size_t grain, F &&fn) {
  if (n == 0) {
    return;
  }
  std::size_t count = chunksOf(n, std::max<std::size_t>(1, grain));
  const int most = width();
  const auto workers = static_cast<std::size_t>(
      othersHoldHalf() ? most : std::max(1, std::min(most, teamWidth())));
  if (count > 1 && workers > 1) {
    count = std::min(n, (count + workers - 1) / workers * workers);
  }
  const std::size_t step = chunksOf(n, count);
  run({chunksOf(n, step)}, [&](std::size_t, std::size_t chunk, int worker) {
    const std::size_t begin = chunk * step;
    fn(begin, std::min(n, begin + step), worker);
  });
}

} // namespace seams::team
