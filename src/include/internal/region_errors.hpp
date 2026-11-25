#pragma once

#include <atomic>
#include <exception>
#include <mutex>

namespace seams {

// An exception that leaves an OpenMP region calls std::terminate, even when
// the region runs on one thread or its if clause is false. Run each
// iteration's body through run(): the first exception is kept, iterations
// that start after it return at once, and rethrow() raises it once the
// region has closed. The skipped iterations leave their outputs unwritten,
// so a caller that catches the rethrow must not read them.
class RegionErrors {
public:
  template <class F> void run(F &&body) noexcept {
    if (failed_.load(std::memory_order_relaxed)) {
      return;
    }
    try {
      body();
    } catch (...) {
      const std::lock_guard<std::mutex> hold(lock_);
      if (!first_) {
        first_ = std::current_exception();
      }
      failed_.store(true, std::memory_order_relaxed);
    }
  }

  bool failed() const noexcept { return failed_.load(std::memory_order_relaxed); }

  // Call after the region, from the thread that opened it.
  void rethrow() const {
    if (first_) {
      std::rethrow_exception(first_);
    }
  }

private:
  std::atomic<bool> failed_{false};
  std::mutex lock_;
  std::exception_ptr first_;
};

} // namespace seams
