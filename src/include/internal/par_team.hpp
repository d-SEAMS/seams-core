#ifndef SEAMS_PAR_TEAM_H_
#define SEAMS_PAR_TEAM_H_

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace seams {

inline int teamThreads() {
#ifdef SEAMS_HAS_OPENMP
  if (omp_in_parallel()) {
    return 1;
  }
  const int n = omp_get_max_threads();
  return n > 0 ? n : 1;
#else
  return 1;
#endif
}

// In-place exclusive prefix: values[i] becomes the sum of the earlier entries.
template <class T> void exclusiveScan(std::vector<T> &values) {
  const int n = static_cast<int>(values.size());
  const int threads = teamThreads();
  if (threads <= 1 || n < threads * 2) {
    T sum = 0;
    for (int i = 0; i < n; i++) {
      const T current = values[static_cast<std::size_t>(i)];
      values[static_cast<std::size_t>(i)] = sum;
      sum += current;
    }
    return;
  }
  const int chunk = (n + threads - 1) / threads;
  std::vector<T> block(static_cast<std::size_t>(threads), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel
#endif
  {
#ifdef SEAMS_HAS_OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int lo = tid * chunk;
    const int hi = std::min(n, lo + chunk);
    T sum = 0;
    for (int i = lo; i < hi; i++) {
      sum += values[static_cast<std::size_t>(i)];
    }
    if (tid < threads) {
      block[static_cast<std::size_t>(tid)] = sum;
    }
  }
  T carried = 0;
  for (int tid = 0; tid < threads; tid++) {
    const T sum = block[static_cast<std::size_t>(tid)];
    block[static_cast<std::size_t>(tid)] = carried;
    carried += sum;
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel
#endif
  {
#ifdef SEAMS_HAS_OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int lo = tid * chunk;
    const int hi = std::min(n, lo + chunk);
    T sum = block[static_cast<std::size_t>(tid)];
    for (int i = lo; i < hi; i++) {
      const T current = values[static_cast<std::size_t>(i)];
      values[static_cast<std::size_t>(i)] = sum;
      sum += current;
    }
  }
}

// Stable merge sort. Ties keep the earlier element, so the order matches
// one serial stable sort at every thread count.
template <class T, class Less>
void parallelStableSort(std::vector<T> &values, Less less) {
  const int n = static_cast<int>(values.size());
  const int threads = teamThreads();
  if (threads <= 1 || n < threads * 2) {
    std::stable_sort(values.begin(), values.end(), less);
    return;
  }
  const int chunk = (n + threads - 1) / threads;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int tid = 0; tid < threads; tid++) {
    const int lo = tid * chunk;
    const int hi = std::min(n, lo + chunk);
    if (hi - lo > 1) {
      std::stable_sort(values.begin() + lo, values.begin() + hi, less);
    }
  }
  std::vector<T> buffer(values.size());
  T *src = values.data();
  T *dst = buffer.data();
  for (int width = chunk; width < n; width *= 2) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int start = 0; start < n; start += width * 2) {
      const int mid = std::min(n, start + width);
      const int stop = std::min(n, start + width * 2);
      int i = start;
      int j = mid;
      int k = start;
      while (i < mid && j < stop) {
        if (less(src[j], src[i])) {
          dst[k++] = src[j++];
        } else {
          dst[k++] = src[i++];
        }
      }
      while (i < mid) {
        dst[k++] = src[i++];
      }
      while (j < stop) {
        dst[k++] = src[j++];
      }
    }
    std::swap(src, dst);
  }
  if (src != values.data()) {
    std::copy(src, src + n, values.begin());
  }
}

} // namespace seams

#endif
