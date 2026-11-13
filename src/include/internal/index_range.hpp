#ifndef SEAMS_INDEX_RANGE_H_
#define SEAMS_INDEX_RANGE_H_

#include <cstdint>
#include <limits>
#include <vector>

// Cell counts and neighbour-list offsets are signed 32-bit indexes.
// A product or a running sum that does not fit is undefined as a signed
// overflow, and the wrapped value then sizes the wrong array.

namespace seams {

inline bool gridFits(int nx, int ny, int nz, int &cells) {
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    return false;
  }
  const std::int64_t prod = std::int64_t{nx} * ny * nz;
  if (prod > std::numeric_limits<int>::max()) {
    return false;
  }
  cells = static_cast<int>(prod);
  return true;
}

// Drop the longest axis until the cell count fits in a signed int and
// does not exceed maxCells. Axes are not taken below floorDim. A coarser
// cell is at least as wide as the cutoff the caller used to size the
// original grid, so a one-cell halo still covers that cutoff.
// 1<<22 is above the orthorhombic bond grid of an 884736-molecule ice
// cell (about 88^3) and below a product that wraps a signed int.
inline void coarsenUntilFits(int &nx, int &ny, int &nz, int floorDim,
                             int maxCells = 1 << 22) {
  if (floorDim < 1) {
    floorDim = 1;
  }
  if (maxCells < 1) {
    maxCells = 1;
  }
  int cells = 0;
  while (!gridFits(nx, ny, nz, cells) || cells > maxCells) {
    int *axis = &nx;
    if (ny > *axis) {
      axis = &ny;
    }
    if (nz > *axis) {
      axis = &nz;
    }
    if (*axis <= floorDim) {
      return;
    }
    --(*axis);
  }
}

// Exclusive of nothing: offsets[0] = 0 and offsets[i+1] is the sum of
// the first i+1 counts. False when a count is negative or the sum does
// not fit in a signed int, which is the width of the offset array.
inline bool prefixCounts(const int *counts, int n, std::vector<int> &offsets) {
  if (n < 0 || (n > 0 && counts == nullptr)) {
    return false;
  }
  offsets.assign(static_cast<std::size_t>(n) + 1, 0);
  std::int64_t sum = 0;
  for (int i = 0; i < n; ++i) {
    if (counts[i] < 0) {
      return false;
    }
    sum += counts[i];
    if (sum > std::numeric_limits<int>::max()) {
      return false;
    }
    offsets[static_cast<std::size_t>(i) + 1] = static_cast<int>(sum);
  }
  return true;
}

} // namespace seams

#endif
