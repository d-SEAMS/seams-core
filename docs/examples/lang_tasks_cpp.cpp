// Neighbours, rings, and CHILL+ through the C++ language front.
// Frame owns the seams_frame pointer and frees it on every exit.
#include "seams_lang.h"

#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Frame {
  seams_frame *ptr;

  explicit Frame(const std::string &path) : ptr(nullptr) {
    ptr = seams_frame_read(path.c_str(), 1, -1);
    if (ptr == nullptr) {
      throw std::runtime_error(seams_last_error());
    }
  }

  Frame(const Frame &) = delete;
  Frame &operator=(const Frame &) = delete;

  ~Frame() { seams_frame_free(ptr); }

  int pairs() const {
    if (seams_frame_neighbours(ptr, 3.5, -1, "cutoff") != 0) {
      throw std::runtime_error(seams_last_error());
    }
    return seams_frame_neighbour_pairs(ptr);
  }

  std::vector<int> rings() const {
    std::vector<int> counts(7, 0);
    if (seams_frame_rings(ptr, 3.5, -1, "cutoff", 6, counts.data(),
                          static_cast<int>(counts.size())) != 0) {
      throw std::runtime_error(seams_last_error());
    }
    return counts;
  }

  std::vector<int> chill() const {
    const int n = seams_frame_nop(ptr);
    std::vector<int> labels(static_cast<std::size_t>(n), 0);
    std::vector<int> counts(9, 0);
    if (seams_frame_chill_plus(ptr, 3.5, -1, "cutoff", labels.data(),
                               counts.data()) != 0) {
      throw std::runtime_error(seams_last_error());
    }
    return counts;
  }
};

}  // namespace

int main() {
  const char *root = std::getenv("DSEAMS_ROOT");
  if (root == nullptr) {
    std::fprintf(stderr, "DSEAMS_ROOT is unset\n");
    return 1;
  }
  try {
    Frame frame(std::string(root) + "/input/traj/ice1h.lammpstrj");
    const int pairs = frame.pairs();
    const std::vector<int> rings = frame.rings();
    const std::vector<int> chill = frame.chill();
    std::printf("pairs %d rings", pairs);
    for (int count : rings) {
      std::printf(" %d", count);
    }
    std::printf(" chill");
    for (int count : chill) {
      std::printf(" %d", count);
    }
    std::printf("\n");
  } catch (const std::exception &err) {
    std::fprintf(stderr, "%s\n", err.what());
    return 1;
  }
  return 0;
}
