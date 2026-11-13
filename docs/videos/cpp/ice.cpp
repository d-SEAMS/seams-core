/* Ice Ih census through seams_lang.h. */
#include "seams_lang.h"

#include <cstdio>
#include <vector>

int main(int argc, char **argv) {
  const char *path = argc > 1 ? argv[1] : "input/traj/ice1h.lammpstrj";
  seams_frame *frame = seams_frame_read(path, 1, 1);
  if (frame == nullptr) {
    std::fprintf(stderr, "read failed\n");
    return 1;
  }
  const int n = seams_frame_nop(frame);
  if (n != 128) {
    std::fprintf(stderr, "nop %d\n", n);
    return 1;
  }
  if (seams_frame_neighbours(frame, 3.5, 1, "cutoff") != 0) {
    return 1;
  }
  const int pairs = seams_frame_neighbour_pairs(frame);
  std::vector<double> q6(static_cast<std::size_t>(n));
  if (seams_frame_steinhardt_graph(frame, 3.5, 1, 6, "cutoff", q6.data(), nullptr) != 0) {
    return 1;
  }
  double mean = 0;
  for (double value : q6) {
    mean += value;
  }
  mean /= n;
  std::vector<int> labels(static_cast<std::size_t>(n));
  int counts[9] = {};
  if (seams_frame_chill_plus(frame, 3.5, 1, "cutoff", labels.data(), counts) != 0) {
    return 1;
  }
  int rings[7] = {};
  if (seams_frame_rings(frame, 3.5, 1, "cutoff", 6, rings, 7) != 0) {
    return 1;
  }
  std::printf("nop %d pairs %d q6 %.6f hexagonal %d rings %d\n", n, pairs, mean, counts[1],
              rings[6]);
  seams_frame_free(frame);
  return 0;
}
