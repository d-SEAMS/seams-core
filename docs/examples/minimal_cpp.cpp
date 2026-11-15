/* Read the GenIce sI cell through the C++ language header. */
#include "seams_lang.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv) {
  const char *path = argc > 1 ? argv[1] : "input/traj/genice_sI.lammpstrj";
  seams_frame *frame = seams_frame_read(path, 1, 0);
  if (frame == nullptr || seams_frame_nop(frame) != 46) {
    std::fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  if (seams_frame_neighbours(frame, 3.5, 0, "cutoff") != 0) {
    return 1;
  }
  if (seams_frame_neighbour_pairs(frame) < 150) {
    return 1;
  }
  double q6[46];
  if (seams_frame_steinhardt(frame, 3.5, 0, 6, q6, nullptr) != 0) {
    return 1;
  }
  double mean = 0;
  for (double value : q6) {
    mean += value;
  }
  mean /= 46.0;
  int labels[46];
  int counts[9];
  if (seams_frame_chill_plus(frame, 3.5, 0, "cutoff", labels, counts) != 0) {
    return 1;
  }
  int rings[7] = {};
  if (seams_frame_rings(frame, 3.5, 0, "cutoff", 6, rings, 7) != 0) {
    return 1;
  }
  int n_cages = 0;
  int n_sig_atoms = 0;
  int n_vertices = 0;
  if (seams_frame_signature(frame, 3.5, 0, "cutoff", "512", &n_cages, &n_sig_atoms) != 0) {
    return 1;
  }
  int n_edge = 0;
  int n_edge_atoms = 0;
  if (seams_frame_signature_edge(frame, 3.5, 0, "cutoff", "512", "none", &n_edge,
                                 &n_edge_atoms) != 0 ||
      n_edge != n_cages) {
    return 1;
  }
  int n_tight = 0;
  int n_tight_atoms = 0;
  if (seams_frame_signature_edge(frame, 3.5, 0, "cutoff", "512", "1,2", &n_tight,
                                 &n_tight_atoms) != 0 ||
      n_tight != 0) {
    return 1;
  }
  if (n_cages > 0 &&
      seams_frame_signature_vertices(frame, 3.5, 0, "cutoff", "512", 0, nullptr, 0,
                                     &n_vertices) != 0) {
    return 1;
  }
  if (n_cages > 0 &&
      seams_frame_signature_vertices_edge(frame, 3.5, 0, "cutoff", "512", "none", 0, nullptr, 0,
                                          &n_vertices) != 0) {
    return 1;
  }
  if (counts[4] != 46 || rings[5] != 48 || rings[6] != 54 || mean <= 0.5) {
    return 1;
  }
  std::printf("%s atoms 46 clathrate %d q6 %.3f rings %d\n", seams_lang_version(), counts[4],
              mean, rings[5] + rings[6]);
  seams_frame_free(frame);
  return 0;
}
