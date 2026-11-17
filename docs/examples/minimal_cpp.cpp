// Copyright (c) d-SEAMS developers
/* Read the GenIce sI cell through the C++ language header. */
#include "seams_lang.h"

#include <cstdio>
#include <cstdlib>
#include <vector>

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
  int n_five = 0;
  int per_atom[46];
  if (seams_frame_chordless(frame, 3.5, 0, "cutoff", 5, nullptr, 0, &n_five, per_atom) != 0) {
    return 1;
  }
  int five[48 * 5];
  if (n_five != 48 ||
      seams_frame_chordless(frame, 3.5, 0, "cutoff", 5, five, 48 * 5, &n_five, nullptr) != 0) {
    return 1;
  }
  int through = 0;
  for (int n : per_atom) {
    through += n;
  }
  if (through != 5 * 48 || five[0] != 0) {
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
  // CHILL+ and rings from one neighbour list, with the phase times.
  int fused[9] = {};
  int fused_rings[7] = {};
  seams_phase_times times{};
  if (seams_frame_chill_rings(frame, 3.5, 0, "cutoff", 6, 0.0, labels, fused, fused_rings, 7,
                              &times) != 0 ||
      fused[4] != counts[4]) {
    return 1;
  }
  int shell[46];
  if (seams_frame_symmetry(frame, 0.12, shell) != 0) {
    return 1;
  }
  if (seams_frame_rdf_nbins(6.0, 0) < 1) {
    return 1;
  }
  const int n_r2 = seams_frame_rdf2d_nbin(3.5, 0.5);
  std::vector<double> r2(static_cast<std::size_t>(n_r2));
  std::vector<double> g2(static_cast<std::size_t>(n_r2));
  if (seams_frame_rdf2d(frame, 3.5, 0.5, nullptr, r2.data(), g2.data()) != 0) {
    return 1;
  }
  int sheet[8] = {};
  if (seams_frame_monolayer(frame, 3.5, 0, 6, 1.0, nullptr, sheet, 8) != 0) {
    return 1;
  }
  // The same positions as a new frame with one molecule per atom.
  double box9[9];
  seams_frame_box(frame, box9);
  std::vector<int> types(46, 1);
  std::vector<int> mols(46);
  for (int i = 0; i < 46; ++i) {
    mols[static_cast<std::size_t>(i)] = i + 1;
  }
  seams_frame *copy = seams_frame_from_species(seams_frame_positions(frame), 46, box9,
                                               types.data(), nullptr, mols.data());
  if (copy == nullptr) {
    return 1;
  }
  double tau = 0.0;
  const int jump = seams_frames_jump_tau90(frame, copy, 1.0, 1, 2, &tau);
  seams_frame_free(copy);
  if (jump != 0) {
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
