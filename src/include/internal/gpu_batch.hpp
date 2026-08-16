#ifndef SEAMS_GPU_BATCH_H_
#define SEAMS_GPU_BATCH_H_

#include <gpu_residency.hpp>

#include <string>
#include <vector>

/** @file gpu_batch.hpp
 *  @brief Run the analysis of a resident frame batch on the device.
 *
 *  Coordinates for every accepted frame are uploaded once. The device
 *  builds a linked cell list the way vesin CUDA and LAMMPS NPairKokkos
 *  do (host bin grid, sort-by-cell, exclusive scan, stencil walk on
 *  coalesced positions, persistent buffers), then the mutual
 *  four-nearest graph, one \(q_{lm}\) write, CHILL, primitive
 *  six-rings and the HC/DDC affiliation predicates. Only labels come
 *  back.
 */

namespace gpu {

struct BatchResult {
  Plan plan;
  std::vector<int> chill;   // nFrames * nAtoms
  std::vector<int> sixCount;
  std::vector<int> atomHc;  // 1 if the atom sits on an HC-affiliated ring
  std::vector<int> atomDdc; // 1 if the atom sits on a DDC-affiliated ring
  std::vector<int> nRings;  // per frame
  int ringsDropped = 0;
  double uploadMs = 0.0;
  double computeMs = 0.0;
  double downloadMs = 0.0;
  std::string error;
};

/** xyz is frame-major, atom, xyz. box is frame-major, three lengths.
 *  rc is the cell-list cutoff; 5.5 matches the host k-NN candidate shell. */
BatchResult analyzeResident(const double *xyz, const double *box, int nAtoms,
                            int nFrames, double rc = 5.5);

} // namespace gpu

#endif
