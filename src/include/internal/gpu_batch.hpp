#ifndef SEAMS_GPU_BATCH_H_
#define SEAMS_GPU_BATCH_H_

#include <gpu_residency.hpp>

#include <string>
#include <vector>

/** @file gpu_batch.hpp
 *  @brief Run the analysis of a resident frame batch on the device.
 *
 *  Coordinates for every accepted frame are uploaded once. The device
 *  path is the TUM ice score: linkcell \(k\)-nearest, the mutual
 *  four-nearest graph, primitive six-rings, and the HC/DDC affiliation
 *  predicates. CHILL and \(q_{lm}\) stay on the host. Only cage labels
 *  come back.
 */

namespace gpu {

struct BatchResult {
  Plan plan;
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
