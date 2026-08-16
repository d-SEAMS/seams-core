#ifndef SEAMS_GPU_BATCH_H_
#define SEAMS_GPU_BATCH_H_

#include <gpu_residency.hpp>

#include <string>
#include <vector>

/** @file gpu_batch.hpp
 *  @brief Run the analysis of a resident frame batch on the device.
 *
 *  Coordinates for every accepted frame are uploaded once. Neighbour
 *  lists, six-ring incidence and CHILL labels are produced on the
 *  device. Only the labels (and a six-ring count per atom) come back.
 */

namespace gpu {

struct BatchResult {
  Plan plan;
  std::vector<int> chill;     // nFrames * nAtoms
  std::vector<int> sixCount;  // nFrames * nAtoms
  double uploadMs = 0.0;
  double computeMs = 0.0;
  double downloadMs = 0.0;
  std::string error;
};

/** xyz is frame-major, atom, xyz. box is frame-major, three lengths. */
BatchResult analyzeResident(const double *xyz, const double *box, int nAtoms,
                            int nFrames, double rc = 3.5);

} // namespace gpu

#endif
