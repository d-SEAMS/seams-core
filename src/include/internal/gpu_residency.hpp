#ifndef SEAMS_GPU_RESIDENCY_H_
#define SEAMS_GPU_RESIDENCY_H_

#include <cstddef>
#include <string>

/** @file gpu_residency.hpp
 *  @brief Decide whether N frames of the analysis working set fit on the GPU.
 *
 *  The offload unit is a batch of frames of the TUM ice-score path.
 *  A gpulite-style runtime probe reports free device memory. The
 *  footprint is positions, the mutual-graph CSR, the six-ring CSR,
 *  affiliation flags and per-atom HC/DDC labels. CHILL and \(q_{lm}\)
 *  are not in this working set.
 */

namespace gpu {

struct DeviceInfo {
  bool available = false;
  std::string name;
  int computeMajor = 0;
  int computeMinor = 0;
  std::size_t totalBytes = 0;
  std::size_t freeBytes = 0;
  std::string reason;
};

struct Footprint {
  int nAtoms = 0;
  int nFrames = 0;
  int kMax = 16;
  int maxSixRingsPerAtom = 16;
  std::size_t xyzBytes = 0;
  std::size_t nlistBytes = 0;
  std::size_t ringsBytes = 0;
  std::size_t labelBytes = 0;
  std::size_t totalBytes = 0;
};

struct Plan {
  bool resident = false;
  int frames = 0;
  Footprint foot;
  DeviceInfo device;
  std::string reason;
};

DeviceInfo probeDevice();

Footprint estimateFootprint(int nAtoms, int nFrames, int kMax = 16,
                            int maxSixRingsPerAtom = 16);

int maxResidentFrames(const DeviceInfo &dev, int nAtoms, int kMax = 16,
                      int maxSixRingsPerAtom = 16, double safety = 0.80);

Plan planBatch(int nAtoms, int requestedFrames, int kMax = 16,
               int maxSixRingsPerAtom = 16, double safety = 0.80);

} // namespace gpu

#endif
