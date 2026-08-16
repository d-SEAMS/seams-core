/*
** Probe whether N frames of the TUM ice-score working set fit on the
** GPU, then time a cold launch and the best of N warm repeats.
**   bench_gpu_batch TRAJ [nFrames] [atomType] [repeats]
*/

#include <gpu_batch.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  const std::string traj = argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int want = argc > 2 ? std::atoi(argv[2]) : 11;
  const int typeI = argc > 3 ? std::atoi(argv[3]) : 1;
  const int repeats = argc > 4 ? std::max(1, std::atoi(argv[4])) : 5;

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, 1, cloud, typeI);
  if (cloud.nop == 0) {
    std::fprintf(stderr, "empty %s\n", traj.c_str());
    return 1;
  }
  const int nAtoms = cloud.nop;
  std::vector<double> xyz(static_cast<std::size_t>(want) *
                          static_cast<std::size_t>(nAtoms) * 3);
  std::vector<double> box(static_cast<std::size_t>(want) * 3);
  int got = 0;
  for (int f = 1; f <= want; ++f) {
    cloud = sinp::readLammpsTrjO(traj, f, cloud, typeI);
    if (cloud.nop != nAtoms) {
      break;
    }
    const auto base = static_cast<std::size_t>(got) *
                      static_cast<std::size_t>(nAtoms) * 3;
    for (int i = 0; i < nAtoms; ++i) {
      xyz[base + static_cast<std::size_t>(i) * 3 + 0] = cloud.pts[i].x;
      xyz[base + static_cast<std::size_t>(i) * 3 + 1] = cloud.pts[i].y;
      xyz[base + static_cast<std::size_t>(i) * 3 + 2] = cloud.pts[i].z;
    }
    box[static_cast<std::size_t>(got) * 3 + 0] = cloud.box[0];
    box[static_cast<std::size_t>(got) * 3 + 1] = cloud.box[1];
    box[static_cast<std::size_t>(got) * 3 + 2] = cloud.box[2];
    ++got;
  }

  const auto plan = gpu::planBatch(nAtoms, got);
  std::printf("device %s\n",
              plan.device.available ? plan.device.name.c_str() : "none");
  std::printf("reason %s\n", plan.reason.c_str());
  std::printf("sm %d.%d\n", plan.device.computeMajor,
              plan.device.computeMinor);
  std::printf("total_bytes %zu\n", plan.device.totalBytes);
  std::printf("free_bytes %zu\n", plan.device.freeBytes);
  std::printf("nAtoms %d\n", nAtoms);
  std::printf("requested_frames %d\n", got);
  std::printf("bytes_per_frame %zu\n",
              gpu::estimateFootprint(nAtoms, 1).totalBytes);
  std::printf("max_resident_frames %d\n",
              gpu::maxResidentFrames(plan.device, nAtoms));
  std::printf("plan_frames %d\n", plan.frames);
  std::printf("resident %s\n", plan.resident ? "yes" : "no");
  if (!plan.resident) {
    return 0;
  }
  const auto printRun = [](const char *tag, const gpu::BatchResult &r) {
    std::printf("%s_upload_ms %.3f\n%s_compute_ms %.3f\n%s_download_ms %.3f\n",
                tag, r.uploadMs, tag, r.computeMs, tag, r.downloadMs);
  };

  const auto cold = gpu::analyzeResident(xyz.data(), box.data(), nAtoms, got);
  if (!cold.error.empty()) {
    std::printf("error %s\n", cold.error.c_str());
    return 1;
  }
  printRun("cold", cold);

  gpu::BatchResult best = cold;
  best.computeMs = std::numeric_limits<double>::max();
  for (int i = 0; i < repeats; ++i) {
    const auto r = gpu::analyzeResident(xyz.data(), box.data(), nAtoms, got);
    if (!r.error.empty()) {
      std::printf("error %s\n", r.error.c_str());
      return 1;
    }
    if (r.computeMs < best.computeMs) {
      best = r;
    }
  }
  printRun("warm", best);
  int nHc = 0;
  int nDdc = 0;
  int nSix = 0;
  for (int v : best.atomHc) {
    nHc += v;
  }
  for (int v : best.atomDdc) {
    nDdc += v;
  }
  for (int v : best.nRings) {
    nSix += v;
  }
  std::printf("ice_hc_atoms %d\nice_ddc_atoms %d\n", nHc, nDdc);
  std::printf("six_rings %d\nrings_dropped %d\n", nSix, best.ringsDropped);
  return 0;
}
