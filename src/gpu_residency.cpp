#include <gpu_residency.hpp>

#include <dlfcn.h>

#include <algorithm>
#include <cstdint>

#ifdef SEAMS_HAS_GPULITE
#include <gpulite/gpulite.hpp>
using gpulite::CUDART;
using gpulite::CUDADriver;
#endif

namespace gpu {
namespace {

constexpr std::size_t kAlign = 256;

std::size_t alignUp(std::size_t n) {
  return (n + kAlign - 1) & ~(kAlign - 1);
}

#ifdef SEAMS_HAS_GPULITE
using cudaMemGetInfo_t = cudaError_t (*)(std::size_t *, std::size_t *);

cudaMemGetInfo_t loadMemGetInfo() {
  static cudaMemGetInfo_t fn = nullptr;
  static bool tried = false;
  if (tried) {
    return fn;
  }
  tried = true;
  static const char *cands[] = {"libcudart.so", "libcudart.so.12",
                                "libcudart.so.13", "libcudart.so.11", nullptr};
  for (int i = 0; cands[i] != nullptr; i++) {
    void *h = dlopen(cands[i], RTLD_NOW | RTLD_NOLOAD);
    if (h == nullptr) {
      h = dlopen(cands[i], RTLD_NOW);
    }
    if (h == nullptr) {
      continue;
    }
    fn = reinterpret_cast<cudaMemGetInfo_t>(dlsym(h, "cudaMemGetInfo"));
    if (fn != nullptr) {
      return fn;
    }
  }
  return nullptr;
}
#endif

} // namespace

Footprint estimateFootprint(int nAtoms, int nFrames, int kMax,
                            int maxSixRingsPerAtom) {
  Footprint f;
  f.nAtoms = nAtoms;
  f.nFrames = nFrames;
  f.kMax = kMax;
  f.maxSixRingsPerAtom = maxSixRingsPerAtom;
  const std::size_t n = static_cast<std::size_t>(nAtoms);
  const std::size_t nf = static_cast<std::size_t>(nFrames);
  const std::size_t k = static_cast<std::size_t>(kMax);
  const std::size_t r = static_cast<std::size_t>(maxSixRingsPerAtom);
  f.xyzBytes = alignUp(nf * n * 3 * sizeof(double));
  const std::size_t boxBytes = alignUp(nf * 3 * sizeof(double));
  f.nlistBytes = alignUp(nf * n * sizeof(int)) +
                 alignUp(nf * n * k * sizeof(int)) +
                 alignUp(nf * n * 4 * sizeof(int));
  const std::size_t maxRings = n * r;
  f.ringsBytes = alignUp(nf * maxRings * 6 * sizeof(int)) +
                 alignUp(nf * sizeof(int)) +
                 alignUp(nf * n * r * sizeof(int)) +
                 alignUp(nf * n * sizeof(int));
  f.labelBytes = alignUp(nf * n * 2 * sizeof(int));
  const std::size_t flags = alignUp(nf * maxRings * 2 * sizeof(int));
  const std::size_t cells = alignUp(nf * n * sizeof(int));
  f.totalBytes = f.xyzBytes + boxBytes + f.nlistBytes + f.ringsBytes +
                 f.labelBytes + flags + cells;
  return f;
}

DeviceInfo probeDevice() {
  DeviceInfo info;
#ifndef SEAMS_HAS_GPULITE
  info.reason = "built without gpulite";
  return info;
#else
  if (!CUDADriver::loaded()) {
    info.reason = "CUDA driver not loaded";
    return info;
  }
  if (!CUDART::loaded()) {
    info.reason = "cudart not loaded";
    return info;
  }
  int nDev = 0;
  if (CUDART::instance().cudaGetDeviceCount(&nDev) != cudaSuccess || nDev < 1) {
    info.reason = "no CUDA device";
    return info;
  }
  CUdevice dev = 0;
  char name[256] = {};
  CUDADriver::instance().cuDeviceGetName(name, 255, dev);
  info.name = name;
  int maj = 0;
  int minv = 0;
  CUDADriver::instance().cuDeviceGetAttribute(
      &maj, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, dev);
  CUDADriver::instance().cuDeviceGetAttribute(
      &minv, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, dev);
  info.computeMajor = maj;
  info.computeMinor = minv;
  std::size_t total = 0;
  CUDADriver::instance().cuDeviceTotalMem(&total, dev);
  info.totalBytes = total;
  info.freeBytes = total;
  auto memInfo = loadMemGetInfo();
  if (memInfo != nullptr) {
    std::size_t freeB = 0;
    std::size_t totB = 0;
    if (memInfo(&freeB, &totB) == cudaSuccess) {
      info.freeBytes = freeB;
      if (totB > 0) {
        info.totalBytes = totB;
      }
    }
  }
  info.available = info.totalBytes > 0;
  if (!info.available) {
    info.reason = "device reported zero memory";
  }
  return info;
#endif
}

int maxResidentFrames(const DeviceInfo &dev, int nAtoms, int kMax,
                      int maxSixRingsPerAtom, double safety) {
  if (!dev.available || nAtoms <= 0) {
    return 0;
  }
  const auto one = estimateFootprint(nAtoms, 1, kMax, maxSixRingsPerAtom);
  if (one.totalBytes == 0) {
    return 0;
  }
  const auto budget =
      static_cast<std::size_t>(static_cast<double>(dev.freeBytes) * safety);
  return static_cast<int>(budget / one.totalBytes);
}

Plan planBatch(int nAtoms, int requestedFrames, int kMax,
               int maxSixRingsPerAtom, double safety) {
  Plan p;
  p.device = probeDevice();
  p.foot = estimateFootprint(nAtoms, requestedFrames, kMax, maxSixRingsPerAtom);
  if (!p.device.available) {
    p.reason = p.device.reason;
    return p;
  }
  const int maxF =
      maxResidentFrames(p.device, nAtoms, kMax, maxSixRingsPerAtom, safety);
  p.frames = std::min(requestedFrames, maxF);
  p.foot = estimateFootprint(nAtoms, p.frames, kMax, maxSixRingsPerAtom);
  if (p.frames <= 0) {
    p.reason = "working set exceeds free device memory";
    return p;
  }
  p.resident = true;
  p.reason = "batch fits under safety margin";
  return p;
}

} // namespace gpu
