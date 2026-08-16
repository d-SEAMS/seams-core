#include <gpu_batch.hpp>

#include <chrono>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <vector>

#ifdef SEAMS_HAS_GPULITE
#include <gpulite/gpulite.hpp>
using gpulite::CUDART;
using gpulite::CUDADriver;
using gpulite::KernelFactory;
using gpulite::NVRTC;
using gpulite::dim3;
#endif

namespace gpu {
namespace {

#ifdef SEAMS_HAS_GPULITE

const char *kKernels = R"CUDA(
extern "C" __global__ void nlist_fill(const double* xyz, const double* box,
    int nAtoms, int nFrames, double rc2, int kMax, int* deg, int* cols) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int nTot = nAtoms * nFrames;
  if (tid >= nTot) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const double bx = box[f * 3 + 0];
  const double by = box[f * 3 + 1];
  const double bz = box[f * 3 + 2];
  const int base = f * nAtoms * 3;
  const double ix = xyz[base + i * 3 + 0];
  const double iy = xyz[base + i * 3 + 1];
  const double iz = xyz[base + i * 3 + 2];
  int found = 0;
  const int row = (f * nAtoms + i) * kMax;
  for (int j = 0; j < nAtoms; ++j) {
    if (j == i) continue;
    double dx = xyz[base + j * 3 + 0] - ix;
    double dy = xyz[base + j * 3 + 1] - iy;
    double dz = xyz[base + j * 3 + 2] - iz;
    dx -= bx * nearbyint(dx / bx);
    dy -= by * nearbyint(dy / by);
    dz -= bz * nearbyint(dz / bz);
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 <= rc2 && r2 > 1.0e-12) {
      if (found < kMax) {
        cols[row + found] = j;
      }
      ++found;
    }
  }
  deg[f * nAtoms + i] = found > kMax ? kMax : found;
}

extern "C" __global__ void chill_and_rings(const double* xyz, const double* box,
    const int* deg, const int* cols, int nAtoms, int nFrames, int kMax,
    int* chill, int* sixCount) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int nTot = nAtoms * nFrames;
  if (tid >= nTot) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const int degI = deg[f * nAtoms + i];
  const int row = (f * nAtoms + i) * kMax;
  int take = degI < 4 ? degI : 4;
  const double bx = box[f * 3 + 0];
  const double by = box[f * 3 + 1];
  const double bz = box[f * 3 + 2];
  const int base = f * nAtoms * 3;
  double ix = xyz[base + i * 3 + 0];
  double iy = xyz[base + i * 3 + 1];
  double iz = xyz[base + i * 3 + 2];
  double qre[4] = {0, 0, 0, 0};
  double qim[4] = {0, 0, 0, 0};
  int nn[8];
  for (int a = 0; a < take; ++a) {
    nn[a] = cols[row + a];
    double dx = xyz[base + nn[a] * 3 + 0] - ix;
    double dy = xyz[base + nn[a] * 3 + 1] - iy;
    double dz = xyz[base + nn[a] * 3 + 2] - iz;
    dx -= bx * nearbyint(dx / bx);
    dy -= by * nearbyint(dy / by);
    dz -= bz * nearbyint(dz / bz);
    const double r = sqrt(dx * dx + dy * dy + dz * dz);
    if (r < 1.0e-12) continue;
    const double z = dz / r;
    const double sphi = sqrt(fmax(0.0, 1.0 - z * z));
    const double cphi = (sphi > 1.0e-12) ? dx / (r * sphi) : 1.0;
    const double spsi = (sphi > 1.0e-12) ? dy / (r * sphi) : 0.0;
    // l=3, m=0..3 real/imag pieces (unnormalized relative cij)
    qre[0] += 0.5 * (5.0 * z * z * z - 3.0 * z);
    qre[1] += -0.5 * 1.5 * sphi * (5.0 * z * z - 1.0) * cphi;
    qim[1] += -0.5 * 1.5 * sphi * (5.0 * z * z - 1.0) * spsi;
    qre[2] += 0.5 * 15.0 * sphi * sphi * z * (cphi * cphi - spsi * spsi);
    qim[2] += 0.5 * 15.0 * sphi * sphi * z * (2.0 * cphi * spsi);
    qre[3] += -0.5 * 15.0 * sphi * sphi * sphi *
              (cphi * cphi * cphi - 3.0 * cphi * spsi * spsi);
    qim[3] += -0.5 * 15.0 * sphi * sphi * sphi *
              (3.0 * cphi * cphi * spsi - spsi * spsi * spsi);
  }
  if (take > 0) {
    const double inv = 1.0 / (double)take;
    for (int m = 0; m < 4; ++m) {
      qre[m] *= inv;
      qim[m] *= inv;
    }
  }
  int S = 0;
  int E = 0;
  for (int a = 0; a < take; ++a) {
    const int j = nn[a];
    double jx = xyz[base + j * 3 + 0];
    double jy = xyz[base + j * 3 + 1];
    double jz = xyz[base + j * 3 + 2];
    double qjre[4] = {0, 0, 0, 0};
    double qjim[4] = {0, 0, 0, 0};
    int jdeg = deg[f * nAtoms + j];
    int jtake = jdeg < 4 ? jdeg : 4;
    const int jrow = (f * nAtoms + j) * kMax;
    for (int b = 0; b < jtake; ++b) {
      const int k = cols[jrow + b];
      double dx = xyz[base + k * 3 + 0] - jx;
      double dy = xyz[base + k * 3 + 1] - jy;
      double dz = xyz[base + k * 3 + 2] - jz;
      dx -= bx * nearbyint(dx / bx);
      dy -= by * nearbyint(dy / by);
      dz -= bz * nearbyint(dz / bz);
      const double r = sqrt(dx * dx + dy * dy + dz * dz);
      if (r < 1.0e-12) continue;
      const double z = dz / r;
      const double sphi = sqrt(fmax(0.0, 1.0 - z * z));
      const double cphi = (sphi > 1.0e-12) ? dx / (r * sphi) : 1.0;
      const double spsi = (sphi > 1.0e-12) ? dy / (r * sphi) : 0.0;
      qjre[0] += 0.5 * (5.0 * z * z * z - 3.0 * z);
      qjre[1] += -0.5 * 1.5 * sphi * (5.0 * z * z - 1.0) * cphi;
      qjim[1] += -0.5 * 1.5 * sphi * (5.0 * z * z - 1.0) * spsi;
      qjre[2] += 0.5 * 15.0 * sphi * sphi * z * (cphi * cphi - spsi * spsi);
      qjim[2] += 0.5 * 15.0 * sphi * sphi * z * (2.0 * cphi * spsi);
      qjre[3] += -0.5 * 15.0 * sphi * sphi * sphi *
                 (cphi * cphi * cphi - 3.0 * cphi * spsi * spsi);
      qjim[3] += -0.5 * 15.0 * sphi * sphi * sphi *
                 (3.0 * cphi * cphi * spsi - spsi * spsi * spsi);
    }
    if (jtake > 0) {
      const double inv = 1.0 / (double)jtake;
      for (int m = 0; m < 4; ++m) {
        qjre[m] *= inv;
        qjim[m] *= inv;
      }
    }
    double num = qre[0] * qjre[0];
    double ni = qre[0] * qre[0];
    double nj = qjre[0] * qjre[0];
    for (int m = 1; m < 4; ++m) {
      num += qre[m] * qjre[m] + qim[m] * qjim[m];
      ni += qre[m] * qre[m] + qim[m] * qim[m];
      nj += qjre[m] * qjre[m] + qjim[m] * qjim[m];
    }
    const double den = sqrt(ni * nj);
    const double cij = den > 1.0e-12 ? num / den : 0.0;
    if (cij <= -0.8) ++S;
    else if (cij >= -0.35 && cij <= 0.25) ++E;
  }
  int lab = 0;
  if (S == 4 && E == 0) lab = 1;
  else if (S == 3 && E == 1) lab = 2;
  else if (take == 4) lab = 3;
  chill[f * nAtoms + i] = lab;

  int cycles = 0;
  if (take >= 2) {
    for (int a = 0; a < take; ++a) {
      for (int b = a + 1; b < take; ++b) {
        const int u = nn[a];
        const int v = nn[b];
        const int udeg = deg[f * nAtoms + u];
        const int vdeg = deg[f * nAtoms + v];
        const int urow = (f * nAtoms + u) * kMax;
        const int vrow = (f * nAtoms + v) * kMax;
        for (int p = 0; p < udeg && p < 4; ++p) {
          const int x = cols[urow + p];
          if (x == i || x == v) continue;
          for (int q = 0; q < vdeg && q < 4; ++q) {
            if (cols[vrow + q] == x) {
              ++cycles;
            }
          }
        }
      }
    }
  }
  sixCount[f * nAtoms + i] = cycles;
}
)CUDA";

void checkCuda(cudaError_t st, const char *what) {
  if (st != cudaSuccess) {
    throw std::runtime_error(std::string(what) + ": " +
                             CUDART::instance().cudaGetErrorString(st));
  }
}

struct DevPtr {
  void *p = nullptr;
  ~DevPtr() {
    if (p != nullptr && CUDART::loaded()) {
      CUDART::instance().cudaFree(p);
    }
  }
};

double msSince(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - t0)
      .count();
}

#endif

} // namespace

BatchResult analyzeResident(const double *xyz, const double *box, int nAtoms,
                            int nFrames, double rc) {
  BatchResult out;
  out.plan = planBatch(nAtoms, nFrames);
#ifndef SEAMS_HAS_GPULITE
  out.error = "built without gpulite";
  return out;
#else
  if (!out.plan.resident) {
    out.error = out.plan.reason;
    return out;
  }
  if (!NVRTC::loaded()) {
    out.error = "nvrtc not loaded";
    return out;
  }
  try {
    auto &rt = CUDART::instance();
    const int kMax = out.plan.foot.kMax;
    const std::size_t n = static_cast<std::size_t>(nAtoms);
    const std::size_t nf = static_cast<std::size_t>(out.plan.frames);
    const std::size_t xyzN = nf * n * 3;
    const std::size_t boxN = nf * 3;
    const std::size_t degN = nf * n;
    const std::size_t colN = nf * n * static_cast<std::size_t>(kMax);
    const std::size_t labN = nf * n;

    DevPtr dxyz, dbox, ddeg, dcols, dchill, dsix;
    auto t0 = std::chrono::steady_clock::now();
    checkCuda(rt.cudaMalloc(&dxyz.p, xyzN * sizeof(double)), "cudaMalloc xyz");
    checkCuda(rt.cudaMalloc(&dbox.p, boxN * sizeof(double)), "cudaMalloc box");
    checkCuda(rt.cudaMalloc(&ddeg.p, degN * sizeof(int)), "cudaMalloc deg");
    checkCuda(rt.cudaMalloc(&dcols.p, colN * sizeof(int)), "cudaMalloc cols");
    checkCuda(rt.cudaMalloc(&dchill.p, labN * sizeof(int)), "cudaMalloc chill");
    checkCuda(rt.cudaMalloc(&dsix.p, labN * sizeof(int)), "cudaMalloc six");
    checkCuda(rt.cudaMemcpy(dxyz.p, xyz, xyzN * sizeof(double),
                            cudaMemcpyHostToDevice),
              "HtoD xyz");
    checkCuda(rt.cudaMemcpy(dbox.p, box, boxN * sizeof(double),
                            cudaMemcpyHostToDevice),
              "HtoD box");
    out.uploadMs = msSince(t0);

    auto &factory = KernelFactory::instance(0);
    auto *knList = factory.create("nlist_fill", kKernels, "batch.cu",
                                  {"-std=c++17"});
    auto *kChill = factory.create("chill_and_rings", kKernels, "batch.cu",
                                  {"-std=c++17"});

    const int nTot = nAtoms * out.plan.frames;
    const int block = 128;
    const int grid = (nTot + block - 1) / block;
    const double rc2 = rc * rc;
    t0 = std::chrono::steady_clock::now();
    {
      int nA = nAtoms;
      int nF = out.plan.frames;
      int kM = kMax;
      std::vector<void *> args = {&dxyz.p, &dbox.p, &nA, &nF, (void *)&rc2,
                                  &kM, &ddeg.p, &dcols.p};
      knList->launch(dim3(grid), dim3(block), 0, nullptr, args, true);
    }
    {
      int nA = nAtoms;
      int nF = out.plan.frames;
      int kM = kMax;
      std::vector<void *> args = {&dxyz.p, &dbox.p, &ddeg.p, &dcols.p,
                                  &nA, &nF, &kM, &dchill.p, &dsix.p};
      kChill->launch(dim3(grid), dim3(block), 0, nullptr, args, true);
    }
    out.computeMs = msSince(t0);

    out.chill.assign(labN, 0);
    out.sixCount.assign(labN, 0);
    t0 = std::chrono::steady_clock::now();
    checkCuda(rt.cudaMemcpy(out.chill.data(), dchill.p, labN * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH chill");
    checkCuda(rt.cudaMemcpy(out.sixCount.data(), dsix.p, labN * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH six");
    out.downloadMs = msSince(t0);
  } catch (const std::exception &ex) {
    out.error = ex.what();
  }
  return out;
#endif
}

} // namespace gpu
