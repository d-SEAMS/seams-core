#include <gpu_batch.hpp>

#include <chrono>
#include <cmath>
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

// Cell-list neighbours, one q_lm pass, then CHILL + 4-NN six-cycles.
// All frames stay in the same device buffers.
const char *kKernels = R"CUDA(
__device__ inline void minImage(double& dx, double& dy, double& dz,
    double bx, double by, double bz) {
  dx -= bx * nearbyint(dx / bx);
  dy -= by * nearbyint(dy / by);
  dz -= bz * nearbyint(dz / bz);
}

__device__ inline void y3(double dx, double dy, double dz,
    double* qre, double* qim) {
  const double r = sqrt(dx * dx + dy * dy + dz * dz);
  if (r < 1.0e-12) return;
  const double z = dz / r;
  const double sphi = sqrt(fmax(0.0, 1.0 - z * z));
  const double cphi = (sphi > 1.0e-12) ? dx / (r * sphi) : 1.0;
  const double spsi = (sphi > 1.0e-12) ? dy / (r * sphi) : 0.0;
  qre[0] += 0.5 * (5.0 * z * z * z - 3.0 * z);
  qre[1] += -0.75 * sphi * (5.0 * z * z - 1.0) * cphi;
  qim[1] += -0.75 * sphi * (5.0 * z * z - 1.0) * spsi;
  qre[2] += 7.5 * sphi * sphi * z * (cphi * cphi - spsi * spsi);
  qim[2] += 7.5 * sphi * sphi * z * (2.0 * cphi * spsi);
  qre[3] += -7.5 * sphi * sphi * sphi *
            (cphi * cphi * cphi - 3.0 * cphi * spsi * spsi);
  qim[3] += -7.5 * sphi * sphi * sphi *
            (3.0 * cphi * cphi * spsi - spsi * spsi * spsi);
}

extern "C" __global__ void bin_atoms(const double* xyz, const double* box,
    int nAtoms, int nFrames, double rc, int* ncell, int* cellOf, int* cellCount) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const double bx = box[f * 3 + 0];
  const double by = box[f * 3 + 1];
  const double bz = box[f * 3 + 2];
  int nx = (int)floor(bx / rc); if (nx < 1) nx = 1;
  int ny = (int)floor(by / rc); if (ny < 1) ny = 1;
  int nz = (int)floor(bz / rc); if (nz < 1) nz = 1;
  while (nx * ny * nz > nAtoms) {
    if (nx >= ny && nx >= nz) --nx;
    else if (ny >= nz) --ny;
    else --nz;
    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    if (nz < 1) nz = 1;
  }
  if (i == 0) {
    ncell[f * 3 + 0] = nx;
    ncell[f * 3 + 1] = ny;
    ncell[f * 3 + 2] = nz;
  }
  const int base = f * nAtoms * 3;
  double fx = xyz[base + i * 3 + 0] / bx;
  double fy = xyz[base + i * 3 + 1] / by;
  double fz = xyz[base + i * 3 + 2] / bz;
  fx -= floor(fx); fy -= floor(fy); fz -= floor(fz);
  int cx = (int)(fx * nx); if (cx < 0) cx = 0; if (cx >= nx) cx = nx - 1;
  int cy = (int)(fy * ny); if (cy < 0) cy = 0; if (cy >= ny) cy = ny - 1;
  int cz = (int)(fz * nz); if (cz < 0) cz = 0; if (cz >= nz) cz = nz - 1;
  const int cid = (cz * ny + cy) * nx + cx;
  cellOf[f * nAtoms + i] = cid;
  atomicAdd(cellCount + f * nAtoms + cid, 1);
}

extern "C" __global__ void prefix_cells(const int* ncell, int* cellCount,
    int* cellOff, int nAtoms, int nFrames) {
  const int f = blockIdx.x * blockDim.x + threadIdx.x;
  if (f >= nFrames) return;
  const int nx = ncell[f * 3 + 0];
  const int ny = ncell[f * 3 + 1];
  const int nz = ncell[f * 3 + 2];
  const int nC = nx * ny * nz;
  int acc = 0;
  for (int c = 0; c < nC; ++c) {
    cellOff[f * (nAtoms + 1) + c] = acc;
    acc += cellCount[f * nAtoms + c];
    cellCount[f * nAtoms + c] = 0;
  }
  cellOff[f * (nAtoms + 1) + nC] = acc;
}

extern "C" __global__ void scatter_atoms(const int* cellOf, int* cellCount,
    const int* cellOff, int* order, int nAtoms, int nFrames) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const int cid = cellOf[f * nAtoms + i];
  const int slot = atomicAdd(cellCount + f * nAtoms + cid, 1);
  order[f * nAtoms + cellOff[f * (nAtoms + 1) + cid] + slot] = i;
}

extern "C" __global__ void nlist_cells(const double* xyz, const double* box,
    const int* ncell, const int* cellOff, const int* order,
    int nAtoms, int nFrames, double rc2, int kMax, int* deg, int* cols) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const double bx = box[f * 3 + 0];
  const double by = box[f * 3 + 1];
  const double bz = box[f * 3 + 2];
  const int nx = ncell[f * 3 + 0];
  const int ny = ncell[f * 3 + 1];
  const int nz = ncell[f * 3 + 2];
  const int base = f * nAtoms * 3;
  const double ix = xyz[base + i * 3 + 0];
  const double iy = xyz[base + i * 3 + 1];
  const double iz = xyz[base + i * 3 + 2];
  double fx = ix / bx; fx -= floor(fx);
  double fy = iy / by; fy -= floor(fy);
  double fz = iz / bz; fz -= floor(fz);
  int cx = (int)(fx * nx); if (cx < 0) cx = 0; if (cx >= nx) cx = nx - 1;
  int cy = (int)(fy * ny); if (cy < 0) cy = 0; if (cy >= ny) cy = ny - 1;
  int cz = (int)(fz * nz); if (cz < 0) cz = 0; if (cz >= nz) cz = nz - 1;
  int found = 0;
  const int row = (f * nAtoms + i) * kMax;
  double bestR2[16];
  int bestJ[16];
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dz = -1; dz <= 1; ++dz) {
        int ncx = cx + dx; int ncy = cy + dy; int ncz = cz + dz;
        ncx = (ncx % nx + nx) % nx;
        ncy = (ncy % ny + ny) % ny;
        ncz = (ncz % nz + nz) % nz;
        const int cid = (ncz * ny + ncy) * nx + ncx;
        const int a0 = cellOff[f * (nAtoms + 1) + cid];
        const int a1 = cellOff[f * (nAtoms + 1) + cid + 1];
        for (int s = a0; s < a1; ++s) {
          const int j = order[f * nAtoms + s];
          if (j == i) continue;
          double rx = xyz[base + j * 3 + 0] - ix;
          double ry = xyz[base + j * 3 + 1] - iy;
          double rz = xyz[base + j * 3 + 2] - iz;
          minImage(rx, ry, rz, bx, by, bz);
          const double r2 = rx * rx + ry * ry + rz * rz;
          if (r2 > rc2 || r2 <= 1.0e-12) continue;
          if (found < kMax) {
            bestR2[found] = r2;
            bestJ[found] = j;
            ++found;
          } else {
            int w = 0;
            for (int t = 1; t < kMax; ++t) {
              if (bestR2[t] > bestR2[w]) w = t;
            }
            if (r2 < bestR2[w]) {
              bestR2[w] = r2;
              bestJ[w] = j;
            }
          }
        }
      }
    }
  }
  for (int a = 0; a < found; ++a) cols[row + a] = bestJ[a];
  deg[f * nAtoms + i] = found;
}

extern "C" __global__ void qlm_l3(const double* xyz, const double* box,
    const int* deg, const int* cols, int nAtoms, int nFrames, int kMax,
    double* qre, double* qim) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const int take0 = deg[f * nAtoms + i];
  const int take = take0 < 4 ? take0 : 4;
  const int row = (f * nAtoms + i) * kMax;
  const int base = f * nAtoms * 3;
  const double bx = box[f * 3 + 0];
  const double by = box[f * 3 + 1];
  const double bz = box[f * 3 + 2];
  const double ix = xyz[base + i * 3 + 0];
  const double iy = xyz[base + i * 3 + 1];
  const double iz = xyz[base + i * 3 + 2];
  double re[4] = {0, 0, 0, 0};
  double im[4] = {0, 0, 0, 0};
  for (int a = 0; a < take; ++a) {
    const int j = cols[row + a];
    double dx = xyz[base + j * 3 + 0] - ix;
    double dy = xyz[base + j * 3 + 1] - iy;
    double dz = xyz[base + j * 3 + 2] - iz;
    minImage(dx, dy, dz, bx, by, bz);
    y3(dx, dy, dz, re, im);
  }
  if (take > 0) {
    const double inv = 1.0 / (double)take;
    for (int m = 0; m < 4; ++m) {
      re[m] *= inv;
      im[m] *= inv;
    }
  }
  const int qrow = (f * nAtoms + i) * 4;
  for (int m = 0; m < 4; ++m) {
    qre[qrow + m] = re[m];
    qim[qrow + m] = im[m];
  }
}

extern "C" __global__ void chill_from_qlm(const int* deg, const int* cols,
    const double* qre, const double* qim, int nAtoms, int nFrames, int kMax,
    int* chill, int* sixCount) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const int take0 = deg[f * nAtoms + i];
  const int take = take0 < 4 ? take0 : 4;
  const int row = (f * nAtoms + i) * kMax;
  const int qi = (f * nAtoms + i) * 4;
  int S = 0;
  int E = 0;
  int nn[8];
  for (int a = 0; a < take; ++a) {
    nn[a] = cols[row + a];
    const int qj = (f * nAtoms + nn[a]) * 4;
    double num = qre[qi] * qre[qj];
    double ni = qre[qi] * qre[qi];
    double nj = qre[qj] * qre[qj];
    for (int m = 1; m < 4; ++m) {
      num += qre[qi + m] * qre[qj + m] + qim[qi + m] * qim[qj + m];
      ni += qre[qi + m] * qre[qi + m] + qim[qi + m] * qim[qi + m];
      nj += qre[qj + m] * qre[qj + m] + qim[qj + m] * qim[qj + m];
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
        const int udeg = deg[f * nAtoms + u] < 4 ? deg[f * nAtoms + u] : 4;
        const int vdeg = deg[f * nAtoms + v] < 4 ? deg[f * nAtoms + v] : 4;
        const int urow = (f * nAtoms + u) * kMax;
        const int vrow = (f * nAtoms + v) * kMax;
        for (int p = 0; p < udeg; ++p) {
          const int x = cols[urow + p];
          if (x == i || x == v) continue;
          for (int q = 0; q < vdeg; ++q) {
            if (cols[vrow + q] == x) ++cycles;
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

void *argp(void *p) { return p; }

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
    const std::size_t nN = nf * n;
    const std::size_t colN = nN * static_cast<std::size_t>(kMax);
    const std::size_t qN = nN * 4;
    const std::size_t offN = nf * (n + 1);

    DevPtr dxyz, dbox, dncell, dcellOf, dcellCount, dcellOff, dorder;
    DevPtr ddeg, dcols, dqre, dqim, dchill, dsix;
    auto t0 = std::chrono::steady_clock::now();
    checkCuda(rt.cudaMalloc(&dxyz.p, xyzN * sizeof(double)), "xyz");
    checkCuda(rt.cudaMalloc(&dbox.p, boxN * sizeof(double)), "box");
    checkCuda(rt.cudaMalloc(&dncell.p, nf * 3 * sizeof(int)), "ncell");
    checkCuda(rt.cudaMalloc(&dcellOf.p, nN * sizeof(int)), "cellOf");
    checkCuda(rt.cudaMalloc(&dcellCount.p, nN * sizeof(int)), "cellCount");
    checkCuda(rt.cudaMalloc(&dcellOff.p, offN * sizeof(int)), "cellOff");
    checkCuda(rt.cudaMalloc(&dorder.p, nN * sizeof(int)), "order");
    checkCuda(rt.cudaMalloc(&ddeg.p, nN * sizeof(int)), "deg");
    checkCuda(rt.cudaMalloc(&dcols.p, colN * sizeof(int)), "cols");
    checkCuda(rt.cudaMalloc(&dqre.p, qN * sizeof(double)), "qre");
    checkCuda(rt.cudaMalloc(&dqim.p, qN * sizeof(double)), "qim");
    checkCuda(rt.cudaMalloc(&dchill.p, nN * sizeof(int)), "chill");
    checkCuda(rt.cudaMalloc(&dsix.p, nN * sizeof(int)), "six");
    checkCuda(rt.cudaMemset(dcellCount.p, 0, nN * sizeof(int)), "zero cells");
    checkCuda(rt.cudaMemcpy(dxyz.p, xyz, xyzN * sizeof(double),
                            cudaMemcpyHostToDevice),
              "HtoD xyz");
    checkCuda(rt.cudaMemcpy(dbox.p, box, boxN * sizeof(double),
                            cudaMemcpyHostToDevice),
              "HtoD box");
    out.uploadMs = msSince(t0);

    auto &factory = KernelFactory::instance(0);
    const std::vector<std::string> opt{"-std=c++17"};
    auto *kBin = factory.create("bin_atoms", kKernels, "batch.cu", opt);
    auto *kPref = factory.create("prefix_cells", kKernels, "batch.cu", opt);
    auto *kScat = factory.create("scatter_atoms", kKernels, "batch.cu", opt);
    auto *kList = factory.create("nlist_cells", kKernels, "batch.cu", opt);
    auto *kQlm = factory.create("qlm_l3", kKernels, "batch.cu", opt);
    auto *kChill = factory.create("chill_from_qlm", kKernels, "batch.cu", opt);

    const int nTot = nAtoms * out.plan.frames;
    const int block = 128;
    const int grid = (nTot + block - 1) / block;
    const double rc2 = rc * rc;
    t0 = std::chrono::steady_clock::now();
    int nA = nAtoms;
    int nF = out.plan.frames;
    int kM = kMax;
    {
      std::vector<void *> a = {&dxyz.p, &dbox.p, &nA, &nF, (void *)&rc,
                               &dncell.p, &dcellOf.p, &dcellCount.p};
      kBin->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&dncell.p, &dcellCount.p, &dcellOff.p, &nA,
                               &nF};
      kPref->launch(dim3((nF + 31) / 32 == 0 ? 1 : (nF + 31) / 32), dim3(32),
                    0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&dcellOf.p, &dcellCount.p, &dcellOff.p,
                               &dorder.p, &nA, &nF};
      kScat->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&dxyz.p, &dbox.p, &dncell.p, &dcellOff.p,
                               &dorder.p, &nA, &nF, (void *)&rc2, &kM,
                               &ddeg.p, &dcols.p};
      kList->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&dxyz.p, &dbox.p, &ddeg.p, &dcols.p, &nA, &nF,
                               &kM, &dqre.p, &dqim.p};
      kQlm->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&ddeg.p, &dcols.p, &dqre.p, &dqim.p, &nA, &nF,
                               &kM, &dchill.p, &dsix.p};
      kChill->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    out.computeMs = msSince(t0);

    out.chill.assign(nN, 0);
    out.sixCount.assign(nN, 0);
    t0 = std::chrono::steady_clock::now();
    checkCuda(rt.cudaMemcpy(out.chill.data(), dchill.p, nN * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH chill");
    checkCuda(rt.cudaMemcpy(out.sixCount.data(), dsix.p, nN * sizeof(int),
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
