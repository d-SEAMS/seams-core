#include <gpu_batch.hpp>

#ifdef SEAMS_HAS_LINKCELL
#include <linkcell_gpu.hpp>
#endif

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

// Resident batch. k-nearest comes from linkcell::gpu (fold, bin,
// Chebyshev shells). Then mutual 4-graph, one q_lm, CHILL, primitive
// six-rings, HC/DDC.
const char *kKernels = R"CUDA(
__device__ inline void minImage(double& dx, double& dy, double& dz,
    double bx, double by, double bz) {
  dx -= bx * nearbyint(dx / bx);
  dy -= by * nearbyint(dy / by);
  dz -= bz * nearbyint(dz / bz);
}

// l=3 Y_lm, m = -3..3, matching seams::steinhardt::ylmAll (phi = atan2(dx, dy)).
__device__ inline void ylm3(double dx, double dy, double dz,
    double* re, double* im) {
  const double r = sqrt(dx * dx + dy * dy + dz * dz);
  if (r < 1.0e-12) return;
  const double z = dz / r;
  const double phi = atan2(dx, dy);
  const double sinT = sqrt(fmax(0.0, 1.0 - z * z));
  const double cosT = z;
  const double cphi = cos(phi);
  const double sphi = sin(phi);
  double s[4], c[4], pr[4], pi[4];
  s[0] = 1.0; c[0] = 1.0; pr[0] = 1.0; pi[0] = 0.0;
  for (int k = 1; k <= 3; ++k) {
    s[k] = s[k - 1] * sinT;
    c[k] = c[k - 1] * cosT;
    const double nr = pr[k - 1] * cphi - pi[k - 1] * sphi;
    const double ni = pr[k - 1] * sphi + pi[k - 1] * cphi;
    pr[k] = nr;
    pi[k] = ni;
  }
  const double piC = 3.14159265358979323846;
  const double amp[4] = {
      0.25 * sqrt(7.0 / piC) * (5.0 * c[3] - 3.0 * c[1]),
      0.125 * sqrt(21.0 / piC) * s[1] * (5.0 * c[2] - 1.0),
      0.25 * sqrt(105.0 / (2.0 * piC)) * s[2] * c[1],
      0.125 * sqrt(35.0 / piC) * s[3]};
  for (int absM = 0; absM <= 3; ++absM) {
    const double a = amp[absM];
    re[3 - absM] += a * pr[absM];
    im[3 - absM] += -a * pi[absM];
    const double sign = (absM % 2 == 0) ? 1.0 : -1.0;
    re[3 + absM] += sign * a * pr[absM];
    im[3 + absM] += sign * a * pi[absM];
  }
}

__device__ inline bool bonded(const int* deg, const int* cols,
    int f, int nAtoms, int kMax, int a, int b) {
  if (a < 0 || b < 0) return false;
  const int d = deg[f * nAtoms + a];
  const int row = (f * nAtoms + a) * kMax;
  for (int t = 0; t < d; ++t) {
    if (cols[row + t] == b) return true;
  }
  return false;
}

__device__ inline bool inSix(const int* r, int atom) {
  for (int t = 0; t < 6; ++t) {
    if (r[t] == atom) return true;
  }
  return false;
}

__device__ inline bool shareAtoms(const int* a, const int* b) {
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i])) return true;
  }
  return false;
}

__device__ inline int commonCount(const int* a, const int* b) {
  int n = 0;
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i])) ++n;
  }
  return n;
}

__device__ inline bool commonInThree(const int* a, const int* b, const int* c) {
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i]) && inSix(c, a[i])) return true;
  }
  return false;
}

__device__ inline bool shareNeigh(const int* deg, const int* cols,
    int f, int nAtoms, int kMax, int a, int b) {
  const int da = deg[f * nAtoms + a];
  const int ra = (f * nAtoms + a) * kMax;
  for (int t = 0; t < da; ++t) {
    if (bonded(deg, cols, f, nAtoms, kMax, b, cols[ra + t])) return true;
  }
  return false;
}

// Franzblau SP on a 6-cycle of the mutual graph: no chords, opposite
// vertices at graph distance 3 (no common neighbour).
__device__ inline bool primitiveSix(const int* r, const int* deg,
    const int* cols, int f, int nAtoms, int kMax) {
  const int d2[6][2] = {{0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}, {5, 1}};
  for (int t = 0; t < 6; ++t) {
    if (bonded(deg, cols, f, nAtoms, kMax, r[d2[t][0]], r[d2[t][1]])) {
      return false;
    }
  }
  const int opp[3][2] = {{0, 3}, {1, 4}, {2, 5}};
  for (int t = 0; t < 3; ++t) {
    const int u = r[opp[t][0]];
    const int v = r[opp[t][1]];
    if (bonded(deg, cols, f, nAtoms, kMax, u, v)) return false;
    if (shareNeigh(deg, cols, f, nAtoms, kMax, u, v)) return false;
  }
  return true;
}

__device__ inline bool basalNeighbours(const int* deg, const int* cols,
    int f, int nAtoms, int kMax, int n1, int n2, int atomOne, int atomTwo) {
  const bool n1one = bonded(deg, cols, f, nAtoms, kMax, atomOne, n1);
  const bool n1two = bonded(deg, cols, f, nAtoms, kMax, atomTwo, n1);
  if (!n1one && !n1two) return false;
  if (n1one) {
    return bonded(deg, cols, f, nAtoms, kMax, atomTwo, n2);
  }
  return bonded(deg, cols, f, nAtoms, kMax, atomOne, n2);
}

__device__ inline bool notNeighboursOfRing(const int* deg, const int* cols,
    int f, int nAtoms, int kMax, const int* trip, const int* ring) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      if (bonded(deg, cols, f, nAtoms, kMax, ring[j], trip[i])) return false;
    }
  }
  return true;
}

__device__ inline bool basalConditions(const int* deg, const int* cols,
    int f, int nAtoms, int kMax, const int* b1, const int* b2) {
  int kIndex = -1;
  int compare1 = 0, compare2 = 0;
  bool l1n = false, l2n = false;
  const int l1 = b1[0];
  const int l2 = b1[1];
  for (int k = 0; k < 6; ++k) {
    const int mk = b2[k];
    if (bonded(deg, cols, f, nAtoms, kMax, l1, mk)) {
      compare1 = b1[2];
      compare2 = b1[4];
      kIndex = k;
      l1n = true;
      break;
    }
    if (bonded(deg, cols, f, nAtoms, kMax, l2, mk)) {
      compare1 = b1[3];
      compare2 = b1[5];
      kIndex = k;
      l2n = true;
      break;
    }
  }
  if (!l1n && !l2n) return false;
  int evenT[3];
  int oddT[3];
  int ie = 0, io = 0;
  for (int k = 0; k <= 5; ++k) {
    int ck = kIndex + k;
    if (ck >= 6) ck -= 6;
    if (k % 2 == 0) evenT[ie++] = b2[ck];
    else oddT[io++] = b2[ck];
  }
  if (!basalNeighbours(deg, cols, f, nAtoms, kMax, evenT[1], evenT[2],
                       compare1, compare2)) {
    return false;
  }
  return notNeighboursOfRing(deg, cols, f, nAtoms, kMax, oddT, b1);
}

__device__ inline int firstRingThrough(const int* A, int nA, const int* B,
    int nB, const int* C, int nC, int skipA, int skipB) {
  int i = 0, j = 0, k = 0;
  while (i < nA && j < nB && k < nC) {
    const int x = A[i], y = B[j], z = C[k];
    if (x == y && y == z) {
      if (x != skipA && x != skipB) return x;
      ++i; ++j; ++k;
      continue;
    }
    int lo = x;
    if (y < lo) lo = y;
    if (z < lo) lo = z;
    if (x == lo) ++i;
    if (y == lo) ++j;
    if (z == lo) ++k;
  }
  return -1;
}

__device__ inline int ringsThrough(const int* A, int nA, const int* B, int nB,
    const int* C, int nC, int skipA, int skipB, int* out, int cap) {
  int i = 0, j = 0, k = 0, n = 0;
  while (i < nA && j < nB && k < nC) {
    const int x = A[i], y = B[j], z = C[k];
    if (x == y && y == z) {
      if (x != skipA && x != skipB && n < cap) out[n++] = x;
      ++i; ++j; ++k;
      continue;
    }
    int lo = x;
    if (y < lo) lo = y;
    if (z < lo) lo = z;
    if (x == lo) ++i;
    if (y == lo) ++j;
    if (z == lo) ++k;
  }
  return n;
}

extern "C" __global__ void mutual_knn(int* deg, const int* cols,
    int nAtoms, int nFrames, int kMax, int* mdeg, int* mcols) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const int row = (f * nAtoms + i) * kMax;
  int found = 0;
  for (int a = 0; a < kMax; ++a) {
    if (cols[row + a] < 0) break;
    ++found;
  }
  deg[f * nAtoms + i] = found;
  const int take = found < 4 ? found : 4;
  const int mrow = (f * nAtoms + i) * 4;
  int kept = 0;
  for (int a = 0; a < take && kept < 4; ++a) {
    const int j = cols[row + a];
    bool back = false;
    const int jrow = (f * nAtoms + j) * kMax;
    for (int t = 0; t < kMax; ++t) {
      if (cols[jrow + t] < 0) break;
      if (cols[jrow + t] == i) {
        back = true;
        break;
      }
    }
    if (back) {
      mcols[mrow + kept] = j;
      ++kept;
    }
  }
  mdeg[f * nAtoms + i] = kept;
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
  double re[7] = {0, 0, 0, 0, 0, 0, 0};
  double im[7] = {0, 0, 0, 0, 0, 0, 0};
  for (int a = 0; a < take; ++a) {
    const int j = cols[row + a];
    double dx = xyz[base + j * 3 + 0] - ix;
    double dy = xyz[base + j * 3 + 1] - iy;
    double dz = xyz[base + j * 3 + 2] - iz;
    minImage(dx, dy, dz, bx, by, bz);
    ylm3(dx, dy, dz, re, im);
  }
  if (take > 0) {
    const double inv = 1.0 / (double)take;
    for (int m = 0; m < 7; ++m) {
      re[m] *= inv;
      im[m] *= inv;
    }
  }
  const int qrow = (f * nAtoms + i) * 7;
  for (int m = 0; m < 7; ++m) {
    qre[qrow + m] = re[m];
    qim[qrow + m] = im[m];
  }
}

extern "C" __global__ void chill_from_qlm(const int* deg, const int* cols,
    const double* qre, const double* qim, int nAtoms, int nFrames, int kMax,
    int* chill) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const int take0 = deg[f * nAtoms + i];
  const int take = take0 < 4 ? take0 : 4;
  const int row = (f * nAtoms + i) * kMax;
  const int qi = (f * nAtoms + i) * 7;
  int S = 0;
  int E = 0;
  for (int a = 0; a < take; ++a) {
    const int j = cols[row + a];
    const int qj = (f * nAtoms + j) * 7;
    double num = 0.0, ni = 0.0, nj = 0.0;
    for (int m = 0; m < 7; ++m) {
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
}

extern "C" __global__ void enum_six(const int* mdeg, const int* mcols,
    int nAtoms, int nFrames, int maxRings, int* nRings, int* ringAtoms,
    int* dropped) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  const int f = tid / nAtoms;
  const int i = tid % nAtoms;
  const int di = mdeg[f * nAtoms + i];
  if (di < 2) return;
  const int irow = (f * nAtoms + i) * 4;
  for (int ia = 0; ia < di; ++ia) {
    const int a = mcols[irow + ia];
    for (int ib = ia + 1; ib < di; ++ib) {
      const int b = mcols[irow + ib];
      const int da = mdeg[f * nAtoms + a];
      const int db = mdeg[f * nAtoms + b];
      const int arow = (f * nAtoms + a) * 4;
      const int brow = (f * nAtoms + b) * 4;
      for (int ix = 0; ix < da; ++ix) {
        const int x = mcols[arow + ix];
        if (x == i || x == b) continue;
        for (int iy = 0; iy < db; ++iy) {
          const int y = mcols[brow + iy];
          if (y == i || y == a || y == x) continue;
          const int dx = mdeg[f * nAtoms + x];
          const int xrow = (f * nAtoms + x) * 4;
          for (int iz = 0; iz < dx; ++iz) {
            const int z = mcols[xrow + iz];
            if (z == i || z == a || z == b || z == y) continue;
            if (!bonded(mdeg, mcols, f, nAtoms, 4, z, y)) continue;
            const int cyc[6] = {i, a, x, z, y, b};
            int mn = i;
            for (int t = 1; t < 6; ++t) {
              if (cyc[t] < mn) mn = cyc[t];
            }
            if (mn != i) continue;
            if (!primitiveSix(cyc, mdeg, mcols, f, nAtoms, 4)) continue;
            const int slot = atomicAdd(nRings + f, 1);
            if (slot >= maxRings) {
              atomicAdd(dropped, 1);
              continue;
            }
            const int dest = (f * maxRings + slot) * 6;
            for (int t = 0; t < 6; ++t) ringAtoms[dest + t] = cyc[t];
          }
        }
      }
    }
  }
}

extern "C" __global__ void invert_rings(const int* nRings, const int* ringAtoms,
    int nAtoms, int nFrames, int maxRings, int maxPer, int* throughCount,
    int* through) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int cap = nFrames * maxRings;
  if (tid >= cap) return;
  const int f = tid / maxRings;
  const int r = tid % maxRings;
  if (r >= nRings[f]) return;
  const int* ring = ringAtoms + (f * maxRings + r) * 6;
  for (int t = 0; t < 6; ++t) {
    const int atom = ring[t];
    if (atom < 0 || atom >= nAtoms) continue;
    const int slot = atomicAdd(throughCount + f * nAtoms + atom, 1);
    if (slot < maxPer) {
      through[(f * nAtoms + atom) * maxPer + slot] = r;
    }
  }
}

extern "C" __global__ void sort_through(int* throughCount, int* through,
    int nAtoms, int nFrames, int maxPer) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  int n = throughCount[tid];
  if (n > maxPer) n = maxPer;
  throughCount[tid] = n;
  int* row = through + tid * maxPer;
  for (int a = 1; a < n; ++a) {
    const int key = row[a];
    int p = a;
    while (p > 0 && row[p - 1] > key) {
      row[p] = row[p - 1];
      --p;
    }
    row[p] = key;
  }
}

// Host cageAffiliation: one ordered basal pair, then its prismatics.
extern "C" __global__ void emit_basal(const int* nRings, const int* ringAtoms,
    const int* mdeg, const int* mcols, const int* throughCount,
    const int* through, int nAtoms, int nFrames, int maxRings, int maxPer,
    int maxPairs, int* nPairs, int* pairs) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int cap = nFrames * maxRings;
  if (tid >= cap) return;
  const int f = tid / maxRings;
  const int i = tid % maxRings;
  if (i >= nRings[f]) return;
  const int* bi = ringAtoms + (f * maxRings + i) * 6;
  for (int slot = 0; slot < 2; ++slot) {
    const int anchor = bi[slot];
    const int da = mdeg[f * nAtoms + anchor];
    const int arow = (f * nAtoms + anchor) * 4;
    for (int n = 0; n < da; ++n) {
      const int nb = mcols[arow + n];
      const int nr = throughCount[f * nAtoms + nb];
      const int* row = through + (f * nAtoms + nb) * maxPer;
      for (int t = 0; t < nr; ++t) {
        const int j = row[t];
        if (j == i) continue;
        const int* bj = ringAtoms + (f * maxRings + j) * 6;
        if (shareAtoms(bi, bj)) continue;
        if (!basalConditions(mdeg, mcols, f, nAtoms, 4, bi, bj)) continue;
        const int slotp = atomicAdd(nPairs + f, 1);
        if (slotp < maxPairs) {
          const int p = (f * maxPairs + slotp) * 2;
          pairs[p] = i;
          pairs[p + 1] = j;
        }
      }
    }
  }
}

extern "C" __global__ void apply_hc(const int* nPairs, const int* pairs,
    const int* ringAtoms, const int* throughCount, const int* through,
    int nAtoms, int nFrames, int maxRings, int maxPer, int maxPairs, int* hc) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int cap = nFrames * maxPairs;
  if (tid >= cap) return;
  const int f = tid / maxPairs;
  const int p = tid % maxPairs;
  if (p >= nPairs[f]) return;
  const int i = pairs[(f * maxPairs + p) * 2];
  const int j = pairs[(f * maxPairs + p) * 2 + 1];
  hc[f * maxRings + i] = 1;
  hc[f * maxRings + j] = 1;
  const int* bi = ringAtoms + (f * maxRings + i) * 6;
  const int* bj = ringAtoms + (f * maxRings + j) * 6;
  for (int q = 0; q < 6; ++q) {
    int trip[3];
    for (int m = 0; m < 3; ++m) trip[m] = bi[(q + m) % 6];
    const int nA = throughCount[f * nAtoms + trip[0]];
    const int nB = throughCount[f * nAtoms + trip[1]];
    const int nC = throughCount[f * nAtoms + trip[2]];
    const int* A = through + (f * nAtoms + trip[0]) * maxPer;
    const int* B = through + (f * nAtoms + trip[1]) * maxPer;
    const int* C = through + (f * nAtoms + trip[2]) * maxPer;
    int cand[16];
    const int nc = ringsThrough(A, nA, B, nB, C, nC, i, j, cand, 16);
    for (int c = 0; c < nc; ++c) {
      const int kr = cand[c];
      const int* bk = ringAtoms + (f * maxRings + kr) * 6;
      int rest[3];
      int nrst = 0;
      for (int u = 0; u < 6; ++u) {
        if (!inSix(trip, bk[u]) && nrst < 3) rest[nrst++] = bk[u];
      }
      if (nrst == 3 && commonCount(rest, bj) == 3) {
        hc[f * maxRings + kr] = 1;
      }
    }
  }
}

extern "C" __global__ void ddc_affil(const int* nRings, const int* ringAtoms,
    const int* throughCount, const int* through, const int* hc,
    int nAtoms, int nFrames, int maxRings, int maxPer, int* ddc) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int cap = nFrames * maxRings;
  if (tid >= cap) return;
  const int f = tid / maxRings;
  const int i = tid % maxRings;
  if (i >= nRings[f]) return;
  if (hc[f * maxRings + i]) return;
  const int* bi = ringAtoms + (f * maxRings + i) * 6;
  int peri[32];
  int nPeri = 0;
  for (int m = 0; m < 6; ++m) {
    const int atom = bi[m];
    const int nr = throughCount[f * nAtoms + atom];
    const int* row = through + (f * nAtoms + atom) * maxPer;
    int common = 0;
    for (int t = 0; t < nr; ++t) {
      if (row[t] == i) continue;
      ++common;
      if (nPeri < 32) peri[nPeri++] = row[t];
    }
    if (common < 3) return;
  }
  int newP[6];
  for (int k = 0; k < 6; ++k) {
    int trip[3];
    for (int t = 0; t < 3; ++t) trip[t] = bi[(k + t) % 6];
    const int nA = throughCount[f * nAtoms + trip[0]];
    const int nB = throughCount[f * nAtoms + trip[1]];
    const int nC = throughCount[f * nAtoms + trip[2]];
    const int* A = through + (f * nAtoms + trip[0]) * maxPer;
    const int* B = through + (f * nAtoms + trip[1]) * maxPer;
    const int* C = through + (f * nAtoms + trip[2]) * maxPer;
    const int j = firstRingThrough(A, nA, B, nB, C, nC, i, -1);
    if (j < 0) return;
    newP[k] = j;
  }
  const int* p0 = ringAtoms + (f * maxRings + newP[0]) * 6;
  const int* p1 = ringAtoms + (f * maxRings + newP[1]) * 6;
  const int* p2 = ringAtoms + (f * maxRings + newP[2]) * 6;
  const int* p3 = ringAtoms + (f * maxRings + newP[3]) * 6;
  const int* p4 = ringAtoms + (f * maxRings + newP[4]) * 6;
  const int* p5 = ringAtoms + (f * maxRings + newP[5]) * 6;
  if (!commonInThree(p0, p2, p4)) return;
  if (!commonInThree(p1, p3, p5)) return;
  const int* pairs[4][2] = {{p0, p2}, {p1, p3}, {p2, p4}, {p3, p5}};
  for (int t = 0; t < 4; ++t) {
    if (commonCount(pairs[t][0], pairs[t][1]) < 3) return;
  }
  ddc[f * maxRings + i] = 1;
  for (int t = 0; t < 6; ++t) ddc[f * maxRings + newP[t]] = 1;
}

extern "C" __global__ void atom_ice(const int* nRings, const int* ringAtoms,
    const int* hc, const int* ddc, int nAtoms, int nFrames, int maxRings,
    int* atomHc, int* atomDdc, int* sixCount) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  const int cap = nFrames * maxRings;
  if (tid >= cap) return;
  const int f = tid / maxRings;
  const int r = tid % maxRings;
  if (r >= nRings[f]) return;
  const int* ring = ringAtoms + (f * maxRings + r) * 6;
  const int isHc = hc[f * maxRings + r];
  const int isDdc = ddc[f * maxRings + r];
  for (int t = 0; t < 6; ++t) {
    const int a = ring[t];
    if (a < 0 || a >= nAtoms) continue;
    atomicAdd(sixCount + f * nAtoms + a, 1);
    if (isHc) atomHc[f * nAtoms + a] = 1;
    if (isDdc) atomDdc[f * nAtoms + a] = 1;
  }
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
  void reset() {
    if (p != nullptr && CUDART::loaded()) {
      CUDART::instance().cudaFree(p);
    }
    p = nullptr;
  }
};

// LAMMPS DualView / vesin CudaNeighborListExtras: grow once, reuse.
struct Workspace {
  int capN = 0;
  int capF = 0;
  int capK = 0;
  int capPer = 0;
  DevPtr dxyz, dbox;
  DevPtr ddeg, dcols, dmdeg, dmcols, dqre, dqim, dchill;
  DevPtr dnRings, dringAtoms, ddropped, dthroughCount, dthrough;
  DevPtr dhc, dddc, datomHc, datomDdc, dsix, dnPairs, dpairs;

  void growOne(DevPtr &d, std::size_t bytes, const char *what) {
    d.reset();
    checkCuda(CUDART::instance().cudaMalloc(&d.p, bytes), what);
  }

  void ensure(int nAtoms, int nF, int kMax, int maxPer) {
    if (capN == nAtoms && capF == nF && capK == kMax && capPer == maxPer &&
        capN > 0) {
      return;
    }
    // Exact extents so kernel strides (nAtoms+1, nAtoms) match the allocation.
    // Reuse across calls of the same batch shape (LAMMPS DualView grow).
    capN = nAtoms;
    capF = nF;
    capK = kMax;
    capPer = maxPer;
    const std::size_t n = static_cast<std::size_t>(capN);
    const std::size_t nf = static_cast<std::size_t>(capF);
    const std::size_t nN = nf * n;
    const std::size_t maxRings = n * static_cast<std::size_t>(capPer);
    growOne(dxyz, nN * 3 * sizeof(double), "xyz");
    growOne(dbox, nf * 3 * sizeof(double), "box");
    growOne(ddeg, nN * sizeof(int), "deg");
    growOne(dcols, nN * static_cast<std::size_t>(capK) * sizeof(int), "cols");
    growOne(dmdeg, nN * sizeof(int), "mdeg");
    growOne(dmcols, nN * 4 * sizeof(int), "mcols");
    growOne(dqre, nN * 7 * sizeof(double), "qre");
    growOne(dqim, nN * 7 * sizeof(double), "qim");
    growOne(dchill, nN * sizeof(int), "chill");
    growOne(dnRings, nf * sizeof(int), "nRings");
    growOne(dringAtoms, nf * maxRings * 6 * sizeof(int), "rings");
    growOne(ddropped, sizeof(int), "dropped");
    growOne(dthroughCount, nN * sizeof(int), "tcount");
    growOne(dthrough, nN * static_cast<std::size_t>(capPer) * sizeof(int),
            "through");
    growOne(dhc, nf * maxRings * sizeof(int), "hc");
    growOne(dddc, nf * maxRings * sizeof(int), "ddc");
    const std::size_t maxPairs = maxRings * 8;
    growOne(dnPairs, nf * sizeof(int), "nPairs");
    growOne(dpairs, nf * maxPairs * 2 * sizeof(int), "pairs");
    growOne(datomHc, nN * sizeof(int), "atomHc");
    growOne(datomDdc, nN * sizeof(int), "atomDdc");
    growOne(dsix, nN * sizeof(int), "six");
  }
};

Workspace &workspace() {
  static Workspace w;
  return w;
}

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
#if !defined(SEAMS_HAS_GPULITE)
  out.error = "built without gpulite";
  return out;
#elif !defined(SEAMS_HAS_LINKCELL) || !defined(LINKCELL_HAS_GPULITE)
  out.error = "linkcell device k-nearest not built";
  return out;
#else
  if (!out.plan.resident) {
    out.error = out.plan.reason;
    return out;
  }
  if (!NVRTC::loaded() || !linkcell::gpu::available()) {
    out.error = "nvrtc or linkcell device path not loaded";
    return out;
  }
  try {
    auto &rt = CUDART::instance();
    const int kMax = out.plan.foot.kMax;
    const int maxPer = out.plan.foot.maxSixRingsPerAtom;
    const int maxRings = nAtoms * maxPer;
    const std::size_t n = static_cast<std::size_t>(nAtoms);
    const std::size_t nf = static_cast<std::size_t>(out.plan.frames);
    const std::size_t xyzN = nf * n * 3;
    const std::size_t boxN = nf * 3;
    const std::size_t nN = nf * n;
    const std::size_t ringN = nf * static_cast<std::size_t>(maxRings);

    auto &ws = workspace();
    ws.ensure(nAtoms, out.plan.frames, kMax, maxPer);
    auto t0 = std::chrono::steady_clock::now();
    checkCuda(rt.cudaMemset(ws.dnRings.p, 0, nf * sizeof(int)), "zero nRings");
    checkCuda(rt.cudaMemset(ws.ddropped.p, 0, sizeof(int)), "zero dropped");
    checkCuda(rt.cudaMemset(ws.dthroughCount.p, 0, nN * sizeof(int)), "zero thru");
    checkCuda(rt.cudaMemset(ws.dhc.p, 0, ringN * sizeof(int)), "zero hc");
    checkCuda(rt.cudaMemset(ws.dnPairs.p, 0, nf * sizeof(int)), "zero nPairs");
    checkCuda(rt.cudaMemset(ws.dddc.p, 0, ringN * sizeof(int)), "zero ddc");
    checkCuda(rt.cudaMemset(ws.datomHc.p, 0, nN * sizeof(int)), "zero atomHc");
    checkCuda(rt.cudaMemset(ws.datomDdc.p, 0, nN * sizeof(int)), "zero atomDdc");
    checkCuda(rt.cudaMemset(ws.dsix.p, 0, nN * sizeof(int)), "zero six");
    checkCuda(rt.cudaMemcpy(ws.dxyz.p, xyz, xyzN * sizeof(double),
                            cudaMemcpyHostToDevice),
              "HtoD xyz");
    checkCuda(rt.cudaMemcpy(ws.dbox.p, box, boxN * sizeof(double),
                            cudaMemcpyHostToDevice),
              "HtoD box");
    out.uploadMs = msSince(t0);

    auto &factory = KernelFactory::instance(0);
    const std::vector<std::string> opt{"-std=c++17"};
    auto *kMut = factory.create("mutual_knn", kKernels, "batch.cu", opt);
    auto *kQlm = factory.create("qlm_l3", kKernels, "batch.cu", opt);
    auto *kChill = factory.create("chill_from_qlm", kKernels, "batch.cu", opt);
    auto *kSix = factory.create("enum_six", kKernels, "batch.cu", opt);
    auto *kInv = factory.create("invert_rings", kKernels, "batch.cu", opt);
    auto *kSort = factory.create("sort_through", kKernels, "batch.cu", opt);
    auto *kEmit = factory.create("emit_basal", kKernels, "batch.cu", opt);
    auto *kHc = factory.create("apply_hc", kKernels, "batch.cu", opt);
    auto *kDdc = factory.create("ddc_affil", kKernels, "batch.cu", opt);
    auto *kAtom = factory.create("atom_ice", kKernels, "batch.cu", opt);

    const int nTot = nAtoms * out.plan.frames;
    const int block = 128;
    const int grid = (nTot + block - 1) / block;
    const int ringTot = out.plan.frames * maxRings;
    const int ringGrid = (ringTot + block - 1) / block;
    t0 = std::chrono::steady_clock::now();
    int nA = nAtoms;
    int nF = out.plan.frames;
    int kM = kMax;
    int mR = maxRings;
    int mP = maxPer;
    {
      static linkcell::gpu::Workspace knn;
      auto *xyzDev = static_cast<double *>(ws.dxyz.p);
      auto *colsDev = static_cast<int *>(ws.dcols.p);
      const std::size_t kSz = static_cast<std::size_t>(kMax);
      const std::size_t nSz = static_cast<std::size_t>(nAtoms);
      const std::size_t fSz = static_cast<std::size_t>(nF);
      const linkcell::Cell cell0 = linkcell::Cell::ortho(
          box[0], box[1], box[2]);
      bool sameBox = true;
      for (int f = 1; f < nF; ++f) {
        if (box[static_cast<std::size_t>(f) * 3 + 0] != box[0] ||
            box[static_cast<std::size_t>(f) * 3 + 1] != box[1] ||
            box[static_cast<std::size_t>(f) * 3 + 2] != box[2]) {
          sameBox = false;
          break;
        }
      }
      if (sameBox) {
        knn.knearest_into_many(xyzDev, nSz, fSz, cell0, kSz, colsDev,
                               fSz * nSz * kSz, nullptr, rc);
      } else {
        for (int f = 0; f < nF; ++f) {
          const linkcell::Cell cell = linkcell::Cell::ortho(
              box[static_cast<std::size_t>(f) * 3 + 0],
              box[static_cast<std::size_t>(f) * 3 + 1],
              box[static_cast<std::size_t>(f) * 3 + 2]);
          knn.knearest_into(xyzDev + static_cast<std::size_t>(f) * nSz * 3, nSz,
                            cell, kSz,
                            colsDev + static_cast<std::size_t>(f) * nSz * kSz,
                            nSz * kSz, nullptr, rc);
        }
      }
    }
    {
      std::vector<void *> a = {&ws.ddeg.p, &ws.dcols.p, &nA, &nF, &kM,
                               &ws.dmdeg.p, &ws.dmcols.p};
      kMut->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&ws.dxyz.p, &ws.dbox.p, &ws.ddeg.p, &ws.dcols.p,
                               &nA, &nF, &kM, &ws.dqre.p, &ws.dqim.p};
      kQlm->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&ws.ddeg.p, &ws.dcols.p, &ws.dqre.p, &ws.dqim.p,
                               &nA, &nF, &kM, &ws.dchill.p};
      kChill->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&ws.dmdeg.p, &ws.dmcols.p, &nA, &nF, &mR,
                               &ws.dnRings.p, &ws.dringAtoms.p, &ws.ddropped.p};
      kSix->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&ws.dnRings.p, &ws.dringAtoms.p, &nA, &nF, &mR,
                               &mP, &ws.dthroughCount.p, &ws.dthrough.p};
      kInv->launch(dim3(ringGrid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&ws.dthroughCount.p, &ws.dthrough.p, &nA, &nF,
                               &mP};
      kSort->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
    }
    {
      int maxPairs = maxRings * 8;
      const int pairGrid = (nF * maxPairs + block - 1) / block;
      std::vector<void *> a = {&ws.dnRings.p, &ws.dringAtoms.p, &ws.dmdeg.p,
                               &ws.dmcols.p, &ws.dthroughCount.p, &ws.dthrough.p,
                               &nA, &nF, &mR, &mP, &maxPairs, &ws.dnPairs.p,
                               &ws.dpairs.p};
      kEmit->launch(dim3(ringGrid), dim3(block), 0, nullptr, a, true);
      std::vector<void *> b = {&ws.dnPairs.p, &ws.dpairs.p, &ws.dringAtoms.p,
                               &ws.dthroughCount.p, &ws.dthrough.p, &nA, &nF,
                               &mR, &mP, &maxPairs, &ws.dhc.p};
      kHc->launch(dim3(pairGrid), dim3(block), 0, nullptr, b, true);
    }
    {
      std::vector<void *> a = {&ws.dnRings.p, &ws.dringAtoms.p,
                               &ws.dthroughCount.p, &ws.dthrough.p, &ws.dhc.p,
                               &nA, &nF, &mR, &mP, &ws.dddc.p};
      kDdc->launch(dim3(ringGrid), dim3(block), 0, nullptr, a, true);
    }
    {
      std::vector<void *> a = {&ws.dnRings.p, &ws.dringAtoms.p, &ws.dhc.p,
                               &ws.dddc.p, &nA, &nF, &mR, &ws.datomHc.p,
                               &ws.datomDdc.p, &ws.dsix.p};
      kAtom->launch(dim3(ringGrid), dim3(block), 0, nullptr, a, true);
    }
    out.computeMs = msSince(t0);

    out.chill.assign(nN, 0);
    out.sixCount.assign(nN, 0);
    out.atomHc.assign(nN, 0);
    out.atomDdc.assign(nN, 0);
    out.nRings.assign(nf, 0);
    t0 = std::chrono::steady_clock::now();
    checkCuda(rt.cudaMemcpy(out.chill.data(), ws.dchill.p, nN * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH chill");
    checkCuda(rt.cudaMemcpy(out.sixCount.data(), ws.dsix.p, nN * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH six");
    checkCuda(rt.cudaMemcpy(out.atomHc.data(), ws.datomHc.p, nN * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH hc");
    checkCuda(rt.cudaMemcpy(out.atomDdc.data(), ws.datomDdc.p, nN * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH ddc");
    checkCuda(rt.cudaMemcpy(out.nRings.data(), ws.dnRings.p, nf * sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH nRings");
    checkCuda(rt.cudaMemcpy(&out.ringsDropped, ws.ddropped.p, sizeof(int),
                            cudaMemcpyDeviceToHost),
              "DtoH dropped");
    for (int &v : out.nRings) {
      if (v > maxRings) {
        v = maxRings;
      }
    }
    out.downloadMs = msSince(t0);
  } catch (const std::exception &ex) {
    out.error = ex.what();
  }
  return out;
#endif
}

} // namespace gpu
