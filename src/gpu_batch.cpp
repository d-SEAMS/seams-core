#include <gpu_batch.hpp>
#include <seams_config.hpp>

#ifdef SEAMS_HAS_LINKCELL
#include <linkcell_gpu.hpp>
#include <neighbours.hpp>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#ifdef SEAMS_HAS_GPULITE
#include <gpulite/gpulite.hpp>
#include <tum_offload.hpp>
using gpulite::CUDART;
using gpulite::KernelFactory;
using gpulite::NVRTC;
#endif

namespace gpu {
namespace {

#ifdef SEAMS_HAS_GPULITE

// kRowCap of the kernel source below
constexpr int kRowCap = 32;

// TUM ice score only: linkcell k-NN, mutual 4-graph, primitive
// six-rings, HC/DDC affiliation. No CHILL, no q_lm.
const char *kKernels = R"CUDA(
// Most ring ids per atom row the local buffers hold; ddc_affil keeps one
// bit per row entry in an unsigned.
constexpr int kRowCap = 32;

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

__device__ inline bool notNeighboursOfRing(const int* deg, const int* cols,
    int f, int nAtoms, int kMax, const int* trip, const int* ring) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      if (bonded(deg, cols, f, nAtoms, kMax, ring[j], trip[i])) return false;
    }
  }
  return true;
}

// beta(b1, b2) for rings sharing no atom: some alternating triplet of b1
// bonds one-to-one onto an alternating triplet {m_q, m_q+2, m_q+4} of b2,
// and the other three atoms of b2 bond to nothing in b1. All 2 x 2 x 6
// alignments, as tum::device::basalConditions.
__device__ inline bool basalConditions(const int* deg, const int* cols,
    int f, int nAtoms, int kMax, const int* b1, const int* b2) {
  const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                          {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (int q = 0; q < 2; ++q) {
    const int oddT[3] = {b2[q + 1], b2[q + 3], b2[(q + 5) % 6]};
    if (!notNeighboursOfRing(deg, cols, f, nAtoms, kMax, oddT, b1)) continue;
    for (int p = 0; p < 2; ++p) {
      for (int s = 0; s < 6; ++s) {
        bool all = true;
        for (int t = 0; t < 3 && all; ++t) {
          all = bonded(deg, cols, f, nAtoms, kMax, b1[p + 2 * t],
                       b2[q + 2 * perm[s][t]]);
        }
        if (all) return true;
      }
    }
  }
  return false;
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

extern "C" __global__ void mutual_knn(const int* cols,
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
  const int take = found < 4 ? found : 4;
  const int mrow = (f * nAtoms + i) * 4;
  int kept = 0;
  for (int a = 0; a < take && kept < 4; ++a) {
    const int j = cols[row + a];
    bool back = false;
    const int jrow = (f * nAtoms + j) * kMax;
    for (int t = 0; t < 4; ++t) {
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

// nRings counts every ring of the frame; only the first maxRings are stored
extern "C" __global__ void enum_six(const int* mdeg, const int* mcols,
    int nAtoms, int nFrames, int maxRings, int* nRings, int* ringAtoms) {
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
            // Slot order depends on scheduling; no HC or DDC test reads it.
            const int slot = atomicAdd(nRings + f, 1);
            if (slot >= maxRings) continue;
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

// rowMax[f] is the longest row of frame f before the clamp
extern "C" __global__ void sort_through(int* throughCount, int* through,
    int nAtoms, int nFrames, int maxPer, int* rowMax) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= nAtoms * nFrames) return;
  int n = throughCount[tid];
  atomicMax(rowMax + tid / nAtoms, n);
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

// True when ring j is in the through row of a bond neighbour of bi[0] or
// bi[1] visited before (slot, n), so each candidate is tested once.
__device__ inline bool seenEarlier(const int* bi, const int* mdeg,
    const int* mcols, const int* throughCount, const int* through, int f,
    int nAtoms, int maxPer, int slot, int n, int j) {
  for (int s = 0; s <= slot; ++s) {
    const int anchor = bi[s];
    const int arow = (f * nAtoms + anchor) * 4;
    const int stop = s < slot ? mdeg[f * nAtoms + anchor] : n;
    for (int m = 0; m < stop; ++m) {
      const int nb = mcols[arow + m];
      const int nr = throughCount[f * nAtoms + nb];
      const int* row = through + (f * nAtoms + nb) * maxPer;
      for (int t = 0; t < nr; ++t) {
        if (row[t] == j) return true;
      }
    }
  }
  return false;
}

// Host cageAffiliation: ordered basal pairs, then their prismatics. Every
// pair passing beta has a bond from bi[0] or bi[1] into bj, so the
// neighbours of those two atoms reach every candidate.
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
        if (seenEarlier(bi, mdeg, mcols, throughCount, through, f, nAtoms,
                        maxPer, slot, n, j)) {
          continue;
        }
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
    int cand[kRowCap];
    const int nc = ringsThrough(A, nA, B, nB, C, nC, i, j, cand, kRowCap);
    for (int c = 0; c < nc; ++c) {
      const int kr = cand[c];
      const int* bk = ringAtoms + (f * maxRings + kr) * 6;
      int rest = 0, restInJ = 0;
      for (int u = 0; u < 6; ++u) {
        const int a = bk[u];
        if (a != trip[0] && a != trip[1] && a != trip[2]) {
          ++rest;
          if (inSix(bj, a)) ++restInJ;
        }
      }
      if (rest == 3 && restInJ == 3) {
        hc[f * maxRings + kr] = 1;
      }
    }
  }
}

// epsilon(i) for a ring that is not HC, as tum::device::ddcFrom: every atom
// in 3 other rings, and pairwise distinct P_k in T_k with P_0, P_2, P_4
// sharing an atom, P_1, P_3, P_5 sharing an atom and, for some starting atom
// and direction of i, the alternate pairs sharing 3 atoms. On a pass marks i
// and every P_k of every valid assignment.
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
  const int* frameRings = ringAtoms + f * maxRings * 6;
  for (int m = 0; m < 6; ++m) {
    const int atom = bi[m];
    const int nr = throughCount[f * nAtoms + atom];
    const int* row = through + (f * nAtoms + atom) * maxPer;
    int common = 0;
    for (int t = 0; t < nr; ++t) {
      if (row[t] != i) ++common;
    }
    if (common < 3) return;
  }
  // T[k]: every ring but i through bi[k], bi[k+1], bi[k+2]
  int T[6][kRowCap];
  int nT[6];
  for (int k = 0; k < 6; ++k) {
    int trip[3];
    for (int t = 0; t < 3; ++t) trip[t] = bi[(k + t) % 6];
    const int nA = throughCount[f * nAtoms + trip[0]];
    const int nB = throughCount[f * nAtoms + trip[1]];
    const int nC = throughCount[f * nAtoms + trip[2]];
    const int* A = through + (f * nAtoms + trip[0]) * maxPer;
    const int* B = through + (f * nAtoms + trip[1]) * maxPer;
    const int* C = through + (f * nAtoms + trip[2]) * maxPer;
    nT[k] = ringsThrough(A, nA, B, nB, C, nC, i, -1, T[k], kRowCap);
    if (nT[k] == 0) return;
  }
  // Bit t of used[k] marks T[k][t] as P_k of some valid assignment
  unsigned used[6] = {0u, 0u, 0u, 0u, 0u, 0u};
  bool any = false;
  for (int a = 0; a < nT[0]; ++a) {
    const int r0 = T[0][a];
    const int *p0 = frameRings + r0 * 6;
    for (int c = 0; c < nT[2]; ++c) {
      const int r2 = T[2][c];
      const int *p2 = frameRings + r2 * 6;
      if (r2 == r0) {
        continue;
      }
      const bool e02 = commonCount(p0, p2) >= 3;
      for (int e = 0; e < nT[4]; ++e) {
        const int r4 = T[4][e];
        const int *p4 = frameRings + r4 * 6;
        if (r4 == r0 || r4 == r2 || !commonInThree(p0, p2, p4)) {
          continue;
        }
        const bool e24 = commonCount(p2, p4) >= 3;
        const bool e40 = commonCount(p4, p0) >= 3;
        // Bit t: the ring on triplet 2t shares 3 atoms with both others
        const unsigned mE = (e40 && e02 ? 1u : 0u) |
                            (e02 && e24 ? 2u : 0u) |
                            (e24 && e40 ? 4u : 0u);
        if (mE == 0u) {
          continue;
        }
        // Some start and direction of i puts an even middle 2t next to an
        // odd middle 2t+1 or 2t-1
        unsigned partners = 0u;
        for (int t = 0; t < 3; ++t) {
          if ((mE >> t) & 1u) {
            partners |= (1u << t) | (1u << ((t + 2) % 3));
          }
        }
        for (int b = 0; b < nT[1]; ++b) {
          const int r1 = T[1][b];
          const int *p1 = frameRings + r1 * 6;
          if (r1 == r0 || r1 == r2 || r1 == r4) {
            continue;
          }
          for (int d = 0; d < nT[3]; ++d) {
            const int r3 = T[3][d];
            const int *p3 = frameRings + r3 * 6;
            if (r3 == r0 || r3 == r2 || r3 == r4 || r3 == r1) {
              continue;
            }
            const bool o13 = commonCount(p1, p3) >= 3;
            for (int g = 0; g < nT[5]; ++g) {
              const int r5 = T[5][g];
              const int *p5 = frameRings + r5 * 6;
              if (r5 == r0 || r5 == r2 || r5 == r4 || r5 == r1 || r5 == r3 ||
                  !commonInThree(p1, p3, p5)) {
                continue;
              }
              const bool o35 = commonCount(p3, p5) >= 3;
              const bool o51 = commonCount(p5, p1) >= 3;
              const unsigned mO = (o51 && o13 ? 1u : 0u) |
                                  (o13 && o35 ? 2u : 0u) |
                                  (o35 && o51 ? 4u : 0u);
              if ((partners & mO) == 0u) {
                continue;
              }
              any = true;
              used[0] |= 1u << a;
              used[1] |= 1u << b;
              used[2] |= 1u << c;
              used[3] |= 1u << d;
              used[4] |= 1u << e;
              used[5] |= 1u << g;
            }
          }
        }
      }
    }
  }
  if (!any) return;
  ddc[f * maxRings + i] = 1;
  for (int k = 0; k < 6; ++k) {
    for (int t = 0; t < nT[k]; ++t) {
      if ((used[k] >> t) & 1u) ddc[f * maxRings + T[k][t]] = 1;
    }
  }
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
  int capRings = 0;
  int capPairs = 0;
  DevPtr dxyz;
  DevPtr dcols, dmdeg, dmcols;
  DevPtr dnRings, dringAtoms, drowMax, dthroughCount, dthrough;
  DevPtr dhc, dddc, datomHc, datomDdc, dsix, dnPairs, dpairs;

  void growOne(DevPtr &d, std::size_t bytes, const char *what) {
    d.reset();
    checkCuda(CUDART::instance().cudaMalloc(&d.p, bytes), what);
  }

  // maxRings and maxPairs are per frame. A failed allocation leaves the
  // workspace empty, so the next call allocates afresh.
  void ensure(int nAtoms, int nF, int kMax, int maxPer, int maxRings,
              int maxPairs) {
    if (capN == nAtoms && capF == nF && capK == kMax && capPer == maxPer &&
        capRings == maxRings && capPairs == maxPairs && capN > 0) {
      return;
    }
    // Exact extents so kernel strides (nAtoms+1, nAtoms) match the allocation.
    // Reuse across calls of the same batch shape (LAMMPS DualView grow).
    capN = 0;
    const std::size_t n = static_cast<std::size_t>(nAtoms);
    const std::size_t nf = static_cast<std::size_t>(nF);
    const std::size_t nN = nf * n;
    const std::size_t rings = nf * static_cast<std::size_t>(maxRings);
    const std::size_t pairs = nf * static_cast<std::size_t>(maxPairs);
    growOne(dxyz, nN * 3 * sizeof(double), "xyz");
    growOne(dcols, nN * static_cast<std::size_t>(kMax) * sizeof(int), "cols");
    growOne(dmdeg, nN * sizeof(int), "mdeg");
    growOne(dmcols, nN * 4 * sizeof(int), "mcols");
    growOne(dnRings, nf * sizeof(int), "nRings");
    growOne(dringAtoms, rings * 6 * sizeof(int), "rings");
    growOne(drowMax, nf * sizeof(int), "rowMax");
    growOne(dthroughCount, nN * sizeof(int), "tcount");
    growOne(dthrough, nN * static_cast<std::size_t>(maxPer) * sizeof(int),
            "through");
    growOne(dhc, rings * sizeof(int), "hc");
    growOne(dddc, rings * sizeof(int), "ddc");
    growOne(dnPairs, nf * sizeof(int), "nPairs");
    growOne(dpairs, pairs * 2 * sizeof(int), "pairs");
    growOne(datomHc, nN * sizeof(int), "atomHc");
    growOne(datomDdc, nN * sizeof(int), "atomDdc");
    growOne(dsix, nN * sizeof(int), "six");
    capN = nAtoms;
    capF = nF;
    capK = kMax;
    capPer = maxPer;
    capRings = maxRings;
    capPairs = maxPairs;
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
                            int nFrames, double rc, const double *boxLow,
                            int nBox) {
  BatchResult out;
  const auto cfg = seams::cfg::load();
  out.plan = planBatch(nAtoms, nFrames, 16, 16, cfg.resident);
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
    const int nF = out.plan.frames;
    const std::size_t n = static_cast<std::size_t>(nAtoms);
    const std::size_t nf = static_cast<std::size_t>(nF);
    const std::size_t xyzN = nf * n * 3;
    const std::size_t nN = nf * n;
    // The planned sizes; a frame that overflows one reruns with it grown,
    // and a frame with a row past kRowCap is evaluated on the host.
    int maxPer = std::min(out.plan.foot.maxSixRingsPerAtom, kRowCap);
    int maxRings = nAtoms * out.plan.foot.maxSixRingsPerAtom;
    int maxPairs = maxRings * 8;

    auto &factory = KernelFactory::instance(0);
    const std::vector<std::string> opt{"-std=c++17"};
    auto *kMut = factory.create("mutual_knn", kKernels, "batch.cu", opt);
    auto *kSix = factory.create("enum_six", kKernels, "batch.cu", opt);
    auto *kInv = factory.create("invert_rings", kKernels, "batch.cu", opt);
    auto *kSort = factory.create("sort_through", kKernels, "batch.cu", opt);
    auto *kEmit = factory.create("emit_basal", kKernels, "batch.cu", opt);
    auto *kHc = factory.create("apply_hc", kKernels, "batch.cu", opt);
    auto *kDdc = factory.create("ddc_affil", kKernels, "batch.cu", opt);
    auto *kAtom = factory.create("atom_ice", kKernels, "batch.cu", opt);
    const int nTot = nAtoms * nF;
    const int block = cfg.block > 0 ? cfg.block : 128;
    const int grid = (nTot + block - 1) / block;
    static linkcell::gpu::Workspace knn;
    auto launchArgs = [](void **raw, std::size_t m) {
      return std::vector<void *>(raw, raw + m);
    };

    out.sixCount.assign(nN, 0);
    out.atomHc.assign(nN, 0);
    out.atomDdc.assign(nN, 0);
    out.nRings.assign(nf, 0);
    std::vector<int> nRings(nf), nPairs(nf), rowMax(nf);
    std::vector<int> six(nN), hcAtom(nN), ddcAtom(nN);
    std::vector<char> done(nf, 0), onHost(nf, 0);
    // Mutual graph of each frame still pending after a pass, kept on the
    // host in case the next allocation fails
    std::vector<std::vector<int>> hostDeg(nf), hostCols(nf);
    auto &ws = workspace();
    for (;;) {
      if (out.retries == 0) {
        ws.ensure(nAtoms, nF, kMax, maxPer, maxRings, maxPairs);
      } else {
        try {
          ws.ensure(nAtoms, nF, kMax, maxPer, maxRings, maxPairs);
        } catch (const std::runtime_error &) {
          // No room for the grown buffers: the pending frames go to the host
          for (std::size_t f = 0; f < nf; ++f) {
            onHost[f] = done[f] ? 0 : 1;
          }
          break;
        }
      }
      const std::size_t ringN = nf * static_cast<std::size_t>(maxRings);
      auto t0 = std::chrono::steady_clock::now();
      checkCuda(rt.cudaMemset(ws.dnRings.p, 0, nf * sizeof(int)), "zero nRings");
      checkCuda(rt.cudaMemset(ws.drowMax.p, 0, nf * sizeof(int)), "zero rowMax");
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
      checkCuda(rt.cudaDeviceSynchronize(), "upload");
      out.uploadMs += msSince(t0);

      const int ringTot = nF * maxRings;
      const int ringGrid = (ringTot + block - 1) / block;
      const int pairGrid = (nF * maxPairs + block - 1) / block;
      t0 = std::chrono::steady_clock::now();
      int nA = nAtoms;
      int nFr = nF;
      int kM = kMax;
      int mR = maxRings;
      int mP = maxPer;
      int mPairs = maxPairs;
      auto *xyzDev = static_cast<double *>(ws.dxyz.p);
      auto *colsDev = static_cast<int *>(ws.dcols.p);
      const std::size_t kSz = static_cast<std::size_t>(kMax);
      linkcell::Cell cell0 = linkcell::Cell::ortho(box[0], box[1], box[2]);
      const double *frameLens = box;
      nneigh::residentFrameCell(box, boxLow, nBox, cell0, frameLens);
      void *q = knn.queue();
      // One multi-frame launch on the workspace stream. Per-frame box
      // walks serialize on the same bin buffers and leave the SMs idle.
      knn.knearest_into_many(xyzDev, n, nf, cell0, kSz, colsDev, nf * n * kSz,
                             nullptr, cfg.cell, false, frameLens);
      {
        void *raw[] = {&ws.dcols.p, &nA, &nFr, &kM, &ws.dmdeg.p, &ws.dmcols.p};
        auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
        kMut->launch(dim3(static_cast<unsigned>(grid)), dim3(static_cast<unsigned>(block)), 0, q, a, false);
      }
      {
        void *raw[] = {&ws.dmdeg.p, &ws.dmcols.p, &nA, &nFr, &mR,
                       &ws.dnRings.p, &ws.dringAtoms.p};
        auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
        kSix->launch(dim3(static_cast<unsigned>(grid)), dim3(static_cast<unsigned>(block)), 0, q, a, false);
      }
      {
        void *raw[] = {&ws.dnRings.p, &ws.dringAtoms.p, &nA, &nFr, &mR, &mP,
                       &ws.dthroughCount.p, &ws.dthrough.p};
        auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
        kInv->launch(dim3(static_cast<unsigned>(ringGrid)), dim3(static_cast<unsigned>(block)), 0, q, a, false);
      }
      {
        void *raw[] = {&ws.dthroughCount.p, &ws.dthrough.p, &nA, &nFr, &mP,
                       &ws.drowMax.p};
        auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
        kSort->launch(dim3(static_cast<unsigned>(grid)), dim3(static_cast<unsigned>(block)), 0, q, a, false);
      }
      {
        void *emitRaw[] = {&ws.dnRings.p,        &ws.dringAtoms.p, &ws.dmdeg.p,
                           &ws.dmcols.p,         &ws.dthroughCount.p,
                           &ws.dthrough.p,       &nA,              &nFr,
                           &mR,                  &mP,              &mPairs,
                           &ws.dnPairs.p,        &ws.dpairs.p};
        auto a = launchArgs(emitRaw, sizeof(emitRaw) / sizeof(emitRaw[0]));
        kEmit->launch(dim3(static_cast<unsigned>(ringGrid)), dim3(static_cast<unsigned>(block)), 0, q, a, false);
        void *hcRaw[] = {&ws.dnPairs.p,  &ws.dpairs.p,        &ws.dringAtoms.p,
                         &ws.dthroughCount.p, &ws.dthrough.p, &nA,
                         &nFr,           &mR,                 &mP,
                         &mPairs,        &ws.dhc.p};
        auto b = launchArgs(hcRaw, sizeof(hcRaw) / sizeof(hcRaw[0]));
        kHc->launch(dim3(static_cast<unsigned>(pairGrid)), dim3(static_cast<unsigned>(block)), 0, q, b, false);
      }
      {
        void *raw[] = {&ws.dnRings.p,        &ws.dringAtoms.p, &ws.dthroughCount.p,
                       &ws.dthrough.p,       &ws.dhc.p,        &nA,
                       &nFr,                 &mR,              &mP,
                       &ws.dddc.p};
        auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
        kDdc->launch(dim3(static_cast<unsigned>(ringGrid)), dim3(static_cast<unsigned>(block)), 0, q, a, false);
      }
      {
        void *raw[] = {&ws.dnRings.p,  &ws.dringAtoms.p, &ws.dhc.p, &ws.dddc.p,
                       &nA,            &nFr,             &mR,       &ws.datomHc.p,
                       &ws.datomDdc.p, &ws.dsix.p};
        auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
        kAtom->launch(dim3(static_cast<unsigned>(ringGrid)), dim3(static_cast<unsigned>(block)), 0, q, a, false);
      }
      knn.wait();
      out.computeMs += msSince(t0);

      t0 = std::chrono::steady_clock::now();
      auto down = [&](void *dst, const void *src, std::size_t bytes,
                      const char *what) {
        checkCuda(rt.cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost), what);
      };
      down(nRings.data(), ws.dnRings.p, nf * sizeof(int), "DtoH nRings");
      down(nPairs.data(), ws.dnPairs.p, nf * sizeof(int), "DtoH nPairs");
      down(rowMax.data(), ws.drowMax.p, nf * sizeof(int), "DtoH rowMax");
      down(six.data(), ws.dsix.p, nN * sizeof(int), "DtoH six");
      down(hcAtom.data(), ws.datomHc.p, nN * sizeof(int), "DtoH hc");
      down(ddcAtom.data(), ws.datomDdc.p, nN * sizeof(int), "DtoH ddc");
      int needPer = maxPer;
      int needRings = maxRings;
      int needPairs = maxPairs;
      bool pending = false;
      for (std::size_t f = 0; f < nf; ++f) {
        if (done[f] || onHost[f]) {
          continue;
        }
        if (rowMax[f] > kRowCap) {
          onHost[f] = 1;
        } else if (nRings[f] <= maxRings && rowMax[f] <= maxPer &&
                   nPairs[f] <= maxPairs) {
          const auto at = static_cast<std::ptrdiff_t>(f * n);
          const auto len = static_cast<std::ptrdiff_t>(n);
          std::copy(six.begin() + at, six.begin() + at + len,
                    out.sixCount.begin() + at);
          std::copy(hcAtom.begin() + at, hcAtom.begin() + at + len,
                    out.atomHc.begin() + at);
          std::copy(ddcAtom.begin() + at, ddcAtom.begin() + at + len,
                    out.atomDdc.begin() + at);
          out.nRings[f] = nRings[f];
          done[f] = 1;
          continue;
        } else {
          needPer = std::max(needPer, rowMax[f]);
          needRings = std::max(needRings, nRings[f]);
          needPairs = std::max(needPairs, nPairs[f]);
          pending = true;
        }
        hostDeg[f].resize(n);
        hostCols[f].resize(n * 4);
        down(hostDeg[f].data(), static_cast<int *>(ws.dmdeg.p) + f * n,
             n * sizeof(int), "DtoH mdeg");
        down(hostCols[f].data(), static_cast<int *>(ws.dmcols.p) + f * n * 4,
             n * 4 * sizeof(int), "DtoH mcols");
      }
      out.downloadMs += msSince(t0);
      if (!pending) {
        break;
      }
      maxPer = needPer;
      maxRings = needRings;
      maxPairs = needPairs;
      ++out.retries;
    }

    // Host labels of each frame the device could not hold, on the mutual
    // graph the device built for it
    const auto t0 = std::chrono::steady_clock::now();
    for (std::size_t f = 0; f < nf; ++f) {
      if (!onHost[f]) {
        continue;
      }
      std::vector<std::vector<int>> nList(n);
      for (std::size_t i = 0; i < n; ++i) {
        nList[i].push_back(static_cast<int>(i));
        for (int t = 0; t < hostDeg[f][i]; ++t) {
          nList[i].push_back(hostCols[f][i * 4 + static_cast<std::size_t>(t)]);
        }
      }
      const auto host = tum::hostCageCounts(nList);
      const auto at = static_cast<std::ptrdiff_t>(f * n);
      std::copy(host.atomSix.begin(), host.atomSix.end(),
                out.sixCount.begin() + at);
      std::copy(host.atomHc.begin(), host.atomHc.end(), out.atomHc.begin() + at);
      std::copy(host.atomDdc.begin(), host.atomDdc.end(),
                out.atomDdc.begin() + at);
      out.nRings[f] = host.nSix;
      out.hostFrames.push_back(static_cast<int>(f));
    }
    out.hostMs = msSince(t0);
  } catch (const std::exception &ex) {
    out.error = ex.what();
  }
  return out;
#endif
}

} // namespace gpu
