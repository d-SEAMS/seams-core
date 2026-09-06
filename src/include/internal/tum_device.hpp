#ifndef SEAMS_TUM_DEVICE_H_
#define SEAMS_TUM_DEVICE_H_

/** @file tum_device.hpp
 *  @brief Device-safe TUM ice-score pieces: hop-bound primitive six-rings
 *  and claim-free HC/DDC affiliation.
 *
 *  No STL containers, no function-local statics. Host OpenMP and OpenMP
 *  target regions call the same functions. CHILL and q_lm are not here.
 *
 *  A six-cycle is primitive when the hop-bounded Franzblau test holds:
 *  every pair of non-adjacent members is at least as far through the
 *  graph as around the ring. On a 4-regular bond graph that is no
 *  chords (ring hops 2) and no common neighbour of opposite vertices
 *  (ring hops 3).
 */

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp declare target
#endif

namespace tum {
namespace device {

inline bool bonded(const int *deg, const int *cols, int nAtoms, int kMax,
                   int a, int b) {
  if (a < 0 || b < 0 || a >= nAtoms || b >= nAtoms) {
    return false;
  }
  const int d = deg[a];
  const int row = a * kMax;
  for (int t = 0; t < d; ++t) {
    if (cols[row + t] == b) {
      return true;
    }
  }
  return false;
}

inline bool inSix(const int *r, int atom) {
  for (int t = 0; t < 6; ++t) {
    if (r[t] == atom) {
      return true;
    }
  }
  return false;
}

inline bool shareAtoms(const int *a, const int *b) {
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i])) {
      return true;
    }
  }
  return false;
}

inline int commonCount(const int *a, const int *b) {
  int n = 0;
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i])) {
      ++n;
    }
  }
  return n;
}

inline bool commonInThree(const int *a, const int *b, const int *c) {
  for (int i = 0; i < 6; ++i) {
    if (inSix(b, a[i]) && inSix(c, a[i])) {
      return true;
    }
  }
  return false;
}

inline bool shareNeigh(const int *deg, const int *cols, int nAtoms, int kMax,
                       int a, int b) {
  const int da = deg[a];
  const int ra = a * kMax;
  for (int t = 0; t < da; ++t) {
    if (bonded(deg, cols, nAtoms, kMax, b, cols[ra + t])) {
      return true;
    }
  }
  return false;
}

/** Graph distance from a to b when it is at most cap, else -1.
 *  cap is 1 or 2: the Franzblau reject bound for a six-ring. */
inline int hopsAtMost(const int *deg, const int *cols, int nAtoms, int kMax,
                      int a, int b, int cap) {
  if (a == b) {
    return 0;
  }
  if (bonded(deg, cols, nAtoms, kMax, a, b)) {
    return 1;
  }
  if (cap < 2) {
    return -1;
  }
  if (shareNeigh(deg, cols, nAtoms, kMax, a, b)) {
    return 2;
  }
  return -1;
}

/** Hop-bounded Franzblau SP test on a six-cycle of the bond graph. */
inline bool hopBoundPrimitiveSix(const int *r, const int *deg, const int *cols,
                                 int nAtoms, int kMax) {
  const int pairI[9] = {0, 1, 2, 3, 4, 5, 0, 1, 2};
  const int pairJ[9] = {2, 3, 4, 5, 0, 1, 3, 4, 5};
  const int ringHops[9] = {2, 2, 2, 2, 2, 2, 3, 3, 3};
  for (int t = 0; t < 9; ++t) {
    if (hopsAtMost(deg, cols, nAtoms, kMax, r[pairI[t]], r[pairJ[t]],
                   ringHops[t] - 1) >= 0) {
      return false;
    }
  }
  return true;
}

inline bool basalNeighbours(const int *deg, const int *cols, int nAtoms,
                            int kMax, int n1, int n2, int atomOne,
                            int atomTwo) {
  const bool n1one = bonded(deg, cols, nAtoms, kMax, atomOne, n1);
  const bool n1two = bonded(deg, cols, nAtoms, kMax, atomTwo, n1);
  if (!n1one && !n1two) {
    return false;
  }
  if (n1one) {
    return bonded(deg, cols, nAtoms, kMax, atomTwo, n2);
  }
  return bonded(deg, cols, nAtoms, kMax, atomOne, n2);
}

inline bool notNeighboursOfRing(const int *deg, const int *cols, int nAtoms,
                                int kMax, const int *trip, const int *ring) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      if (bonded(deg, cols, nAtoms, kMax, ring[j], trip[i])) {
        return false;
      }
    }
  }
  return true;
}

inline bool basalConditions(const int *deg, const int *cols, int nAtoms,
                            int kMax, const int *b1, const int *b2) {
  int kIndex = -1;
  int compare1 = 0;
  int compare2 = 0;
  bool l1n = false;
  bool l2n = false;
  const int l1 = b1[0];
  const int l2 = b1[1];
  for (int k = 0; k < 6; ++k) {
    const int mk = b2[k];
    if (bonded(deg, cols, nAtoms, kMax, l1, mk)) {
      compare1 = b1[2];
      compare2 = b1[4];
      kIndex = k;
      l1n = true;
      break;
    }
    if (bonded(deg, cols, nAtoms, kMax, l2, mk)) {
      compare1 = b1[3];
      compare2 = b1[5];
      kIndex = k;
      l2n = true;
      break;
    }
  }
  if (!l1n && !l2n) {
    return false;
  }
  int evenT[3];
  int oddT[3];
  int ie = 0;
  int io = 0;
  for (int k = 0; k <= 5; ++k) {
    int ck = kIndex + k;
    if (ck >= 6) {
      ck -= 6;
    }
    if (k % 2 == 0) {
      evenT[ie++] = b2[ck];
    } else {
      oddT[io++] = b2[ck];
    }
  }
  if (!basalNeighbours(deg, cols, nAtoms, kMax, evenT[1], evenT[2], compare1,
                       compare2)) {
    return false;
  }
  return notNeighboursOfRing(deg, cols, nAtoms, kMax, oddT, b1);
}

inline int firstRingThrough(const int *A, int nA, const int *B, int nB,
                            const int *C, int nC, int skipA, int skipB) {
  int i = 0;
  int j = 0;
  int k = 0;
  while (i < nA && j < nB && k < nC) {
    const int x = A[i];
    const int y = B[j];
    const int z = C[k];
    if (x == y && y == z) {
      if (x != skipA && x != skipB) {
        return x;
      }
      ++i;
      ++j;
      ++k;
      continue;
    }
    int lo = x;
    if (y < lo) {
      lo = y;
    }
    if (z < lo) {
      lo = z;
    }
    if (x == lo) {
      ++i;
    }
    if (y == lo) {
      ++j;
    }
    if (z == lo) {
      ++k;
    }
  }
  return -1;
}

inline int ringsThrough(const int *A, int nA, const int *B, int nB,
                        const int *C, int nC, int skipA, int skipB, int *out,
                        int cap) {
  int i = 0;
  int j = 0;
  int k = 0;
  int n = 0;
  while (i < nA && j < nB && k < nC) {
    const int x = A[i];
    const int y = B[j];
    const int z = C[k];
    if (x == y && y == z) {
      if (x != skipA && x != skipB && n < cap) {
        out[n++] = x;
      }
      ++i;
      ++j;
      ++k;
      continue;
    }
    int lo = x;
    if (y < lo) {
      lo = y;
    }
    if (z < lo) {
      lo = z;
    }
    if (x == lo) {
      ++i;
    }
    if (y == lo) {
      ++j;
    }
    if (z == lo) {
      ++k;
    }
  }
  return n;
}

inline int fetchAdd(int *p) {
  int old;
#if defined(_OPENMP)
#pragma omp atomic capture
  old = (*p)++;
#else
  old = (*p)++;
#endif
  return old;
}

/** Enumerate hop-bound primitive six-rings whose lowest vertex is i. */
inline void enumSixFrom(int i, const int *deg, const int *cols, int nAtoms,
                        int kMax, int maxRings, int *nRings, int *ringAtoms,
                        int *dropped) {
  const int di = deg[i];
  if (di < 2) {
    return;
  }
  const int irow = i * kMax;
  for (int ia = 0; ia < di; ++ia) {
    const int a = cols[irow + ia];
    for (int ib = ia + 1; ib < di; ++ib) {
      const int b = cols[irow + ib];
      const int da = deg[a];
      const int db = deg[b];
      const int arow = a * kMax;
      const int brow = b * kMax;
      for (int ix = 0; ix < da; ++ix) {
        const int x = cols[arow + ix];
        if (x == i || x == b) {
          continue;
        }
        for (int iy = 0; iy < db; ++iy) {
          const int y = cols[brow + iy];
          if (y == i || y == a || y == x) {
            continue;
          }
          const int dx = deg[x];
          const int xrow = x * kMax;
          for (int iz = 0; iz < dx; ++iz) {
            const int z = cols[xrow + iz];
            if (z == i || z == a || z == b || z == y) {
              continue;
            }
            if (!bonded(deg, cols, nAtoms, kMax, z, y)) {
              continue;
            }
            const int cyc[6] = {i, a, x, z, y, b};
            int mn = i;
            for (int t = 1; t < 6; ++t) {
              if (cyc[t] < mn) {
                mn = cyc[t];
              }
            }
            if (mn != i) {
              continue;
            }
            if (!hopBoundPrimitiveSix(cyc, deg, cols, nAtoms, kMax)) {
              continue;
            }
            const int slot = fetchAdd(nRings);
            if (slot >= maxRings) {
              fetchAdd(dropped);
              continue;
            }
            const int dest = slot * 6;
            for (int t = 0; t < 6; ++t) {
              ringAtoms[dest + t] = cyc[t];
            }
          }
        }
      }
    }
  }
}

inline void invertOneRing(int r, const int *nRings, const int *ringAtoms,
                          int nAtoms, int maxPer, int *throughCount,
                          int *through) {
  if (r >= *nRings) {
    return;
  }
  const int *ring = ringAtoms + r * 6;
  for (int t = 0; t < 6; ++t) {
    const int atom = ring[t];
    if (atom < 0 || atom >= nAtoms) {
      continue;
    }
    const int slot = fetchAdd(throughCount + atom);
    if (slot < maxPer) {
      through[atom * maxPer + slot] = r;
    }
  }
}

inline void sortThroughRow(int *row, int n) {
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

inline void emitBasalFrom(int i, const int *nRings, const int *ringAtoms,
                          const int *deg, const int *cols,
                          const int *throughCount, const int *through,
                          int nAtoms, int kMax, int maxPer, int maxPairs,
                          int *nPairs, int *pairs) {
  if (i >= *nRings) {
    return;
  }
  const int *bi = ringAtoms + i * 6;
  for (int slot = 0; slot < 2; ++slot) {
    const int anchor = bi[slot];
    if (anchor < 0 || anchor >= nAtoms) {
      continue;
    }
    const int da = deg[anchor];
    const int arow = anchor * kMax;
    for (int n = 0; n < da; ++n) {
      const int nb = cols[arow + n];
      if (nb < 0 || nb >= nAtoms) {
        continue;
      }
      const int nr = throughCount[nb];
      const int *row = through + nb * maxPer;
      for (int t = 0; t < nr; ++t) {
        const int j = row[t];
        if (j == i) {
          continue;
        }
        const int *bj = ringAtoms + j * 6;
        if (shareAtoms(bi, bj)) {
          continue;
        }
        if (!basalConditions(deg, cols, nAtoms, kMax, bi, bj)) {
          continue;
        }
        const int slotp = fetchAdd(nPairs);
        if (slotp < maxPairs) {
          pairs[slotp * 2] = i;
          pairs[slotp * 2 + 1] = j;
        }
      }
    }
  }
}

inline void applyHcPair(int p, const int *nPairs, const int *pairs,
                        const int *ringAtoms, const int *throughCount,
                        const int *through, int nAtoms, int maxPer, int *hc) {
  if (p >= *nPairs) {
    return;
  }
  const int i = pairs[p * 2];
  const int j = pairs[p * 2 + 1];
  hc[i] = 1;
  hc[j] = 1;
  const int *bi = ringAtoms + i * 6;
  const int *bj = ringAtoms + j * 6;
  for (int q = 0; q < 6; ++q) {
    int trip[3];
    for (int m = 0; m < 3; ++m) {
      trip[m] = bi[(q + m) % 6];
    }
    const int nA = throughCount[trip[0]];
    const int nB = throughCount[trip[1]];
    const int nC = throughCount[trip[2]];
    const int *A = through + trip[0] * maxPer;
    const int *B = through + trip[1] * maxPer;
    const int *C = through + trip[2] * maxPer;
    int cand[16];
    const int nc = ringsThrough(A, nA, B, nB, C, nC, i, j, cand, 16);
    for (int c = 0; c < nc; ++c) {
      const int kr = cand[c];
      const int *bk = ringAtoms + kr * 6;
      int rest[3];
      int nrst = 0;
      for (int u = 0; u < 6; ++u) {
        if (!inSix(trip, bk[u]) && nrst < 3) {
          rest[nrst++] = bk[u];
        }
      }
      if (nrst == 3 && commonCount(rest, bj) == 3) {
        hc[kr] = 1;
      }
    }
  }
}

inline void ddcFrom(int i, const int *nRings, const int *ringAtoms,
                    const int *throughCount, const int *through, const int *hc,
                    int nAtoms, int maxPer, int *ddc) {
  if (i >= *nRings) {
    return;
  }
  if (hc[i]) {
    return;
  }
  const int *bi = ringAtoms + i * 6;
  int peri[32];
  int nPeri = 0;
  for (int m = 0; m < 6; ++m) {
    const int atom = bi[m];
    if (atom < 0 || atom >= nAtoms) {
      return;
    }
    const int nr = throughCount[atom];
    const int *row = through + atom * maxPer;
    int common = 0;
    for (int t = 0; t < nr; ++t) {
      if (row[t] == i) {
        continue;
      }
      ++common;
      if (nPeri < 32) {
        peri[nPeri++] = row[t];
      }
    }
    if (common < 3) {
      return;
    }
  }
  int newP[6];
  for (int k = 0; k < 6; ++k) {
    int trip[3];
    for (int t = 0; t < 3; ++t) {
      trip[t] = bi[(k + t) % 6];
    }
    const int nA = throughCount[trip[0]];
    const int nB = throughCount[trip[1]];
    const int nC = throughCount[trip[2]];
    const int *A = through + trip[0] * maxPer;
    const int *B = through + trip[1] * maxPer;
    const int *C = through + trip[2] * maxPer;
    const int j = firstRingThrough(A, nA, B, nB, C, nC, i, -1);
    if (j < 0) {
      return;
    }
    newP[k] = j;
  }
  const int *p0 = ringAtoms + newP[0] * 6;
  const int *p1 = ringAtoms + newP[1] * 6;
  const int *p2 = ringAtoms + newP[2] * 6;
  const int *p3 = ringAtoms + newP[3] * 6;
  const int *p4 = ringAtoms + newP[4] * 6;
  const int *p5 = ringAtoms + newP[5] * 6;
  if (!commonInThree(p0, p2, p4)) {
    return;
  }
  if (!commonInThree(p1, p3, p5)) {
    return;
  }
  const int *pairs[4][2] = {{p0, p2}, {p1, p3}, {p2, p4}, {p3, p5}};
  for (int t = 0; t < 4; ++t) {
    if (commonCount(pairs[t][0], pairs[t][1]) < 3) {
      return;
    }
  }
  ddc[i] = 1;
  for (int t = 0; t < 6; ++t) {
    ddc[newP[t]] = 1;
  }
}

inline void atomIceFrom(int r, const int *nRings, const int *ringAtoms,
                        const int *hc, const int *ddc, int nAtoms, int *atomHc,
                        int *atomDdc) {
  if (r >= *nRings) {
    return;
  }
  const int *ring = ringAtoms + r * 6;
  const int isHc = hc[r];
  const int isDdc = ddc[r];
  for (int t = 0; t < 6; ++t) {
    const int a = ring[t];
    if (a < 0 || a >= nAtoms) {
      continue;
    }
    if (isHc) {
      atomHc[a] = 1;
    }
    if (isDdc) {
      atomDdc[a] = 1;
    }
  }
}

} // namespace device
} // namespace tum

#ifdef SEAMS_HAS_OFFLOAD
#pragma omp end declare target
#endif

#endif
