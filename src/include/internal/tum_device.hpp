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

/** Longest atom row the kernels accept. Candidate buffers are sized by
 *  it and ddcFrom keeps one bit per entry in an unsigned. A liquid atom
 *  of a 3.5 A cutoff graph can lie in more than 16 six-rings; a graph
 *  with a longer row goes to the host (tum_offload.hpp). */
constexpr int kRowCap = 32;
static_assert(kRowCap <= 32, "ddcFrom marks row entries in 32-bit masks");

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

/** beta(b1, b2) for rings sharing no atom: some alternating triplet
 *  {l_p, l_p+2, l_p+4} of b1 bonds one-to-one onto an alternating triplet
 *  {m_q, m_q+2, m_q+4} of b2, and no other atom of b2 is bonded to b1.
 *  All 2 x 2 x 6 alignments are tried, so rotation, direction and ring
 *  order do not matter. */
inline bool basalConditions(const int *deg, const int *cols, int nAtoms,
                            int kMax, const int *b1, const int *b2) {
  const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                          {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  for (int q = 0; q < 2; ++q) {
    const int oddT[3] = {b2[q + 1], b2[q + 3], b2[(q + 5) % 6]};
    if (!notNeighboursOfRing(deg, cols, nAtoms, kMax, oddT, b1)) {
      continue;
    }
    for (int p = 0; p < 2; ++p) {
      for (int s = 0; s < 6; ++s) {
        bool all = true;
        for (int t = 0; t < 3 && all; ++t) {
          all = bonded(deg, cols, nAtoms, kMax, b1[p + 2 * t],
                       b2[q + 2 * perm[s][t]]);
        }
        if (all) {
          return true;
        }
      }
    }
  }
  return false;
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

/** Enumerate hop-bound primitive six-rings whose lowest vertex is i.
 *  nRings counts every ring found; only the first maxRings are stored. */
inline void enumSixFrom(int i, const int *deg, const int *cols, int nAtoms,
                        int kMax, int maxRings, int *nRings, int *ringAtoms) {
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

/** throughCount counts every stored ring through the atom; only the
 *  first maxPer ids are kept. */
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

/** True when ring j is in the through row of a bond neighbour of
 *  bi[0] or bi[1] visited before (slot, n), so each candidate is tested
 *  once. */
inline bool seenEarlier(const int *bi, const int *deg, const int *cols,
                        const int *throughCount, const int *through,
                        int nAtoms, int kMax, int maxPer, int slot, int n,
                        int j) {
  for (int s = 0; s <= slot; ++s) {
    const int anchor = bi[s];
    if (anchor < 0 || anchor >= nAtoms) {
      continue;
    }
    const int arow = anchor * kMax;
    const int stop = s < slot ? deg[anchor] : n;
    for (int m = 0; m < stop; ++m) {
      const int nb = cols[arow + m];
      if (nb < 0 || nb >= nAtoms) {
        continue;
      }
      const int nr = throughCount[nb];
      const int *row = through + nb * maxPer;
      for (int t = 0; t < nr; ++t) {
        if (row[t] == j) {
          return true;
        }
      }
    }
  }
  return false;
}

/** Ordered basal pairs (i, j). Every pair passing beta has a bond from
 *  bi[0] or bi[1] into bj, so the neighbours of those two atoms reach
 *  every candidate. */
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
        if (seenEarlier(bi, deg, cols, throughCount, through, nAtoms, kMax,
                        maxPer, slot, n, j)) {
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
    int cand[kRowCap];
    const int nc = ringsThrough(A, nA, B, nB, C, nC, i, j, cand, kRowCap);
    for (int c = 0; c < nc; ++c) {
      const int kr = cand[c];
      const int *bk = ringAtoms + kr * 6;
      int rest = 0;
      int restInJ = 0;
      for (int u = 0; u < 6; ++u) {
        const int a = bk[u];
        if (a != trip[0] && a != trip[1] && a != trip[2]) {
          ++rest;
          if (inSix(bj, a)) {
            ++restInJ;
          }
        }
      }
      if (rest == 3 && restInJ == 3) {
        hc[kr] = 1;
      }
    }
  }
}

/** epsilon(i) for a ring that is not HC, as ring::equatorialConditions:
 *  every atom in 3 other rings, and pairwise distinct P_k in T_k with P_0,
 *  P_2, P_4 sharing an atom, P_1, P_3, P_5 sharing an atom and, for some
 *  starting atom and direction of i, the alternate pairs sharing 3 atoms.
 *  On a pass marks i and every P_k of every valid assignment. */
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
  for (int m = 0; m < 6; ++m) {
    const int atom = bi[m];
    if (atom < 0 || atom >= nAtoms) {
      return;
    }
    const int nr = throughCount[atom];
    const int *row = through + atom * maxPer;
    int common = 0;
    for (int t = 0; t < nr; ++t) {
      if (row[t] != i) {
        ++common;
      }
    }
    if (common < 3) {
      return;
    }
  }
  // T[k]: every ring but i through bi[k], bi[k+1], bi[k+2]
  int T[6][kRowCap];
  int nT[6];
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
    nT[k] = ringsThrough(A, nA, B, nB, C, nC, i, -1, T[k], kRowCap);
    if (nT[k] == 0) {
      return;
    }
  }
  // Bit t of used[k] marks T[k][t] as P_k of some valid assignment
  unsigned used[6] = {0U, 0U, 0U, 0U, 0U, 0U};
  bool any = false;
  for (int a = 0; a < nT[0]; ++a) {
    const int r0 = T[0][a];
    const int *p0 = ringAtoms + r0 * 6;
    for (int c = 0; c < nT[2]; ++c) {
      const int r2 = T[2][c];
      const int *p2 = ringAtoms + r2 * 6;
      if (r2 == r0) {
        continue;
      }
      const bool e02 = commonCount(p0, p2) >= 3;
      for (int e = 0; e < nT[4]; ++e) {
        const int r4 = T[4][e];
        const int *p4 = ringAtoms + r4 * 6;
        if (r4 == r0 || r4 == r2 || !commonInThree(p0, p2, p4)) {
          continue;
        }
        const bool e24 = commonCount(p2, p4) >= 3;
        const bool e40 = commonCount(p4, p0) >= 3;
        // Bit t: the ring on triplet 2t shares 3 atoms with both others
        const unsigned mE = (e40 && e02 ? 1U : 0U) |
                            (e02 && e24 ? 2U : 0U) |
                            (e24 && e40 ? 4U : 0U);
        if (mE == 0U) {
          continue;
        }
        // Some start and direction of i puts an even middle 2t next to an
        // odd middle 2t+1 or 2t-1
        unsigned partners = 0U;
        for (int t = 0; t < 3; ++t) {
          if ((mE >> t) & 1U) {
            partners |= (1U << t) | (1U << ((t + 2) % 3));
          }
        }
        for (int b = 0; b < nT[1]; ++b) {
          const int r1 = T[1][b];
          const int *p1 = ringAtoms + r1 * 6;
          if (r1 == r0 || r1 == r2 || r1 == r4) {
            continue;
          }
          for (int d = 0; d < nT[3]; ++d) {
            const int r3 = T[3][d];
            const int *p3 = ringAtoms + r3 * 6;
            if (r3 == r0 || r3 == r2 || r3 == r4 || r3 == r1) {
              continue;
            }
            const bool o13 = commonCount(p1, p3) >= 3;
            for (int g = 0; g < nT[5]; ++g) {
              const int r5 = T[5][g];
              const int *p5 = ringAtoms + r5 * 6;
              if (r5 == r0 || r5 == r2 || r5 == r4 || r5 == r1 || r5 == r3 ||
                  !commonInThree(p1, p3, p5)) {
                continue;
              }
              const bool o35 = commonCount(p3, p5) >= 3;
              const bool o51 = commonCount(p5, p1) >= 3;
              const unsigned mO = (o51 && o13 ? 1U : 0U) |
                                  (o13 && o35 ? 2U : 0U) |
                                  (o35 && o51 ? 4U : 0U);
              if ((partners & mO) == 0U) {
                continue;
              }
              any = true;
              used[0] |= 1U << a;
              used[1] |= 1U << b;
              used[2] |= 1U << c;
              used[3] |= 1U << d;
              used[4] |= 1U << e;
              used[5] |= 1U << g;
            }
          }
        }
      }
    }
  }
  if (!any) {
    return;
  }
  ddc[i] = 1;
  for (int k = 0; k < 6; ++k) {
    for (int t = 0; t < nT[k]; ++t) {
      if ((used[k] >> t) & 1U) {
        ddc[T[k][t]] = 1;
      }
    }
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
