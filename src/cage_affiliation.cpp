//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <cage_affiliation.hpp>

#include <ring.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

//! Fingerprint of a ring's cycle, the same however the ring was rotated or
//! traversed: ring atoms are distinct, so the lexicographically least
//! reading starts at the least atom and runs toward its lesser neighbour
std::uint64_t cycleHash(const int *atoms, std::size_t k) {
  std::uint64_t h = 1469598103934665603ULL ^ k;
  if (k == 0) {
    return h;
  }
  const std::size_t start =
      static_cast<std::size_t>(std::min_element(atoms, atoms + k) - atoms);
  const bool forward = atoms[(start + 1) % k] <= atoms[(start + k - 1) % k];
  for (std::size_t m = 0; m < k; m++) {
    h ^= static_cast<std::uint32_t>(
        atoms[forward ? (start + m) % k : (start + k - m) % k]);
    h *= 1099511628211ULL;
  }
  return h;
}

bool sameCycle(const int *a, std::size_t ka, const int *b, std::size_t kb) {
  if (ka != kb) {
    return false;
  }
  if (ka == 0) {
    return true;
  }
  const std::size_t sa =
      static_cast<std::size_t>(std::min_element(a, a + ka) - a);
  const std::size_t sb =
      static_cast<std::size_t>(std::min_element(b, b + kb) - b);
  const bool fa = a[(sa + 1) % ka] <= a[(sa + ka - 1) % ka];
  const bool fb = b[(sb + 1) % kb] <= b[(sb + kb - 1) % kb];
  for (std::size_t m = 0; m < ka; m++) {
    if (a[fa ? (sa + m) % ka : (sa + ka - m) % ka] !=
        b[fb ? (sb + m) % kb : (sb + kb - m) % kb]) {
      return false;
    }
  }
  return true;
}

//! Shared per-frame machinery: the inverted index, sorted ring copies for
//! overlap tests, and cached basal/adjacency lists (one sort per ring).
struct FrameContext {
  // Pointers so a caller can keep one context and rebind it each frame.
  const std::vector<std::vector<int>> *rings = nullptr;
  const std::vector<std::vector<int>> *nList = nullptr;
  ring::RingSearchIndex index;
  std::vector<std::vector<int>> sortedRings;
  mutable   std::vector<std::vector<int>> basalCand;
  mutable std::vector<char> basalReady;
  mutable std::vector<std::vector<int>> adjRings;
  mutable std::vector<char> adjReady;
  ring::NeighBond bonds;
  // Scratch for findPrismatic, which appends to a list and stamps a type
  // vector neither of which affiliation reads
  std::vector<int> scratchListHC;
  std::vector<ring::strucType> scratchType;
  // Ring atoms in one block, six per ring, when every ring has length 6.
  // The basal search reads these atoms for every candidate.
  std::vector<int> flat6;
  bool allHex = false;

  FrameContext() = default;

  FrameContext(const std::vector<std::vector<int>> &rings_,
               const std::vector<std::vector<int>> &nList_) {
    prepare(rings_, nList_);
  }

  void prepare(const std::vector<std::vector<int>> &rings_,
               const std::vector<std::vector<int>> &nList_) {
    rings = &rings_;
    nList = &nList_;
    int maxAtom = 0;
    for (const auto &r : rings_) {
      for (const int a : r) {
        maxAtom = std::max(maxAtom, a);
      }
    }
    index = ring::buildRingSearchIndex(rings_, maxAtom + 1);
    bonds.build(nList_);
    const int nR = static_cast<int>(rings_.size());
    sortedRings.resize(static_cast<std::size_t>(nR));
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel for schedule(static) if (nR >= 1024 && !omp_in_parallel())
#endif
    for (int i = 0; i < nR; i++) {
      sortedRings[static_cast<std::size_t>(i)] =
          rings_[static_cast<std::size_t>(i)];
      std::sort(sortedRings[static_cast<std::size_t>(i)].begin(),
                sortedRings[static_cast<std::size_t>(i)].end());
    }
    scratchType.assign(static_cast<std::size_t>(nR),
                       ring::strucType::unclassified);
    allHex = true;
    for (int i = 0; i < nR; i++) {
      if (rings_[static_cast<std::size_t>(i)].size() != 6) {
        allHex = false;
        break;
      }
    }
    flat6.clear();
    if (allHex) {
      flat6.resize(static_cast<std::size_t>(nR) * 6);
      for (int i = 0; i < nR; i++) {
        const auto &ring = rings_[static_cast<std::size_t>(i)];
        std::copy(ring.begin(), ring.end(),
                  flat6.begin() + static_cast<std::ptrdiff_t>(i) * 6);
      }
    }
    basalCand.resize(static_cast<std::size_t>(nR));
    adjRings.resize(static_cast<std::size_t>(nR));
    adjReady.assign(static_cast<std::size_t>(nR), 0);
    // Left clear. basalCandidates fills a ring the first time a sweep asks.
    basalReady.assign(static_cast<std::size_t>(nR), 0);
  }

  [[nodiscard]] bool shareAtoms(int i, int j) const {
    const auto &a = sortedRings[i];
    const auto &b = sortedRings[j];
    size_t p = 0, q = 0;
    while (p < a.size() && q < b.size()) {
      if (a[p] == b[q]) {
        return true;
      }
      (a[p] < b[q]) ? p++ : q++;
    }
    return false;
  }

  struct RingRow {
    const int *p = nullptr;
    int n = 0;
  };

  [[nodiscard]] RingRow ringsOf(int atom) const {
    if (atom < 0 || index.ringOff.empty() ||
        static_cast<std::size_t>(atom) + 1 >= index.ringOff.size()) {
      return {};
    }
    const int begin = index.ringOff[static_cast<std::size_t>(atom)];
    const int end = index.ringOff[static_cast<std::size_t>(atom) + 1];
    return {index.ringCol.data() + begin, end - begin};
  }

  [[nodiscard]] const std::vector<int> &ringsThroughAtom(int atom) const {
    static const std::vector<int> empty;
    if (atom < 0 ||
        atom >= static_cast<int>(index.ringsContainingAtom.size())) {
      return empty;
    }
    return index.ringsContainingAtom[atom];
  }

  [[nodiscard]] std::vector<int> collectAdjacentRings(int r) const {
    std::vector<int> out;
    for (const int atom : (*rings)[r]) {
      for (const int s : ringsThroughAtom(atom)) {
        if (s != r) {
          out.push_back(s);
        }
      }
      if (atom >= 0 && atom < static_cast<int>(nList->size())) {
        for (size_t n = 1; n < (*nList)[atom].size(); n++) {
          for (const int s : ringsThroughAtom((*nList)[atom][n])) {
            if (s != r) {
              out.push_back(s);
            }
          }
        }
      }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  //! Ordered basal-pair predicate, exactly the findHC acceptance test
  [[nodiscard]] bool basalPair(int i, int j) const {
    if (shareAtoms(i, j)) {
      return false;
    }
    return ring::basalConditions(bonds, (*rings)[static_cast<std::size_t>(i)],
                                 (*rings)[static_cast<std::size_t>(j)]);
  }

  [[nodiscard]] const std::vector<int> &adjacentRings(int r) {
    static const std::vector<int> empty;
    if (r < 0 || r >= static_cast<int>(adjRings.size())) {
      return empty;
    }
    if (!adjReady[static_cast<std::size_t>(r)]) {
      adjRings[static_cast<std::size_t>(r)] = collectAdjacentRings(r);
      adjReady[static_cast<std::size_t>(r)] = 1;
    }
    return adjRings[static_cast<std::size_t>(r)];
  }

  //! Rings that can pass basalConditions with ring i as the first argument.
  //! A passing partner holds a neighbour of each atom of some alternate
  //! triplet of i. Other rings are not candidates.
  struct BasalScratch {
    std::vector<int> onRing;
    std::vector<int> arm1;
    std::vector<int> arm2;
    std::vector<int> seen;
    // Unsigned, so the generation counter wraps by the language rules.
    // A signed increment that passed INT_MAX would be undefined, and the
    // reset under that test was not a path the compiler could take.
    std::uint8_t tickOn = 0;
    std::uint8_t tick1 = 0;
    std::uint8_t tick2 = 0;
    std::uint8_t tickSeen = 0;
    void prepare(int atoms, int nRings) {
      if (static_cast<int>(onRing.size()) < atoms) {
        onRing.assign(static_cast<size_t>(atoms), 0);
        arm1.assign(static_cast<size_t>(atoms), 0);
        arm2.assign(static_cast<size_t>(atoms), 0);
      }
      if (static_cast<int>(seen.size()) < nRings) {
        seen.assign(static_cast<size_t>(nRings), 0);
      }
    }
    int next(std::uint8_t &tick, std::vector<int> &stamp) {
      ++tick;
      if (tick == 0) {
        std::fill(stamp.begin(), stamp.end(), 0);
        tick = 1;
      }
      return tick;
    }
  };

  [[nodiscard]] std::vector<int> collectBasalCandidates(int i,
                                                        BasalScratch &scratch) const {
    std::vector<int> out;
    const std::vector<int> &basal1 = (*rings)[static_cast<size_t>(i)];
    if (basal1.size() != 6) {
      return out;
    }
    const int nAtoms = static_cast<int>(nList->size());
    const int nRings = static_cast<int>(rings->size());
    scratch.prepare(nAtoms, nRings);
    const int on = scratch.next(scratch.tickOn, scratch.onRing);
    for (const int atom : basal1) {
      if (atom >= 0 && atom < nAtoms) {
        scratch.onRing[static_cast<size_t>(atom)] = on;
      }
    }
    auto neighbourSpan = [&](int atom) {
      int span = 0;
      if (atom < 0 || atom >= nAtoms) {
        return span;
      }
      const auto &row = (*nList)[static_cast<size_t>(atom)];
      for (size_t n = 1; n < row.size(); n++) {
        const int nb = row[n];
        if (nb < 0 || nb >= nAtoms ||
            scratch.onRing[static_cast<size_t>(nb)] == on) {
          continue;
        }
        span += ringsOf(nb).n;
      }
      return span;
    };
    auto stampArm = [&](int atom, std::vector<int> &stamp, int tick) {
      if (atom < 0 || atom >= nAtoms) {
        return;
      }
      const auto &row = (*nList)[static_cast<size_t>(atom)];
      for (size_t n = 1; n < row.size(); n++) {
        const int nb = row[n];
        if (nb >= 0 && nb < nAtoms &&
            scratch.onRing[static_cast<size_t>(nb)] != on) {
          stamp[static_cast<size_t>(nb)] = tick;
        }
      }
    };
    for (int parity = 0; parity < 2; parity++) {
      int arm = 0;
      int best = 0;
      for (int t = 0; t < 3; t++) {
        const int span =
            neighbourSpan(basal1[static_cast<size_t>((parity + 2 * t) % 6)]);
        if (t == 0 || span < best) {
          best = span;
          arm = t;
        }
      }
      const int tick1 = scratch.next(scratch.tick1, scratch.arm1);
      const int tick2 = scratch.next(scratch.tick2, scratch.arm2);
      stampArm(basal1[static_cast<size_t>((parity + 2 * ((arm + 1) % 3)) % 6)],
               scratch.arm1, tick1);
      stampArm(basal1[static_cast<size_t>((parity + 2 * ((arm + 2) % 3)) % 6)],
               scratch.arm2, tick2);
      const int seenTick = scratch.next(scratch.tickSeen, scratch.seen);
      const int pivot = basal1[static_cast<size_t>((parity + 2 * arm) % 6)];
      if (pivot < 0 || pivot >= nAtoms) {
        continue;
      }
      const auto &row = (*nList)[static_cast<size_t>(pivot)];
      for (size_t n = 1; n < row.size(); n++) {
        const int nb = row[n];
        if (nb < 0 || nb >= nAtoms ||
            scratch.onRing[static_cast<size_t>(nb)] == on) {
          continue;
        }
        const RingRow through = ringsOf(nb);
        for (int c = 0; c < through.n; c++) {
          const int candidate = through.p[c];
          if (scratch.seen[static_cast<size_t>(candidate)] == seenTick) {
            continue;
          }
          scratch.seen[static_cast<size_t>(candidate)] = seenTick;
          bool share = false;
          bool hit1 = false;
          bool hit2 = false;
          if (allHex) {
            const int *atoms = flat6.data() + static_cast<std::ptrdiff_t>(candidate) * 6;
            for (int t = 0; t < 6; t++) {
              const int atom = atoms[t];
              if (static_cast<unsigned>(atom) >= static_cast<unsigned>(nAtoms)) {
                continue;
              }
              if (scratch.onRing[static_cast<size_t>(atom)] == on) {
                share = true;
                break;
              }
              hit1 = hit1 || scratch.arm1[static_cast<size_t>(atom)] == tick1;
              hit2 = hit2 || scratch.arm2[static_cast<size_t>(atom)] == tick2;
            }
          } else {
            for (const int atom : (*rings)[static_cast<size_t>(candidate)]) {
              if (atom < 0 || atom >= nAtoms) {
                continue;
              }
              if (scratch.onRing[static_cast<size_t>(atom)] == on) {
                share = true;
                break;
              }
              hit1 = hit1 || scratch.arm1[static_cast<size_t>(atom)] == tick1;
              hit2 = hit2 || scratch.arm2[static_cast<size_t>(atom)] == tick2;
            }
          }
          if (!share && hit1 && hit2) {
            out.push_back(candidate);
          }
        }
      }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  mutable BasalScratch lazyScratch;

  [[nodiscard]] const std::vector<int> &basalCandidates(int i) const {
    if (!basalReady[static_cast<std::size_t>(i)]) {
      basalCand[static_cast<std::size_t>(i)] =
          collectBasalCandidates(i, lazyScratch);
      basalReady[static_cast<std::size_t>(i)] = 1;
    }
    return basalCand[static_cast<std::size_t>(i)];
  }

  //! Prismatic rings of the ordered passing pair (i, j)
  std::vector<int> prismaticOf(int i, int j) {
    return prismaticOf(i, j, scratchListHC, scratchType);
  }

  //! As above, on caller-owned scratch, so threads can share the context
  std::vector<int> prismaticOf(int i, int j, std::vector<int> &listHC,
                               std::vector<ring::strucType> &types) const {
    std::vector<int> prismatic;
    listHC.clear();
    (void)ring::findPrismatic(*rings, listHC, types, i, j, prismatic, index);
    return prismatic;
  }

  //! The equatorial-conditions test of findDDC without the claim skip; on a
  //! pass, fills every ring that is a peripheral in some valid assignment
  bool equatorialPass(int i, bool hcAffiliated,
                      std::vector<int> &peripherals) const {
    peripherals.clear();
    if (hcAffiliated) {
      return false;
    }
    return ring::equatorialConditions(*rings, i, index, &peripherals);
  }
};

} // namespace

namespace {

ring::CageAffiliation batchSweep(const FrameContext &ctx) {
  const std::vector<std::vector<int>> &rings = *ctx.rings;
  ring::CageAffiliation result;
  const size_t nRings = rings.size();
  result.hc.assign(nRings, false);
  result.ddc.assign(nRings, false);
  if (nRings == 0) {
    return result;
  }
  const int nR = static_cast<int>(nRings);
  // nvc++ -mp=gpu (SEAMS_HAS_OFFLOAD) SIGSEGVs host parallel-fors in libnvomp
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  const bool threaded = nR >= 1024 && !omp_in_parallel();
#else
  const bool threaded = false;
#endif

  // A flag only turns on. Once every ring is on, later pairs cannot change
  // the result, and on a ring-rich graph those pairs are most of the sweep.
  {
    std::vector<char> on(static_cast<std::size_t>(nR), 0);
    int left = nR;
    auto mark = [&](int m) {
      char &slot = on[static_cast<std::size_t>(m)];
      if (slot != 0) {
        return;
      }
      slot = 1;
      result.hc[static_cast<std::size_t>(m)] = true;
      --left;
    };
    std::vector<int> listHC;
    std::vector<ring::strucType> types(nRings, ring::strucType::unclassified);
    for (int i = 0; i < nR && left > 0; ++i) {
      for (const int j : ctx.basalCandidates(i)) {
        if (left == 0) {
          break;
        }
        if (!ctx.basalPair(i, j)) {
          continue;
        }
        mark(i);
        mark(j);
        for (const int k : ctx.prismaticOf(i, j, listHC, types)) {
          mark(k);
        }
      }
    }
  }

  // DDC affiliation: every ring passing the equatorial conditions marks
  // itself and its six peripherals
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel if (threaded)
#endif
  {
    std::vector<int> marks;
    std::vector<int> peripherals;
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp for schedule(dynamic, 256) nowait
#endif
    for (int i = 0; i < nR; i++) {
      if (ctx.equatorialPass(i, result.hc[static_cast<std::size_t>(i)],
                             peripherals)) {
        marks.push_back(i);
        marks.insert(marks.end(), peripherals.begin(), peripherals.end());
      }
    }
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp critical(seams_cage_affiliation)
#endif
    for (const int m : marks) {
      result.ddc[static_cast<std::size_t>(m)] = true;
    }
  }
  return result;
}

} // namespace

ring::CageAffiliation
ring::cageAffiliation(const std::vector<std::vector<int>> &rings,
                      const std::vector<std::vector<int>> &nList) {
  if (rings.empty()) {
    return batchSweep(FrameContext(rings, nList));
  }
  const FrameContext ctx(rings, nList);
  return batchSweep(ctx);
}

ring::StackingPlanes
ring::stackingPlanes(const std::vector<std::vector<int>> &rings,
                     const std::vector<std::vector<int>> &nList) {
  StackingPlanes out;
  const size_t nRings = rings.size();
  out.basal.assign(nRings, false);
  out.equatorial.assign(nRings, false);
  if (nRings == 0) {
    return out;
  }
  FrameContext ctx(rings, nList);
  std::vector<bool> hc(nRings, false);
  for (size_t i = 0; i < nRings; i++) {
    const std::vector<int> &candidates = ctx.basalCandidates(static_cast<int>(i));
    for (const int j : candidates) {
      if (!ctx.basalPair(static_cast<int>(i), j)) {
        continue;
      }
      out.basal[i] = true;
      out.basal[static_cast<std::size_t>(j)] = true;
      hc[i] = true;
      hc[static_cast<std::size_t>(j)] = true;
      for (const int k : ctx.prismaticOf(static_cast<int>(i), j)) {
        if (k >= 0 && static_cast<std::size_t>(k) < nRings) {
          hc[static_cast<std::size_t>(k)] = true;
        }
      }
    }
  }
  std::vector<int> peripherals;
  for (size_t i = 0; i < nRings; i++) {
    if (ctx.equatorialPass(static_cast<int>(i), hc[i], peripherals)) {
      out.equatorial[i] = true;
    }
  }
  return out;
}

// ---------------------------------------------------------------------------

/**
 * @details Implementation state of ring::AffiliationUpdater.
 *
 *  What each answer reads. The ordered basal pair (i, j) of batchSweep, and
 *  the prismatic rings it marks, read the rings i and j, the neighbour rows
 *  of the atoms of i, and the rings sharing an atom with i; a partner j
 *  bonds an atom of i. The equatorial test of i reads the rings sharing an
 *  atom with i, and the DDC sweep runs it only on rings without an HC flag.
 *  The flags are unions over the passing pairs and the passing equatorial
 *  rings. So the update stores, per ring, its passing pairs with their
 *  prismatic rings and its equatorial result, and evaluates again only the
 *  pairs of rings that are new, hold an atom whose row changed, or hold or
 *  bond an atom of an added or removed ring, and the equatorial test of
 *  rings that are new or hold an atom of an added or removed ring. Every
 *  other stored answer reads the same input in both frames and carries
 *  over; the flags are then the same unions, which equal the batch sweep.
 */
struct ring::AffiliationUpdater::Impl {
  bool primed = false;
  // Previous frame: rings flattened, neighbour rows sorted, and the answers
  std::vector<int> prevAtoms;
  std::vector<std::size_t> prevOffset;
  bool prevSorted = false;
  std::vector<std::vector<int>> prevSortedRows;
  std::vector<std::vector<int>> incoming;
  // False after a whole-frame cageAffiliation that did not store pairs.
  bool answersValid = false;
  // Once a frame has crossed into the whole-frame sweep, the next frame
  // skips the row diff. The diff plus a sweep was slower than the sweep
  // alone on the mid-radius shells.
  bool batchHot = false;
  // Reused across whole-frame passes. The cold cageAffiliation builds a
  // fresh context, so this only speeds the updater.
  std::unique_ptr<FrameContext> warm;

  void fastAffiliation(const std::vector<std::vector<int>> &rings,
                       const std::vector<std::vector<int>> &nList) {
    if (!warm) {
      warm = std::make_unique<FrameContext>();
    }
    warm->prepare(rings, nList);
    current = batchSweep(*warm);
    lastReclassified = static_cast<int>(rings.size());
  }
  // Passing pairs of ring i as runs [j, m, k_1 .. k_m] (partner, then its m
  // prismatic rings) in pairData[pairOffset[i]] up to pairOffset[i + 1]
  std::vector<int> pairData;
  std::vector<std::size_t> pairOffset;
  // Equatorial result per ring: -1 not evaluated (an HC ring), 0 fails,
  // 1 passes with its peripherals in periData[periOffset[i]..]
  std::vector<signed char> eq;
  std::vector<int> periData;
  std::vector<std::size_t> periOffset;
  CageAffiliation current;
  int lastReclassified = 0;
  double batchFraction = 0.25;
  // The back-off of RingUpdater: after a frame that went to batch on churn,
  // the next `skip` frames run cageAffiliation with no diff state, the wait
  // doubling from 1 to kMaxBackoff while frames keep going to batch
  int backoff = 0;
  int skip = 0;
  static constexpr int kMaxBackoff = 8;

  void fellBack() {
    backoff = std::min(std::max(1, 2 * backoff), kMaxBackoff);
    skip = backoff;
  }

  void release() {
    std::vector<std::vector<int>>().swap(prevSortedRows);
    std::vector<int>().swap(prevAtoms);
    std::vector<std::size_t>().swap(prevOffset);
    std::vector<int>().swap(pairData);
    std::vector<std::size_t>().swap(pairOffset);
    std::vector<signed char>().swap(eq);
    std::vector<int>().swap(periData);
    std::vector<std::size_t>().swap(periOffset);
    primed = false;
  }

  // Filled in place. The outer vector and each row keep their buffers.
  void fillSorted(const std::vector<std::vector<int>> &nList,
                  std::vector<std::vector<int>> &out) const {
    const int n = static_cast<int>(nList.size());
    if (static_cast<int>(out.size()) != n) {
      out.assign(static_cast<std::size_t>(n), {});
    }
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel for schedule(static) if (n >= 1024 && !omp_in_parallel())
#endif
    for (int a = 0; a < n; a++) {
      auto &row = out[static_cast<std::size_t>(a)];
      row.assign(nList[static_cast<std::size_t>(a)].begin(),
                 nList[static_cast<std::size_t>(a)].end());
      std::sort(row.begin(), row.end());
    }
  }

  // The passing pairs (i, j) of batchSweep with their prismatic rings
  static void pairsOf(const FrameContext &ctx, int i, std::vector<int> &listHC,
                      std::vector<ring::strucType> &types,
                      std::vector<int> &out, FrameContext::BasalScratch &scratch) {
    out.clear();
    for (const int j : ctx.collectBasalCandidates(i, scratch)) {
      if (!ctx.basalPair(i, j)) {
        continue;
      }
      const std::vector<int> prismatic = ctx.prismaticOf(i, j, listHC, types);
      out.push_back(j);
      out.push_back(static_cast<int>(prismatic.size()));
      out.insert(out.end(), prismatic.begin(), prismatic.end());
    }
  }

  // Evaluate the pairs of the rings in `which` into per-ring runs
  static std::vector<std::vector<int>>
  evaluatePairs(const FrameContext &ctx, const std::vector<int> &which) {
    const int n = static_cast<int>(which.size());
    std::vector<std::vector<int>> out(which.size());
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel if (!omp_in_parallel())
#endif
    {
      std::vector<int> listHC;
      std::vector<ring::strucType> types(ctx.rings->size(),
                                         ring::strucType::unclassified);
      FrameContext::BasalScratch scratch;
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp for schedule(dynamic, 64)
#endif
      for (int k = 0; k < n; k++) {
        pairsOf(ctx, which[static_cast<std::size_t>(k)], listHC, types,
                out[static_cast<std::size_t>(k)], scratch);
      }
    }
    return out;
  }

  // HC flags from the stored pairs; then the equatorial test of every ring
  // without one that has no stored result, and the DDC flags
  void flags(const FrameContext &ctx, int &evaluated) {
    const std::size_t nR = ctx.rings->size();
    current.hc.assign(nR, false);
    current.ddc.assign(nR, false);
    for (std::size_t i = 0; i < nR; i++) {
      std::size_t p = pairOffset[i];
      const std::size_t end = pairOffset[i + 1];
      if (p < end) {
        current.hc[i] = true;
      }
      while (p < end) {
        current.hc[static_cast<std::size_t>(pairData[p])] = true;
        const int m = pairData[p + 1];
        for (int t = 0; t < m; t++) {
          current.hc[static_cast<std::size_t>(pairData[p + 2 + t])] = true;
        }
        p += 2 + static_cast<std::size_t>(m);
      }
    }
    std::vector<int> pending;
    for (std::size_t i = 0; i < nR; i++) {
      if (!current.hc[i] && eq[i] < 0) {
        pending.push_back(static_cast<int>(i));
      }
    }
    const int n = static_cast<int>(pending.size());
    std::vector<std::vector<int>> found(pending.size());
    std::vector<signed char> pass(pending.size(), 0);
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
#pragma omp parallel for schedule(guided) if (!omp_in_parallel())
#endif
    for (int k = 0; k < n; k++) {
      pass[static_cast<std::size_t>(k)] =
          ctx.equatorialPass(pending[static_cast<std::size_t>(k)], false,
                             found[static_cast<std::size_t>(k)])
              ? 1
              : 0;
    }
    evaluated += n;
    if (n > 0) {
      // Splice the new peripheral lists in ring order
      std::vector<int> data;
      std::vector<std::size_t> offset(nR + 1, 0);
      std::size_t k = 0;
      for (std::size_t i = 0; i < nR; i++) {
        if (k < pending.size() &&
            static_cast<std::size_t>(pending[k]) == i) {
          eq[i] = pass[k];
          if (pass[k] != 0) {
            data.insert(data.end(), found[k].begin(), found[k].end());
          }
          ++k;
        } else {
          data.insert(data.end(),
                      periData.begin() +
                          static_cast<std::ptrdiff_t>(periOffset[i]),
                      periData.begin() +
                          static_cast<std::ptrdiff_t>(periOffset[i + 1]));
        }
        offset[i + 1] = data.size();
      }
      periData.swap(data);
      periOffset.swap(offset);
    }
    for (std::size_t i = 0; i < nR; i++) {
      if (current.hc[i] || eq[i] != 1) {
        continue;
      }
      current.ddc[i] = true;
      for (std::size_t p = periOffset[i]; p < periOffset[i + 1]; p++) {
        current.ddc[static_cast<std::size_t>(periData[p])] = true;
      }
    }
  }

  // Whole-frame evaluation that also fills the stored answers
  void batch(const std::vector<std::vector<int>> &rings,
             const std::vector<std::vector<int>> &nList,
             std::vector<std::vector<int>> rows) {
    const FrameContext ctx(rings, nList);
    const int nR = static_cast<int>(rings.size());
    std::vector<int> all(rings.size());
    std::iota(all.begin(), all.end(), 0);
    const auto fresh = evaluatePairs(ctx, all);
    pairOffset.assign(rings.size() + 1, 0);
    pairData.clear();
    for (int i = 0; i < nR; i++) {
      pairData.insert(pairData.end(), fresh[static_cast<std::size_t>(i)].begin(),
                      fresh[static_cast<std::size_t>(i)].end());
      pairOffset[static_cast<std::size_t>(i) + 1] = pairData.size();
    }
    eq.assign(rings.size(), -1);
    periData.clear();
    periOffset.assign(rings.size() + 1, 0);
    int evaluated = 0;
    flags(ctx, evaluated);
    lastReclassified = nR;
    remember(rings, std::move(rows));
    primed = true;
  }

  void remember(const std::vector<std::vector<int>> &rings,
                std::vector<std::vector<int>> rows) {
    prevOffset.assign(rings.size() + 1, 0);
    for (std::size_t i = 0; i < rings.size(); i++) {
      prevOffset[i + 1] = prevOffset[i] + rings[i].size();
    }
    prevAtoms.resize(prevOffset.back());
    for (std::size_t i = 0; i < rings.size(); i++) {
      std::copy(rings[i].begin(), rings[i].end(),
                prevAtoms.begin() + static_cast<std::ptrdiff_t>(prevOffset[i]));
    }
    prevSorted = std::is_sorted(rings.begin(), rings.end());
    prevSortedRows = std::move(rows);
  }

  // prevOf[i]: the previous index of ring i, or -1 for an added ring. Two
  // sorted lists (the order RingUpdater returns) are matched by a merge on
  // equal atom sequences; otherwise by a fingerprint of the canonical cycle,
  // confirmed atom by atom on a match.
  std::vector<int> matchRings(const std::vector<std::vector<int>> &rings) const {
    const std::size_t nPrev = prevOffset.size() - 1;
    std::vector<int> prevOf(rings.size(), -1);
    if (prevSorted && std::is_sorted(rings.begin(), rings.end())) {
      std::size_t p = 0;
      for (std::size_t q = 0; q < rings.size() && p < nPrev;) {
        const int *a = prevAtoms.data() + prevOffset[p];
        const int *aEnd = prevAtoms.data() + prevOffset[p + 1];
        const auto &b = rings[q];
        if (std::lexicographical_compare(a, aEnd, b.begin(), b.end())) {
          ++p;
        } else if (std::lexicographical_compare(b.begin(), b.end(), a, aEnd)) {
          ++q;
        } else {
          prevOf[q] = static_cast<int>(p);
          ++p;
          ++q;
        }
      }
      return prevOf;
    }
    std::vector<std::pair<std::uint64_t, int>> prevFP(nPrev);
    std::vector<std::pair<std::uint64_t, int>> newFP(rings.size());
    for (std::size_t i = 0; i < nPrev; i++) {
      prevFP[i] = {cycleHash(prevAtoms.data() + prevOffset[i],
                             prevOffset[i + 1] - prevOffset[i]),
                   static_cast<int>(i)};
    }
    for (std::size_t i = 0; i < rings.size(); i++) {
      newFP[i] = {cycleHash(rings[i].data(), rings[i].size()),
                  static_cast<int>(i)};
    }
    std::sort(prevFP.begin(), prevFP.end());
    std::sort(newFP.begin(), newFP.end());
    std::vector<char> prevMatched(nPrev, 0);
    for (std::size_t p = 0, q = 0; p < nPrev && q < newFP.size();) {
      if (prevFP[p].first < newFP[q].first) {
        ++p;
        continue;
      }
      if (newFP[q].first < prevFP[p].first) {
        ++q;
        continue;
      }
      std::size_t pEnd = p;
      std::size_t qEnd = q;
      while (pEnd < nPrev && prevFP[pEnd].first == prevFP[p].first) {
        ++pEnd;
      }
      while (qEnd < newFP.size() && newFP[qEnd].first == newFP[q].first) {
        ++qEnd;
      }
      for (std::size_t y = q; y < qEnd; y++) {
        const auto &r = rings[static_cast<std::size_t>(newFP[y].second)];
        for (std::size_t x = p; x < pEnd; x++) {
          const int i = prevFP[x].second;
          if (!prevMatched[static_cast<std::size_t>(i)] &&
              sameCycle(r.data(), r.size(),
                        prevAtoms.data() + prevOffset[i],
                        prevOffset[i + 1] - prevOffset[i])) {
            prevMatched[static_cast<std::size_t>(i)] = 1;
            prevOf[static_cast<std::size_t>(newFP[y].second)] = i;
            break;
          }
        }
      }
      p = pEnd;
      q = qEnd;
    }
    return prevOf;
  }
};

ring::AffiliationUpdater::AffiliationUpdater()
    : impl_(std::make_unique<Impl>()) {}
ring::AffiliationUpdater::~AffiliationUpdater() = default;
ring::AffiliationUpdater::AffiliationUpdater(AffiliationUpdater &&) noexcept =
    default;
ring::AffiliationUpdater &
ring::AffiliationUpdater::operator=(AffiliationUpdater &&) noexcept = default;

int ring::AffiliationUpdater::lastReclassified() const {
  return impl_->lastReclassified;
}

void ring::AffiliationUpdater::setBatchFraction(double fraction) {
  impl_->batchFraction = fraction;
}

const ring::CageAffiliation &
ring::AffiliationUpdater::update(const std::vector<std::vector<int>> &rings,
                                 const std::vector<std::vector<int>> &nList) {
  Impl &st = *impl_;
  if (rings.empty()) {
    // No six-ring to classify. A full context build still walks the
    // neighbour list; the flags are empty either way. Drop the diff state
    // so the next non-empty frame starts from a clean batch.
    st.current.hc.clear();
    st.current.ddc.clear();
    st.lastReclassified = 0;
    st.batchHot = false;
    st.primed = false;
    st.answersValid = false;
    return st.current;
  }
  if (st.skip > 0) {
    // A plain affiliation, keeping the sorted-row buffers for the next diff.
    st.skip--;
    st.current = ring::cageAffiliation(rings, nList);
    st.lastReclassified = static_cast<int>(rings.size());
    st.answersValid = false;
    st.batchHot = false;
    st.fillSorted(nList, st.incoming);
    st.prevSortedRows.swap(st.incoming);
    return st.current;
  }
  if (st.batchHot && st.prevSortedRows.size() == nList.size()) {
    st.fastAffiliation(rings, nList);
    return st.current;
  }
  st.fillSorted(nList, st.incoming);
  auto &rows = st.incoming;
  if (!st.primed || st.prevSortedRows.size() != nList.size()) {
    st.batchHot = false;
    st.batch(rings, nList, rows);
    st.answersValid = true;
    return st.current;
  }
  const int nR = static_cast<int>(rings.size());
  const int nAtoms = static_cast<int>(nList.size());
  const auto batchLimit = static_cast<std::size_t>(st.batchFraction * nR);

  std::vector<char> changedRow(static_cast<std::size_t>(nAtoms), 0);
  std::size_t nChanged = 0;
  for (int a = 0; a < nAtoms; a++) {
    if (rows[static_cast<std::size_t>(a)] !=
        st.prevSortedRows[static_cast<std::size_t>(a)]) {
      changedRow[static_cast<std::size_t>(a)] = 1;
      ++nChanged;
    }
  }
  // Every ring through a changed row is evaluated again: past the limit, go
  // to batch before matching rings or building any index
  std::size_t throughChanged = 0;
  if (static_cast<double>(nChanged) * nR <=
      static_cast<double>(batchLimit) * nAtoms) {
    for (const auto &r : rings) {
      for (const int a : r) {
        if (a < 0 || a >= nAtoms || changedRow[static_cast<std::size_t>(a)]) {
          ++throughChanged;
          break;
        }
      }
      if (throughChanged > batchLimit) {
        break;
      }
    }
  }
  if (static_cast<double>(nChanged) * nR >
          static_cast<double>(batchLimit) * nAtoms ||
      throughChanged > batchLimit) {
    // The batch pass stores every pair. On a frame this large the stored
    // pairs are not reused, and building them is slower than cageAffiliation.
    st.fastAffiliation(rings, nList);
    st.prevSortedRows.swap(st.incoming);
    st.answersValid = false;
    st.batchHot = true;
    return st.current;
  }
  if (!st.answersValid) {
    st.batch(rings, nList, rows);
    st.answersValid = true;
    return st.current;
  }

  const std::vector<int> prevOf = st.matchRings(rings);
  const std::size_t nPrev = st.prevOffset.size() - 1;
  std::vector<int> newOf(nPrev, -1);
  for (int i = 0; i < nR; i++) {
    if (prevOf[static_cast<std::size_t>(i)] >= 0) {
      newOf[static_cast<std::size_t>(prevOf[static_cast<std::size_t>(i)])] = i;
    }
  }
  if (nChanged == 0 && nPrev == rings.size() &&
      std::find(prevOf.begin(), prevOf.end(), -1) == prevOf.end()) {
    bool same = true;
    for (int i = 0; i < nR && same; i++) {
      same = prevOf[static_cast<std::size_t>(i)] == i;
    }
    if (same) {
      st.lastReclassified = 0;
      st.backoff = 0;
      return st.current;
    }
  }

  // Atoms of added and removed rings
  std::vector<char> moved(static_cast<std::size_t>(nAtoms), 0);
  auto markAtom = [&](int a) {
    if (a >= 0 && a < nAtoms) {
      moved[static_cast<std::size_t>(a)] = 1;
    }
  };
  for (int i = 0; i < nR; i++) {
    if (prevOf[static_cast<std::size_t>(i)] < 0) {
      for (const int a : rings[static_cast<std::size_t>(i)]) {
        markAtom(a);
      }
    }
  }
  for (std::size_t p = 0; p < nPrev; p++) {
    if (newOf[p] < 0) {
      for (std::size_t k = st.prevOffset[p]; k < st.prevOffset[p + 1]; k++) {
        markAtom(st.prevAtoms[k]);
      }
    }
  }

  std::vector<int> dirtyPairs;
  std::vector<char> dirtyEq(static_cast<std::size_t>(nR), 0);
  std::vector<char> isDirty(static_cast<std::size_t>(nR), 0);
  for (int i = 0; i < nR; i++) {
    const auto &r = rings[static_cast<std::size_t>(i)];
    bool pairs = prevOf[static_cast<std::size_t>(i)] < 0;
    bool equatorial = pairs;
    for (const int a : r) {
      if (a < 0 || a >= nAtoms) {
        pairs = equatorial = true;
        break;
      }
      if (moved[static_cast<std::size_t>(a)]) {
        pairs = equatorial = true;
        break;
      }
      if (changedRow[static_cast<std::size_t>(a)]) {
        pairs = true;
        continue;
      }
      if (!pairs) {
        const auto &row = nList[static_cast<std::size_t>(a)];
        for (std::size_t m = 1; m < row.size(); m++) {
          const int b = row[m];
          if (b >= 0 && b < nAtoms && moved[static_cast<std::size_t>(b)]) {
            pairs = true;
            break;
          }
        }
      }
    }
    if (pairs) {
      dirtyPairs.push_back(i);
      isDirty[static_cast<std::size_t>(i)] = 1;
    }
    dirtyEq[static_cast<std::size_t>(i)] = equatorial ? 1 : 0;
  }
  if (dirtyPairs.size() > batchLimit) {
    // The next frames skip the diff. The call has to run before the return.
    st.fellBack();
    st.fastAffiliation(rings, nList);
    st.prevSortedRows.swap(st.incoming);
    st.answersValid = false;
    st.batchHot = true;
    return st.current;
  }

  // Carried answers name previous ring indices: map them, and evaluate a
  // ring again if one of them has no image
  std::vector<int> pairData;
  std::vector<std::size_t> pairOffset(rings.size() + 1, 0);
  for (int i = 0; i < nR; i++) {
    if (isDirty[static_cast<std::size_t>(i)]) {
      pairOffset[static_cast<std::size_t>(i) + 1] = pairData.size();
      continue;
    }
    const auto was = static_cast<std::size_t>(prevOf[static_cast<std::size_t>(i)]);
    for (std::size_t p = st.pairOffset[was]; p < st.pairOffset[was + 1];) {
      const int j = newOf[static_cast<std::size_t>(st.pairData[p])];
      const int m = st.pairData[p + 1];
      pairData.push_back(j);
      pairData.push_back(m);
      for (int t = 0; t < m; t++) {
        const int k = newOf[static_cast<std::size_t>(st.pairData[p + 2 + t])];
        pairData.push_back(k);
      }
      p += 2 + static_cast<std::size_t>(m);
    }
    pairOffset[static_cast<std::size_t>(i) + 1] = pairData.size();
  }

  const FrameContext ctx(rings, nList);
  const auto fresh = Impl::evaluatePairs(ctx, dirtyPairs);
  {
    // Merge the fresh runs of the dirty rings into ring order
    std::vector<int> data;
    data.reserve(pairData.size() + 16 * dirtyPairs.size());
    std::vector<std::size_t> offset(rings.size() + 1, 0);
    std::size_t k = 0;
    for (int i = 0; i < nR; i++) {
      if (k < dirtyPairs.size() && dirtyPairs[k] == i) {
        data.insert(data.end(), fresh[k].begin(), fresh[k].end());
        ++k;
      } else {
        data.insert(data.end(),
                    pairData.begin() + static_cast<std::ptrdiff_t>(
                                           pairOffset[static_cast<std::size_t>(i)]),
                    pairData.begin() + static_cast<std::ptrdiff_t>(
                                           pairOffset[static_cast<std::size_t>(i) + 1]));
      }
      offset[static_cast<std::size_t>(i) + 1] = data.size();
    }
    st.pairData.swap(data);
    st.pairOffset.swap(offset);
  }

  // Equatorial answers carry over unless the ring is dirty for them
  std::vector<signed char> eq(rings.size(), -1);
  std::vector<int> periData;
  std::vector<std::size_t> periOffset(rings.size() + 1, 0);
  for (int i = 0; i < nR; i++) {
    const auto ii = static_cast<std::size_t>(i);
    if (!dirtyEq[ii]) {
      const auto was = static_cast<std::size_t>(prevOf[ii]);
      eq[ii] = st.eq[was];
      for (std::size_t p = st.periOffset[was]; p < st.periOffset[was + 1]; p++) {
        const int k = newOf[static_cast<std::size_t>(st.periData[p])];
        periData.push_back(k);
      }
    }
    periOffset[ii + 1] = periData.size();
  }
  st.eq.swap(eq);
  st.periData.swap(periData);
  st.periOffset.swap(periOffset);

  int evaluated = 0;
  st.flags(ctx, evaluated);
  // The equatorial flag is set only together with the pairs flag,
  // so the pairs flag already counts every ring that test will rerun.
  std::size_t touched = 0;
  for (int i = 0; i < nR; i++) {
    touched += isDirty[static_cast<std::size_t>(i)] ? 1 : 0;
  }
  st.lastReclassified = static_cast<int>(touched);
  st.remember(rings, rows);
  st.answersValid = true;
  st.backoff = 0;
  return st.current;
}

ring::SeededAtomLabels ring::seededCageAffiliation(
    const std::vector<std::vector<int>> &strictRings,
    const std::vector<std::vector<int>> &strictNList,
    const std::vector<std::vector<int>> &permissiveRings,
    const std::vector<std::vector<int>> &permissiveNList,
    bool ringAdjacentCompletion) {
  return seededCageAffiliation(
      strictRings, cageAffiliation(strictRings, strictNList), permissiveRings,
      cageAffiliation(permissiveRings, permissiveNList), permissiveNList,
      ringAdjacentCompletion);
}

ring::SeededAtomLabels ring::seededCageAffiliation(
    const std::vector<std::vector<int>> &strictRings,
    const CageAffiliation &strictAffiliation,
    const std::vector<std::vector<int>> &permissiveRings,
    const CageAffiliation &permissiveAffiliation,
    const std::vector<std::vector<int>> &permissiveNList,
    bool ringAdjacentCompletion) {
  const int nAtoms = static_cast<int>(permissiveNList.size());
  SeededAtomLabels out;
  out.hc.assign(nAtoms, false);
  out.ddc.assign(nAtoms, false);
  if (nAtoms == 0) {
    return out;
  }

  const auto toAtoms = [nAtoms](const std::vector<std::vector<int>> &rings,
                                const CageAffiliation &affiliation) {
    std::pair<std::vector<bool>, std::vector<bool>> flags{
        std::vector<bool>(nAtoms, false), std::vector<bool>(nAtoms, false)};
    for (size_t i = 0; i < rings.size(); i++) {
      for (const int atom : rings[i]) {
        if (atom < 0 || atom >= nAtoms) {
          continue;
        }
        if (affiliation.hc[i]) {
          flags.first[atom] = true;
        }
        if (affiliation.ddc[i]) {
          flags.second[atom] = true;
        }
      }
    }
    return flags;
  };

  const auto [hcS, ddcS] = toAtoms(strictRings, strictAffiliation);
  const auto [hcP, ddcP] = toAtoms(permissiveRings, permissiveAffiliation);

  // One flood per cage type: an HC seed must not keep a DDC-only atom
  // that happens to sit in the same permissive H-bond component.
  auto flood = [&](const std::vector<bool> &seeds,
                   const std::vector<bool> &affiliated) {
    std::vector<bool> kept(nAtoms, false);
    std::vector<int> stack;
    for (int a = 0; a < nAtoms; a++) {
      if (!seeds[a] || kept[a]) {
        continue;
      }
      stack.push_back(a);
      kept[a] = true;
      while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        if (u < 0 || static_cast<size_t>(u) >= permissiveNList.size()) {
          continue;
        }
        for (size_t m = 1; m < permissiveNList[u].size(); m++) {
          const int v = permissiveNList[u][m];
          if (v >= 0 && v < nAtoms && affiliated[v] && !kept[v]) {
            kept[v] = true;
            stack.push_back(v);
          }
        }
      }
    }
    return kept;
  };

  std::vector<bool> affiliatedHC(nAtoms, false);
  std::vector<bool> affiliatedDDC(nAtoms, false);
  for (int a = 0; a < nAtoms; a++) {
    affiliatedHC[a] = hcS[a] || hcP[a];
    affiliatedDDC[a] = ddcS[a] || ddcP[a];
  }
  const auto keptHC = flood(hcS, affiliatedHC);
  const auto keptDDC = flood(ddcS, affiliatedDDC);

  for (int a = 0; a < nAtoms; a++) {
    if (keptHC[a]) {
      const bool fromStrict = hcS[a];
      const bool fromPermissive = hcP[a];
      if (fromStrict || fromPermissive) {
        out.hc[a] = true;
      }
    }
    if (keptDDC[a]) {
      const bool fromStrict = ddcS[a];
      const bool fromPermissive = ddcP[a];
      if (fromStrict || fromPermissive) {
        out.ddc[a] = true;
      }
    }
  }
  if (ringAdjacentCompletion) {
    return ring::ringAdjacentCompletion(out, permissiveRings);
  }
  return out;
}

ring::SeededAtomLabels
ring::ringAdjacentCompletion(const SeededAtomLabels &labels,
                             const std::vector<std::vector<int>> &rings) {
  SeededAtomLabels out = labels;
  const int nAtoms = static_cast<int>(out.hc.size());
  const int nRings = static_cast<int>(rings.size());
  if (nAtoms == 0 || nRings == 0) {
    return out;
  }
  // Rings through each atom, so that accepting one vertex only revisits
  // the rings it belongs to
  std::vector<std::vector<int>> ringsOf(static_cast<std::size_t>(nAtoms));
  for (int r = 0; r < nRings; r++) {
    for (const int a : rings[r]) {
      if (a >= 0 && a < nAtoms) {
        ringsOf[static_cast<std::size_t>(a)].push_back(r);
      }
    }
  }
  // A ring with every vertex but one labelled is a cage ring with a
  // vacancy in the label, not a liquid ring at the interface: a liquid
  // six-ring touching a nucleus has at most a few labelled vertices.
  // Requiring all but one keeps the completion from walking into the
  // liquid, where a single shared edge would let it flood.
  auto complete = [&](std::vector<bool> &flag) {
    std::vector<int> frontier;
    frontier.reserve(static_cast<std::size_t>(nRings));
    for (int r = 0; r < nRings; r++) {
      frontier.push_back(r);
    }
    std::vector<char> queued(static_cast<std::size_t>(nRings), 1);
    while (!frontier.empty()) {
      const int r = frontier.back();
      frontier.pop_back();
      queued[static_cast<std::size_t>(r)] = 0;
      int missing = -1;
      int nMissing = 0;
      for (const int a : rings[r]) {
        if (a < 0 || a >= nAtoms) {
          nMissing = 2;
          break;
        }
        if (!flag[static_cast<std::size_t>(a)]) {
          ++nMissing;
          missing = a;
        }
      }
      if (nMissing != 1) {
        continue;
      }
      flag[static_cast<std::size_t>(missing)] = true;
      for (const int s : ringsOf[static_cast<std::size_t>(missing)]) {
        if (!queued[static_cast<std::size_t>(s)]) {
          queued[static_cast<std::size_t>(s)] = 1;
          frontier.push_back(s);
        }
      }
    }
  };
  complete(out.hc);
  complete(out.ddc);
  return out;
}
