//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <tum_offload.hpp>

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <tum_device.hpp>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <vector>

#ifdef SEAMS_HAS_OFFLOAD
#include <omp.h>
#endif

namespace {

struct FlatGraph {
  int nAtoms = 0;
  int kMax = 0;
  std::vector<int> deg;
  std::vector<int> cols;
};

FlatGraph flattenIndexList(const std::vector<std::vector<int>> &nList) {
  FlatGraph g;
  g.nAtoms = static_cast<int>(nList.size());
  int kMax = 0;
  for (const auto &row : nList) {
    const int deg = row.empty() ? 0 : static_cast<int>(row.size()) - 1;
    kMax = std::max(kMax, deg);
  }
  g.kMax = std::max(kMax, 1);
  g.deg.assign(static_cast<std::size_t>(g.nAtoms), 0);
  g.cols.assign(static_cast<std::size_t>(g.nAtoms) *
                    static_cast<std::size_t>(g.kMax),
                -1);
  for (int i = 0; i < g.nAtoms; ++i) {
    const auto &row = nList[static_cast<std::size_t>(i)];
    int kept = 0;
    for (std::size_t t = 1; t < row.size() && kept < g.kMax; ++t) {
      g.cols[static_cast<std::size_t>(i * g.kMax + kept)] = row[t];
      ++kept;
    }
    g.deg[static_cast<std::size_t>(i)] = kept;
  }
  return g;
}

void tallyAtoms(tum::CageCounts &out) {
  out.nHcAtoms = 0;
  out.nDdcAtoms = 0;
  for (int v : out.atomHc) {
    out.nHcAtoms += v;
  }
  for (int v : out.atomDdc) {
    out.nDdcAtoms += v;
  }
}

// Buffer sizes of one pass.
struct Caps {
  int rowCap = 0;
  int ringCap = 0;
  int pairCap = 0;
};

// What a pass needed: rings and basal pairs as counted, and the longest
// atom row before the clamp. A stage after the first overflow is skipped
// and reports 0.
struct Need {
  int nRings = 0;
  int rowMax = 0;
  int nPairs = 0;
};

bool fits(const Need &need, const Caps &caps) {
  return need.nRings <= caps.ringCap && need.rowMax <= caps.rowCap &&
         need.nPairs <= caps.pairCap;
}

int capProduct(long long a, long long b, long long limit) {
  const long long v = std::max(1LL, a * b);
  return static_cast<int>(std::min(v, limit));
}

// One pass over fixed buffers. out holds the labels only when the
// returned Need fits caps.
Need onePass(const FlatGraph &g, bool target, const Caps &caps,
             tum::CageCounts &out) {
  const int nAtoms = g.nAtoms;
  const int kMax = g.kMax;
  const int maxPer = caps.rowCap;
  const int maxRings = caps.ringCap;
  const int maxPairs = caps.pairCap;
  std::vector<int> nRings(1, 0);
  std::vector<int> nPairs(1, 0);
  std::vector<int> ringAtoms(static_cast<std::size_t>(maxRings) * 6, -1);
  std::vector<int> throughCount(static_cast<std::size_t>(nAtoms), 0);
  std::vector<int> through(static_cast<std::size_t>(nAtoms) *
                               static_cast<std::size_t>(maxPer),
                           -1);
  std::vector<int> hc(static_cast<std::size_t>(maxRings), 0);
  std::vector<int> ddc(static_cast<std::size_t>(maxRings), 0);
  std::vector<int> pairs(static_cast<std::size_t>(maxPairs) * 2, -1);
  out.atomHc.assign(static_cast<std::size_t>(nAtoms), 0);
  out.atomDdc.assign(static_cast<std::size_t>(nAtoms), 0);

  const int *deg = g.deg.data();
  const int *cols = g.cols.data();
  int *nRingsP = nRings.data();
  int *ringAtomsP = ringAtoms.data();
  int *throughCountP = throughCount.data();
  int *throughP = through.data();
  int *hcP = hc.data();
  int *ddcP = ddc.data();
  int *nPairsP = nPairs.data();
  int *pairsP = pairs.data();
  int *atomHcP = out.atomHc.data();
  int *atomDdcP = out.atomDdc.data();
  int rowMax = 0;

#ifdef SEAMS_HAS_OFFLOAD
  if (target) {
    const int degN = nAtoms;
    const int colN = nAtoms * kMax;
    const int ringN = maxRings * 6;
    const int thruN = nAtoms * maxPer;
    const int pairN = maxPairs * 2;
#pragma omp target data map(to : deg[0 : degN], cols[0 : colN], nAtoms, kMax,  \
                                maxRings, maxPer, maxPairs)                    \
    map(tofrom : nRingsP[0 : 1], nPairsP[0 : 1], ringAtomsP[0 : ringN],        \
            throughCountP[0 : nAtoms], throughP[0 : thruN], hcP[0 : maxRings], \
            ddcP[0 : maxRings], pairsP[0 : pairN], atomHcP[0 : nAtoms],        \
            atomDdcP[0 : nAtoms])
    {
#pragma omp target teams distribute parallel for
      for (int i = 0; i < nAtoms; ++i) {
        tum::device::enumSixFrom(i, deg, cols, nAtoms, kMax, maxRings, nRingsP,
                                 ringAtomsP);
      }
      // Kernels no-op past nRings/nPairs. The host reads a counter only
      // after a target update.
#pragma omp target update from(nRingsP[0 : 1])
      bool more = nRingsP[0] <= maxRings;
      if (more) {
#pragma omp target teams distribute parallel for
        for (int r = 0; r < maxRings; ++r) {
          tum::device::invertOneRing(r, nRingsP, ringAtomsP, nAtoms, maxPer,
                                     throughCountP, throughP);
        }
#pragma omp target teams distribute parallel for reduction(max : rowMax)     \
    map(tofrom : rowMax)
        for (int a = 0; a < nAtoms; ++a) {
          int n = throughCountP[a];
          rowMax = n > rowMax ? n : rowMax;
          if (n > maxPer) {
            n = maxPer;
          }
          throughCountP[a] = n;
          tum::device::sortThroughRow(throughP + a * maxPer, n);
        }
        more = rowMax <= maxPer;
      }
      if (more) {
#pragma omp target teams distribute parallel for
        for (int r = 0; r < maxRings; ++r) {
          tum::device::emitBasalFrom(r, nRingsP, ringAtomsP, deg, cols,
                                     throughCountP, throughP, nAtoms, kMax,
                                     maxPer, maxPairs, nPairsP, pairsP);
        }
#pragma omp target update from(nPairsP[0 : 1])
        more = nPairsP[0] <= maxPairs;
      }
      if (more) {
#pragma omp target teams distribute parallel for
        for (int p = 0; p < maxPairs; ++p) {
          tum::device::applyHcPair(p, nPairsP, pairsP, ringAtomsP,
                                   throughCountP, throughP, nAtoms, maxPer,
                                   hcP);
        }
#pragma omp target teams distribute parallel for
        for (int r = 0; r < maxRings; ++r) {
          tum::device::ddcFrom(r, nRingsP, ringAtomsP, throughCountP, throughP,
                               hcP, nAtoms, maxPer, ddcP);
        }
#pragma omp target teams distribute parallel for
        for (int r = 0; r < maxRings; ++r) {
          tum::device::atomIceFrom(r, nRingsP, ringAtomsP, hcP, ddcP, nAtoms,
                                   atomHcP, atomDdcP);
        }
      }
    }
  } else
#else
  (void)target;
#endif
  {
    for (int i = 0; i < nAtoms; ++i) {
      tum::device::enumSixFrom(i, deg, cols, nAtoms, kMax, maxRings, nRingsP,
                               ringAtomsP);
    }
    bool more = nRings[0] <= maxRings;
    if (more) {
      for (int r = 0; r < nRings[0]; ++r) {
        tum::device::invertOneRing(r, nRingsP, ringAtomsP, nAtoms, maxPer,
                                   throughCountP, throughP);
      }
      for (int a = 0; a < nAtoms; ++a) {
        int n = throughCount[static_cast<std::size_t>(a)];
        rowMax = std::max(rowMax, n);
        n = std::min(n, maxPer);
        throughCount[static_cast<std::size_t>(a)] = n;
        tum::device::sortThroughRow(throughP + a * maxPer, n);
      }
      more = rowMax <= maxPer;
    }
    if (more) {
      for (int r = 0; r < nRings[0]; ++r) {
        tum::device::emitBasalFrom(r, nRingsP, ringAtomsP, deg, cols,
                                   throughCountP, throughP, nAtoms, kMax,
                                   maxPer, maxPairs, nPairsP, pairsP);
      }
      more = nPairs[0] <= maxPairs;
    }
    if (more) {
      for (int p = 0; p < nPairs[0]; ++p) {
        tum::device::applyHcPair(p, nPairsP, pairsP, ringAtomsP, throughCountP,
                                 throughP, nAtoms, maxPer, hcP);
      }
      for (int r = 0; r < nRings[0]; ++r) {
        tum::device::ddcFrom(r, nRingsP, ringAtomsP, throughCountP, throughP,
                             hcP, nAtoms, maxPer, ddcP);
      }
      for (int r = 0; r < nRings[0]; ++r) {
        tum::device::atomIceFrom(r, nRingsP, ringAtomsP, hcP, ddcP, nAtoms,
                                 atomHcP, atomDdcP);
      }
    }
  }

  const Need need{nRings[0], rowMax, nPairs[0]};
  if (!fits(need, caps)) {
    return need;
  }
  out.nSix = nRings[0];
  out.nHcRings = 0;
  out.nDdcRings = 0;
  for (int r = 0; r < out.nSix; ++r) {
    out.nHcRings += hc[static_cast<std::size_t>(r)];
    out.nDdcRings += ddc[static_cast<std::size_t>(r)];
  }
  out.atomSix = std::move(throughCount);
  tallyAtoms(out);
  return need;
}

// Rerun with every overflowed buffer grown to what the pass counted.
// Caps only grow and never past the true sizes, so the loop ends.
tum::CageCounts runSpecialized(const std::vector<std::vector<int>> &nList,
                               bool target, bool usedDevice,
                               const tum::SpecializedCaps &start) {
  const FlatGraph g = flattenIndexList(nList);
  tum::CageCounts out;
  out.usedDevice = usedDevice;
  if (g.nAtoms <= 0) {
    return out;
  }
  constexpr long long kIntMax = std::numeric_limits<int>::max();
  Caps caps;
  caps.rowCap = std::clamp(start.rowCap, 1, tum::device::kRowCap);
  caps.ringCap = start.ringCap > 0 ? capProduct(start.ringCap, 1, kIntMax / 6)
                                   : capProduct(g.nAtoms, 16, kIntMax / 6);
  caps.pairCap = start.pairCap > 0 ? capProduct(start.pairCap, 1, kIntMax / 2)
                                   : capProduct(caps.ringCap, 8, kIntMax / 2);
  int retries = 0;
  for (;;) {
    out = tum::CageCounts{};
    out.usedDevice = usedDevice;
    const Need need = onePass(g, target, caps, out);
    if (need.rowMax > tum::device::kRowCap) {
      out = tum::hostCageCounts(nList);
      out.retries = retries;
      out.hostFallback = true;
      return out;
    }
    if (fits(need, caps)) {
      out.retries = retries;
      return out;
    }
    caps.ringCap = std::max(caps.ringCap, need.nRings);
    caps.rowCap = std::max(caps.rowCap, need.rowMax);
    caps.pairCap = std::max(caps.pairCap, need.nPairs);
    ++retries;
  }
}

} // namespace

bool tum::preferOffload() {
  const char *env = std::getenv("SEAMS_OFFLOAD");
  return env != nullptr && env[0] != '\0' && env[0] != '0';
}

tum::CageCounts
tum::hostCageCounts(const std::vector<std::vector<int>> &nList) {
  CageCounts out;
  out.usedDevice = false;
  if (nList.empty()) {
    return out;
  }
  const int nAtoms = static_cast<int>(nList.size());
  auto rings = primitive::ringNetwork(nList, 6);
  std::vector<std::vector<int>> six;
  six.reserve(rings.size());
  for (auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(std::move(r));
    }
  }
  const auto aff = ring::cageAffiliation(six, nList);
  out.nSix = static_cast<int>(six.size());
  out.atomSix.assign(static_cast<std::size_t>(nAtoms), 0);
  out.atomHc.assign(static_cast<std::size_t>(nAtoms), 0);
  out.atomDdc.assign(static_cast<std::size_t>(nAtoms), 0);
  for (std::size_t r = 0; r < six.size(); ++r) {
    if (aff.hc[r]) {
      ++out.nHcRings;
    }
    if (aff.ddc[r]) {
      ++out.nDdcRings;
    }
    for (const int a : six[r]) {
      if (a >= 0 && a < nAtoms) {
        ++out.atomSix[static_cast<std::size_t>(a)];
        if (aff.hc[r]) {
          out.atomHc[static_cast<std::size_t>(a)] = 1;
        }
        if (aff.ddc[r]) {
          out.atomDdc[static_cast<std::size_t>(a)] = 1;
        }
      }
    }
  }
  tallyAtoms(out);
  return out;
}

tum::CageCounts
tum::specializedCageCounts(const std::vector<std::vector<int>> &nList) {
  return specializedCageCounts(nList, SpecializedCaps{});
}

tum::CageCounts
tum::specializedCageCounts(const std::vector<std::vector<int>> &nList,
                           const SpecializedCaps &caps) {
#ifdef SEAMS_HAS_OFFLOAD
  const bool device = preferOffload() && omp_get_num_devices() > 0;
  if (device || caps.forceTarget) {
    return runSpecialized(nList, true, device, caps);
  }
#endif
  return runSpecialized(nList, false, false, caps);
}

tum::CageCounts tum::cageCounts(const std::vector<std::vector<int>> &nList) {
  if (preferOffload()) {
    return specializedCageCounts(nList);
  }
  return hostCageCounts(nList);
}
