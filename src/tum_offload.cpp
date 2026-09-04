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

tum::CageCounts runSpecialized(FlatGraph g, bool onDevice) {
  tum::CageCounts out;
  out.usedDevice = onDevice;
  if (g.nAtoms <= 0) {
    return out;
  }
  const int nAtoms = g.nAtoms;
  const int kMax = g.kMax;
  const int maxPer = 16;
  const int maxRings = nAtoms * maxPer;
  const int maxPairs = maxRings * 8;
  std::vector<int> nRings(1, 0);
  std::vector<int> dropped(1, 0);
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

  int *deg = g.deg.data();
  int *cols = g.cols.data();
  int *nRingsP = nRings.data();
  int *droppedP = dropped.data();
  int *ringAtomsP = ringAtoms.data();
  int *throughCountP = throughCount.data();
  int *throughP = through.data();
  int *hcP = hc.data();
  int *ddcP = ddc.data();
  int *nPairsP = nPairs.data();
  int *pairsP = pairs.data();
  int *atomHcP = out.atomHc.data();
  int *atomDdcP = out.atomDdc.data();
  const int degN = nAtoms;
  const int colN = nAtoms * kMax;
  const int ringN = maxRings * 6;
  const int thruN = nAtoms * maxPer;
  const int pairN = maxPairs * 2;

#ifdef SEAMS_HAS_OFFLOAD
  if (onDevice) {
#pragma omp target data map(to : deg[0 : degN], cols[0 : colN], nAtoms, kMax,  \
                                maxRings, maxPer, maxPairs)                    \
    map(tofrom : nRingsP[0 : 1], droppedP[0 : 1], nPairsP[0 : 1],              \
            ringAtomsP[0 : ringN], throughCountP[0 : nAtoms],                  \
            throughP[0 : thruN], hcP[0 : maxRings], ddcP[0 : maxRings],        \
            pairsP[0 : pairN], atomHcP[0 : nAtoms], atomDdcP[0 : nAtoms])
    {
#pragma omp target teams distribute parallel for
      for (int i = 0; i < nAtoms; ++i) {
        tum::device::enumSixFrom(i, deg, cols, nAtoms, kMax, maxRings, nRingsP,
                                 ringAtomsP, droppedP);
      }
      // Kernels no-op when the index is past nRings/nPairs; do not read
      // those counters on the host inside the target data region.
#pragma omp target teams distribute parallel for
      for (int r = 0; r < maxRings; ++r) {
        tum::device::invertOneRing(r, nRingsP, ringAtomsP, nAtoms, maxPer,
                                   throughCountP, throughP);
      }
#pragma omp target teams distribute parallel for
      for (int a = 0; a < nAtoms; ++a) {
        int n = throughCountP[a];
        if (n > maxPer) {
          n = maxPer;
        }
        throughCountP[a] = n;
        tum::device::sortThroughRow(throughP + a * maxPer, n);
      }
#pragma omp target teams distribute parallel for
      for (int r = 0; r < maxRings; ++r) {
        tum::device::emitBasalFrom(r, nRingsP, ringAtomsP, deg, cols,
                                   throughCountP, throughP, nAtoms, kMax,
                                   maxPer, maxPairs, nPairsP, pairsP);
      }
#pragma omp target teams distribute parallel for
      for (int p = 0; p < maxPairs; ++p) {
        tum::device::applyHcPair(p, nPairsP, pairsP, ringAtomsP, throughCountP,
                                 throughP, nAtoms, maxPer, hcP);
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
    out.nSix = nRings[0] < maxRings ? nRings[0] : maxRings;
    out.ringsDropped = dropped[0];
    out.nHcRings = 0;
    out.nDdcRings = 0;
    for (int r = 0; r < out.nSix; ++r) {
      out.nHcRings += hc[static_cast<std::size_t>(r)];
      out.nDdcRings += ddc[static_cast<std::size_t>(r)];
    }
    tallyAtoms(out);
    return out;
  }
#else
  (void)onDevice;
  (void)degN;
  (void)colN;
  (void)ringN;
  (void)thruN;
  (void)pairN;
#endif

  for (int i = 0; i < nAtoms; ++i) {
    tum::device::enumSixFrom(i, deg, cols, nAtoms, kMax, maxRings, nRingsP,
                             ringAtomsP, droppedP);
  }
  const int found = nRings[0] < maxRings ? nRings[0] : maxRings;
  for (int r = 0; r < found; ++r) {
    tum::device::invertOneRing(r, nRingsP, ringAtomsP, nAtoms, maxPer,
                               throughCountP, throughP);
  }
  for (int a = 0; a < nAtoms; ++a) {
    int n = throughCount[static_cast<std::size_t>(a)];
    if (n > maxPer) {
      n = maxPer;
    }
    throughCount[static_cast<std::size_t>(a)] = n;
    tum::device::sortThroughRow(throughP + a * maxPer, n);
  }
  for (int r = 0; r < found; ++r) {
    tum::device::emitBasalFrom(r, nRingsP, ringAtomsP, deg, cols, throughCountP,
                               throughP, nAtoms, kMax, maxPer, maxPairs,
                               nPairsP, pairsP);
  }
  const int np = nPairs[0] < maxPairs ? nPairs[0] : maxPairs;
  for (int p = 0; p < np; ++p) {
    tum::device::applyHcPair(p, nPairsP, pairsP, ringAtomsP, throughCountP,
                             throughP, nAtoms, maxPer, hcP);
  }
  for (int r = 0; r < found; ++r) {
    tum::device::ddcFrom(r, nRingsP, ringAtomsP, throughCountP, throughP, hcP,
                         nAtoms, maxPer, ddcP);
  }
  for (int r = 0; r < found; ++r) {
    tum::device::atomIceFrom(r, nRingsP, ringAtomsP, hcP, ddcP, nAtoms, atomHcP,
                             atomDdcP);
  }
  out.nSix = found;
  out.ringsDropped = dropped[0];
  out.nHcRings = 0;
  out.nDdcRings = 0;
  for (int r = 0; r < found; ++r) {
    out.nHcRings += hc[static_cast<std::size_t>(r)];
    out.nDdcRings += ddc[static_cast<std::size_t>(r)];
  }
  tallyAtoms(out);
  return out;
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
  const FlatGraph g = flattenIndexList(nList);
#ifdef SEAMS_HAS_OFFLOAD
  if (preferOffload() && omp_get_num_devices() > 0) {
    return runSpecialized(g, true);
  }
#endif
  return runSpecialized(g, false);
}

tum::CageCounts tum::cageCounts(const std::vector<std::vector<int>> &nList) {
  if (preferOffload()) {
    return specializedCageCounts(nList);
  }
  return hostCageCounts(nList);
}
