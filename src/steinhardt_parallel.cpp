//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <bop.hpp>
#include <generic.hpp>
#include <sphericart_ylm.hpp>
#include <steinhardt_device.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <vector>

#ifdef SEAMS_HAS_MPI
#include <mpi.h>
#endif

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

constexpr int kParallelThreshold = 256;

struct NeighbourCSR {
  std::vector<double> dr;
  std::vector<int> offsets;
  std::vector<int> cols;
  int nop = 0;
};

NeighbourCSR flatten(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &nList) {
  NeighbourCSR g;
  g.nop = yCloud.nop;
  g.offsets.resize(static_cast<size_t>(g.nop) + 1, 0);

  int maxId = 0;
  for (int i = 0; i < g.nop; i++) {
    if (yCloud.pts[static_cast<size_t>(i)].atomID > maxId) {
      maxId = yCloud.pts[static_cast<size_t>(i)].atomID;
    }
  }
  const bool dense = maxId > 0 && maxId < g.nop * 8 + 8;
  std::vector<int> idToIdx;
  if (dense) {
    idToIdx.assign(static_cast<size_t>(maxId) + 1, -1);
    for (int i = 0; i < g.nop; i++) {
      idToIdx[static_cast<size_t>(yCloud.pts[static_cast<size_t>(i)].atomID)] = i;
    }
  }

  int nnzGuess = 0;
  for (int i = 0; i < g.nop && static_cast<size_t>(i) < nList.size(); i++) {
    if (nList[static_cast<size_t>(i)].size() > 1) {
      nnzGuess += static_cast<int>(nList[static_cast<size_t>(i)].size()) - 1;
    }
  }
  g.cols.reserve(static_cast<size_t>(std::max(nnzGuess, 0)));
  g.dr.reserve(static_cast<size_t>(std::max(nnzGuess, 0)) * 3);

  int nnz = 0;
  for (int i = 0; i < g.nop; i++) {
    if (static_cast<size_t>(i) >= nList.size() ||
        nList[static_cast<size_t>(i)].size() < 2) {
      g.offsets[static_cast<size_t>(i) + 1] = nnz;
      continue;
    }
    for (size_t j = 1; j < nList[static_cast<size_t>(i)].size(); j++) {
      const int atomId = nList[static_cast<size_t>(i)][j];
      int jidx = -1;
      if (dense && atomId >= 0 && atomId <= maxId) {
        jidx = idToIdx[static_cast<size_t>(atomId)];
      } else {
        const auto it = yCloud.idIndexMap.find(atomId);
        if (it != yCloud.idIndexMap.end()) {
          jidx = it->second;
        }
      }
      if (jidx < 0 || jidx >= g.nop) {
        continue;
      }
      g.cols.push_back(jidx);
      const auto d = gen::relDist(yCloud, i, jidx);
      g.dr.push_back(d[0]);
      g.dr.push_back(d[1]);
      g.dr.push_back(d[2]);
      nnz++;
    }
    g.offsets[static_cast<size_t>(i) + 1] = nnz;
  }
  return g;
}

void atomRange(int n, int rank, int nranks, int &begin, int &end) {
  const int base = n / nranks;
  const int rem = n % nranks;
  begin = rank * base + std::min(rank, rem);
  end = begin + base + (rank < rem ? 1 : 0);
}

#ifdef SEAMS_HAS_SPHERICART
void runPass1Sphericart(const NeighbourCSR &g, int orderL, int begin, int end,
                        std::vector<double> &qlm) {
  const int nComp = 2 * orderL + 1;
  int nBonds = 0;
  for (int i = begin; i < end; i++) {
    nBonds += g.offsets[static_cast<size_t>(i) + 1] - g.offsets[static_cast<size_t>(i)];
  }
  std::vector<double> cart(static_cast<size_t>(std::max(nBonds, 0)) * 3, 0.0);
  std::vector<int> owner;
  owner.reserve(static_cast<size_t>(std::max(nBonds, 0)));
  int b = 0;
  for (int i = begin; i < end; i++) {
    for (int p = g.offsets[static_cast<size_t>(i)];
         p < g.offsets[static_cast<size_t>(i) + 1]; p++) {
      const double dx = g.dr[static_cast<size_t>(3 * p)];
      const double dy = g.dr[static_cast<size_t>(3 * p + 1)];
      const double dz = g.dr[static_cast<size_t>(3 * p + 2)];
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 == 0.0) {
        continue;
      }
      const double inv = 1.0 / std::sqrt(r2);
      cart[static_cast<size_t>(3 * b)] = dx * inv;
      cart[static_cast<size_t>(3 * b + 1)] = dy * inv;
      cart[static_cast<size_t>(3 * b + 2)] = dz * inv;
      owner.push_back(i);
      b++;
    }
  }
  std::vector<double> ylm(static_cast<size_t>(b) * nComp * 2, 0.0);
  if (seams::sphericart_ylm::ylmCartesian(orderL, cart.data(), b, ylm.data()) !=
      0) {
    if (orderL == 12) {
      return;
    }
    for (int i = begin; i < end; i++) {
      seams::steinhardt::qlmOneAtomDr(i, orderL, g.dr.data(), g.offsets.data(),
                                      g.cols.data(), qlm.data());
    }
    return;
  }
  for (int k = 0; k < b; k++) {
    const int i = owner[static_cast<size_t>(k)];
    const size_t row = static_cast<size_t>(i) * nComp;
    const double *src = ylm.data() + static_cast<size_t>(k) * nComp * 2;
    for (int m = 0; m < nComp; m++) {
      qlm[2 * (row + static_cast<size_t>(m))] += src[2 * m];
      qlm[2 * (row + static_cast<size_t>(m)) + 1] += src[2 * m + 1];
    }
  }
  for (int i = begin; i < end; i++) {
    int nUsed = 0;
    for (int k = 0; k < b; k++) {
      if (owner[static_cast<size_t>(k)] == i) {
        nUsed++;
      }
    }
    if (nUsed == 0) {
      continue;
    }
    const double inv = 1.0 / static_cast<double>(nUsed);
    const size_t row = static_cast<size_t>(i) * nComp;
    for (int m = 0; m < nComp; m++) {
      qlm[2 * (row + static_cast<size_t>(m))] *= inv;
      qlm[2 * (row + static_cast<size_t>(m)) + 1] *= inv;
    }
  }
}
#endif

void runPass1(const NeighbourCSR &g, int orderL, int begin, int end,
              std::vector<double> &qlm) {
  // l=12 has no closed-form device/host Ylm (arrays cap at l=8).
#ifdef SEAMS_HAS_SPHERICART
  if (orderL == 12) {
    if (seams::sphericart_ylm::available()) {
      runPass1Sphericart(g, orderL, begin, end, qlm);
    }
    return;
  }
  if (seams::sphericart_ylm::available()) {
    runPass1Sphericart(g, orderL, begin, end, qlm);
    return;
  }
#else
  if (orderL == 12) {
    return;
  }
#endif
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = begin; i < end; i++) {
    seams::steinhardt::qlmOneAtomDr(i, orderL, g.dr.data(), g.offsets.data(),
                                    g.cols.data(), qlm.data());
  }
}

void runPass2(const NeighbourCSR &g, int orderL, int begin, int end,
              const std::vector<double> &qlm, std::vector<double> &ql,
              std::vector<double> &qlBar) {
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = begin; i < end; i++) {
    seams::steinhardt::qlOneAtom(i, orderL, qlm.data(), g.offsets.data(),
                                 g.cols.data(), ql.data(), qlBar.data());
  }
}

#ifdef SEAMS_HAS_OFFLOAD
bool wantOffload() {
  const char *env = std::getenv("SEAMS_OFFLOAD");
  if (env != nullptr && env[0] == '0') {
    return false;
  }
  return omp_get_num_devices() > 0;
}

// One device allocation per process. map() on every call was a
// managed malloc (~20 ms) around a ~0.1 ms Ylm kernel.
struct DeviceScratch {
  int device = -1;
  int nCap = 0;
  int nnzCap = 0;
  int qlmCap = 0;
  double *dr = nullptr;
  int *offsets = nullptr;
  int *cols = nullptr;
  double *qlm = nullptr;

  void release() {
    if (device < 0) {
      return;
    }
    if (dr != nullptr) {
      omp_target_free(dr, device);
    }
    if (offsets != nullptr) {
      omp_target_free(offsets, device);
    }
    if (cols != nullptr) {
      omp_target_free(cols, device);
    }
    if (qlm != nullptr) {
      omp_target_free(qlm, device);
    }
    dr = nullptr;
    offsets = nullptr;
    cols = nullptr;
    qlm = nullptr;
    nCap = 0;
    nnzCap = 0;
    qlmCap = 0;
    device = -1;
  }

  bool ensure(int n, int nnz, int qlmN) {
    const int dev = omp_get_default_device();
    if (device == dev && nCap >= n + 1 && nnzCap >= nnz && qlmCap >= qlmN &&
        dr != nullptr && offsets != nullptr && cols != nullptr &&
        qlm != nullptr) {
      return true;
    }
    release();
    device = dev;
    const int nnzUse = std::max(nnz, 1);
    const int qlmUse = std::max(qlmN, 1);
    dr = static_cast<double *>(
        omp_target_alloc(sizeof(double) * 3 * static_cast<size_t>(nnzUse), dev));
    offsets = static_cast<int *>(
        omp_target_alloc(sizeof(int) * static_cast<size_t>(n + 1), dev));
    cols = static_cast<int *>(
        omp_target_alloc(sizeof(int) * static_cast<size_t>(nnzUse), dev));
    qlm = static_cast<double *>(
        omp_target_alloc(sizeof(double) * static_cast<size_t>(qlmUse), dev));
    if (dr == nullptr || offsets == nullptr || cols == nullptr ||
        qlm == nullptr) {
      release();
      return false;
    }
    nCap = n + 1;
    nnzCap = nnzUse;
    qlmCap = qlmUse;
    return true;
  }
};

DeviceScratch gOffloadScratch;

void runOffloadMapped(const NeighbourCSR &g, int orderL, std::vector<double> &qlm,
                      std::vector<double> &ql, std::vector<double> &qlBar) {
  const int n = g.nop;
  const int nnz = static_cast<int>(g.cols.size());
  const double *dr = g.dr.data();
  const int *offsets = g.offsets.data();
  const int *cols = g.cols.data();
  double *qlmP = qlm.data();
  const int drN = 3 * nnz;
  const int offN = n + 1;
  const int qlmN = static_cast<int>(qlm.size());
#pragma omp target data map(to : dr[0 : drN], offsets[0 : offN],                      \
                                cols[0 : nnz], orderL)                                \
    map(from : qlmP[0 : qlmN])
  {
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; i++) {
      seams::steinhardt::qlmOneAtomDr(i, orderL, dr, offsets, cols, qlmP);
    }
  }
  runPass2(g, orderL, 0, n, qlm, ql, qlBar);
}

void runOffload(const NeighbourCSR &g, int orderL, std::vector<double> &qlm,
                std::vector<double> &ql, std::vector<double> &qlBar) {
  const int n = g.nop;
  const int nnz = static_cast<int>(g.cols.size());
  const int qlmN = static_cast<int>(qlm.size());
  if (n <= 0 || nnz <= 0 || qlmN <= 0 ||
      !gOffloadScratch.ensure(n, nnz, qlmN)) {
    runOffloadMapped(g, orderL, qlm, ql, qlBar);
    return;
  }
  const int hostDev = omp_get_initial_device();
  const int dev = gOffloadScratch.device;
  if (omp_target_memcpy(gOffloadScratch.dr, g.dr.data(),
                        sizeof(double) * 3 * static_cast<size_t>(nnz), 0, 0, dev,
                        hostDev) != 0 ||
      omp_target_memcpy(gOffloadScratch.offsets, g.offsets.data(),
                        sizeof(int) * static_cast<size_t>(n + 1), 0, 0, dev,
                        hostDev) != 0 ||
      omp_target_memcpy(gOffloadScratch.cols, g.cols.data(),
                        sizeof(int) * static_cast<size_t>(nnz), 0, 0, dev,
                        hostDev) != 0) {
    runOffloadMapped(g, orderL, qlm, ql, qlBar);
    return;
  }
  double *dDr = gOffloadScratch.dr;
  int *dOff = gOffloadScratch.offsets;
  int *dCols = gOffloadScratch.cols;
  double *dQlm = gOffloadScratch.qlm;
#pragma omp target teams distribute parallel for is_device_ptr(dDr, dOff, dCols, \
                                                               dQlm)
  for (int i = 0; i < n; i++) {
    seams::steinhardt::qlmOneAtomDr(i, orderL, dDr, dOff, dCols, dQlm);
  }
  if (omp_target_memcpy(qlm.data(), dQlm,
                        sizeof(double) * static_cast<size_t>(qlmN), 0, 0,
                        hostDev, dev) != 0) {
    runOffloadMapped(g, orderL, qlm, ql, qlBar);
    return;
  }
  runPass2(g, orderL, 0, n, qlm, ql, qlBar);
}
#endif

#ifdef SEAMS_HAS_MPI
void allgathervDoubles(std::vector<double> &buf, int nAtoms, int perAtom,
                       int rank, int nranks) {
  std::vector<int> counts(nranks);
  std::vector<int> displs(nranks);
  int disp = 0;
  for (int r = 0; r < nranks; r++) {
    int b = 0;
    int e = 0;
    atomRange(nAtoms, r, nranks, b, e);
    counts[r] = (e - b) * perAtom;
    displs[r] = disp;
    disp += counts[r];
  }
  std::vector<double> send(buf.begin() + displs[rank],
                           buf.begin() + displs[rank] + counts[rank]);
  MPI_Allgatherv(send.data(), counts[rank], MPI_DOUBLE, buf.data(),
                 counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
}
#endif

} // namespace

namespace chill {

SteinhardtQl steinhardtQl(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                          const std::vector<std::vector<int>> &nList, int orderL) {
  SteinhardtQl result;
  result.ql.assign(yCloud.nop, 0.0);
  result.qlBar.assign(yCloud.nop, 0.0);

  if (orderL != 3 && orderL != 4 && orderL != 6 && orderL != 8 &&
      orderL != 12) {
    return result;
  }
  if (yCloud.nop <= 0) {
    return result;
  }

  const NeighbourCSR graph = flatten(yCloud, nList);
  const int nComp = 2 * orderL + 1;
  std::vector<double> qlm(static_cast<size_t>(graph.nop) * nComp * 2, 0.0);

#ifdef SEAMS_HAS_OFFLOAD
  // Device Ylm has no l=12; force the host sphericart path.
  if (orderL != 12 && wantOffload()) {
    runOffload(graph, orderL, qlm, result.ql, result.qlBar);
    return result;
  }
#endif

  int rank = 0;
  int nranks = 1;
#ifdef SEAMS_HAS_MPI
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  }
#endif

  int begin = 0;
  int end = graph.nop;
  atomRange(graph.nop, rank, nranks, begin, end);
  runPass1(graph, orderL, begin, end, qlm);

#ifdef SEAMS_HAS_MPI
  if (initialized && nranks > 1) {
    // Lechner-Dellago qlBar averages q_lm over neighbours that may live
    // on another rank, so the first pass is gathered before the second.
    allgathervDoubles(qlm, graph.nop, 2 * nComp, rank, nranks);
  }
#endif
  runPass2(graph, orderL, begin, end, qlm, result.ql, result.qlBar);
#ifdef SEAMS_HAS_MPI
  if (initialized && nranks > 1) {
    allgathervDoubles(result.ql, graph.nop, 1, rank, nranks);
    allgathervDoubles(result.qlBar, graph.nop, 1, rank, nranks);
  }
#endif

  return result;
}

} // namespace chill
