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
#include <sphericart_ylm.hpp>
#include <steinhardt_device.hpp>

#include <algorithm>
#include <cstdlib>
#include <vector>

#ifdef SEAMS_HAS_MPI
#include <mpi.h>
#endif

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

constexpr int kParallelThreshold = 50000;

struct NeighbourCSR {
  std::vector<double> xyz;
  std::vector<int> offsets;
  std::vector<int> cols;
  double box[3];
  int nop = 0;
};

NeighbourCSR flatten(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &nList) {
  NeighbourCSR g;
  g.nop = yCloud.nop;
  g.box[0] = yCloud.box.size() > 0 ? yCloud.box[0] : 0.0;
  g.box[1] = yCloud.box.size() > 1 ? yCloud.box[1] : 0.0;
  g.box[2] = yCloud.box.size() > 2 ? yCloud.box[2] : 0.0;
  g.xyz.resize(static_cast<size_t>(g.nop) * 3);
  g.offsets.resize(static_cast<size_t>(g.nop) + 1, 0);

  for (int i = 0; i < g.nop; i++) {
    g.xyz[static_cast<size_t>(3 * i)] = yCloud.pts[i].x;
    g.xyz[static_cast<size_t>(3 * i + 1)] = yCloud.pts[i].y;
    g.xyz[static_cast<size_t>(3 * i + 2)] = yCloud.pts[i].z;
  }

  int nnz = 0;
  for (int i = 0; i < g.nop; i++) {
    if (static_cast<size_t>(i) >= nList.size() || nList[i].size() < 2) {
      g.offsets[static_cast<size_t>(i) + 1] = nnz;
      continue;
    }
    for (size_t j = 1; j < nList[i].size(); j++) {
      const auto it = yCloud.idIndexMap.find(nList[i][j]);
      if (it == yCloud.idIndexMap.end()) {
        continue;
      }
      if (it->second < 0 || it->second >= g.nop) {
        continue;
      }
      g.cols.push_back(it->second);
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
    const int iOff = 3 * i;
    for (int p = g.offsets[static_cast<size_t>(i)];
         p < g.offsets[static_cast<size_t>(i) + 1]; p++) {
      const int j = g.cols[static_cast<size_t>(p)];
      double dx = g.xyz[static_cast<size_t>(iOff)] - g.xyz[static_cast<size_t>(3 * j)];
      double dy = g.xyz[static_cast<size_t>(iOff + 1)] -
                  g.xyz[static_cast<size_t>(3 * j + 1)];
      double dz = g.xyz[static_cast<size_t>(iOff + 2)] -
                  g.xyz[static_cast<size_t>(3 * j + 2)];
      seams::steinhardt::minImage(dx, dy, dz, g.box[0], g.box[1], g.box[2]);
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
    for (int i = begin; i < end; i++) {
      seams::steinhardt::qlmOneAtom(i, orderL, g.xyz.data(), g.offsets.data(),
                                    g.cols.data(), g.box[0], g.box[1], g.box[2],
                                    qlm.data());
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
#ifdef SEAMS_HAS_SPHERICART
  if (seams::sphericart_ylm::available()) {
    runPass1Sphericart(g, orderL, begin, end, qlm);
    return;
  }
#endif
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = begin; i < end; i++) {
    seams::steinhardt::qlmOneAtom(i, orderL, g.xyz.data(), g.offsets.data(),
                                  g.cols.data(), g.box[0], g.box[1], g.box[2],
                                  qlm.data());
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
void runOffload(const NeighbourCSR &g, int orderL, std::vector<double> &qlm,
                std::vector<double> &ql, std::vector<double> &qlBar) {
  const int n = g.nop;
  const int nnz = static_cast<int>(g.cols.size());
  const double bx = g.box[0];
  const double by = g.box[1];
  const double bz = g.box[2];
  const double *xyz = g.xyz.data();
  const int *offsets = g.offsets.data();
  const int *cols = g.cols.data();
  double *qlmP = qlm.data();
  double *qlP = ql.data();
  double *qlBarP = qlBar.data();
  const int xyzN = 3 * n;
  const int offN = n + 1;
  const int qlmN = static_cast<int>(qlm.size());

#pragma omp target data map(to : xyz[0 : xyzN], offsets[0 : offN],                    \
                                cols[0 : nnz], bx, by, bz, orderL)                    \
    map(alloc : qlmP[0 : qlmN]) map(from : qlP[0 : n], qlBarP[0 : n])
  {
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; i++) {
      seams::steinhardt::qlmOneAtom(i, orderL, xyz, offsets, cols, bx, by, bz,
                                    qlmP);
    }
#pragma omp target teams distribute parallel for
    for (int i = 0; i < n; i++) {
      seams::steinhardt::qlOneAtom(i, orderL, qlmP, offsets, cols, qlP,
                                   qlBarP);
    }
  }
}
#endif

bool wantOffload() {
#ifdef SEAMS_HAS_OFFLOAD
  const char *env = std::getenv("SEAMS_OFFLOAD");
  if (env != nullptr && env[0] == '0') {
    return false;
  }
  return omp_get_num_devices() > 0;
#else
  return false;
#endif
}

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

  if (orderL != 3 && orderL != 4 && orderL != 6) {
    return result;
  }
  if (yCloud.nop <= 0) {
    return result;
  }

  const NeighbourCSR graph = flatten(yCloud, nList);
  const int nComp = 2 * orderL + 1;
  std::vector<double> qlm(static_cast<size_t>(graph.nop) * nComp * 2, 0.0);

#ifdef SEAMS_HAS_OFFLOAD
  if (wantOffload()) {
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
