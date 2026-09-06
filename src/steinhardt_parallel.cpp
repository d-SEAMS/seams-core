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
#include <cmath>
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

// Host-only l=12 Ylm. Associated Legendre with the same Condon-Shortley
// pairing as ylmAllTrig (Y_{l,m} = (-1)^m Y_{l,-m}^*). The device header
// keeps barRe[17] and no l=12 tables.
void ylm12Host(double sinT, double cosT, double cphi, double sphi,
               double *out) {
  constexpr int orderL = 12;
  constexpr int nComp = 25;
  constexpr double pi = 3.14159265358979323846;
  for (int k = 0; k < nComp; k++) {
    out[2 * k] = 0.0;
    out[2 * k + 1] = 0.0;
  }
  double pmm = 1.0;
  double fact = 1.0;
  double pL[13];
  for (int m = 0; m <= orderL; m++) {
    double plm;
    if (m == orderL) {
      plm = pmm;
    } else {
      double pmm1 = cosT * (2.0 * static_cast<double>(m) + 1.0) * pmm;
      if (m + 1 == orderL) {
        plm = pmm1;
      } else {
        double plm2 = pmm;
        double plm1 = pmm1;
        for (int l = m + 2; l <= orderL; l++) {
          const double pl =
              ((2.0 * static_cast<double>(l) - 1.0) * cosT * plm1 -
               (static_cast<double>(l + m) - 1.0) * plm2) /
              static_cast<double>(l - m);
          plm2 = plm1;
          plm1 = pl;
        }
        plm = plm1;
      }
    }
    pL[m] = plm;
    if (m < orderL) {
      pmm *= -fact * sinT;
      fact += 2.0;
    }
  }
  double pr[13];
  double piIm[13];
  pr[0] = 1.0;
  piIm[0] = 0.0;
  for (int k = 1; k <= orderL; k++) {
    pr[k] = pr[k - 1] * cphi - piIm[k - 1] * sphi;
    piIm[k] = pr[k - 1] * sphi + piIm[k - 1] * cphi;
  }
  for (int absM = 0; absM <= orderL; absM++) {
    double ratio = 1.0;
    for (int k = orderL - absM + 1; k <= orderL + absM; k++) {
      ratio /= static_cast<double>(k);
    }
    const double norm =
        std::sqrt((2.0 * orderL + 1.0) / (4.0 * pi) * ratio);
    const double amp = norm * pL[absM];
    const double nre = amp * pr[absM];
    const double nim = -amp * piIm[absM];
    const int ineg = orderL - absM;
    out[2 * ineg] = nre;
    out[2 * ineg + 1] = nim;
    const double sign = (absM % 2 == 0) ? 1.0 : -1.0;
    const int ipos = orderL + absM;
    out[2 * ipos] = sign * amp * pr[absM];
    out[2 * ipos + 1] = sign * amp * piIm[absM];
  }
}

void qlmAddBond12(double dx, double dy, double dz, double *qlmInterleaved,
                  int row, int &nUsed) {
  const double r2 = dx * dx + dy * dy + dz * dz;
  if (r2 == 0.0) {
    return;
  }
  const double r = std::sqrt(r2);
  const double invr = 1.0 / r;
  double cosT = dz * invr;
  if (cosT > 1.0) {
    cosT = 1.0;
  } else if (cosT < -1.0) {
    cosT = -1.0;
  }
  const double rho2 = dx * dx + dy * dy;
  double sinT = 0.0;
  double cphi = 1.0;
  double sphi = 0.0;
  if (rho2 != 0.0) {
    const double rho = std::sqrt(rho2);
    sinT = rho * invr;
    const double invrho = 1.0 / rho;
    cphi = dy * invrho;
    sphi = dx * invrho;
  }
  double ylm[50];
  ylm12Host(sinT, cosT, cphi, sphi, ylm);
  constexpr int nComp = 25;
  for (int m = 0; m < nComp; m++) {
    qlmInterleaved[2 * (row + m)] += ylm[2 * m];
    qlmInterleaved[2 * (row + m) + 1] += ylm[2 * m + 1];
  }
  nUsed++;
}

void runPass1Host12(const NeighbourCSR &g, int begin, int end,
                    std::vector<double> &qlm) {
  constexpr int nComp = 25;
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = begin; i < end; i++) {
    const int row = i * nComp;
    for (int m = 0; m < nComp; m++) {
      qlm[static_cast<size_t>(2 * (row + m))] = 0.0;
      qlm[static_cast<size_t>(2 * (row + m) + 1)] = 0.0;
    }
    const int j0 = g.offsets[static_cast<size_t>(i)];
    const int j1 = g.offsets[static_cast<size_t>(i) + 1];
    int nUsed = 0;
    for (int p = j0; p < j1; p++) {
      qlmAddBond12(g.dr[static_cast<size_t>(3 * p)],
                   g.dr[static_cast<size_t>(3 * p + 1)],
                   g.dr[static_cast<size_t>(3 * p + 2)], qlm.data(), row,
                   nUsed);
    }
    if (nUsed == 0) {
      continue;
    }
    const double inv = 1.0 / static_cast<double>(nUsed);
    for (int m = 0; m < nComp; m++) {
      qlm[static_cast<size_t>(2 * (row + m))] *= inv;
      qlm[static_cast<size_t>(2 * (row + m) + 1)] *= inv;
    }
  }
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
      runPass1Host12(g, begin, end, qlm);
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
  // l=12 is host-only. Prefer sphericart; otherwise associated Legendre.
  // Device arrays still cap at l=8.
  if (orderL == 12) {
#ifdef SEAMS_HAS_SPHERICART
    if (seams::sphericart_ylm::available()) {
      runPass1Sphericart(g, orderL, begin, end, qlm);
      return;
    }
#endif
    runPass1Host12(g, begin, end, qlm);
    return;
  }
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
    seams::steinhardt::qlmOneAtomDr(i, orderL, g.dr.data(), g.offsets.data(),
                                    g.cols.data(), qlm.data());
  }
}

// qlOneAtom keeps barRe[17] for the device l<=8 path. l=12 uses
// heap buffers on the host.
void runPass2Host(const NeighbourCSR &g, int orderL, int begin, int end,
                  const std::vector<double> &qlm, std::vector<double> &ql,
                  std::vector<double> &qlBar) {
  const int nComp = 2 * orderL + 1;
  constexpr double pi = 3.14159265358979323846;
  const double prefactor = 4.0 * pi / (2.0 * orderL + 1.0);
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = begin; i < end; i++) {
    const int row = i * nComp;
    double sumLocal = 0.0;
    for (int m = 0; m < nComp; m++) {
      const double re = qlm[static_cast<size_t>(2 * (row + m))];
      const double im = qlm[static_cast<size_t>(2 * (row + m) + 1)];
      sumLocal += re * re + im * im;
    }
    ql[static_cast<size_t>(i)] = std::sqrt(prefactor * sumLocal);
    const int j0 = g.offsets[static_cast<size_t>(i)];
    const int j1 = g.offsets[static_cast<size_t>(i) + 1];
    if (j0 == j1) {
      qlBar[static_cast<size_t>(i)] = ql[static_cast<size_t>(i)];
      continue;
    }
    std::vector<double> barRe(static_cast<size_t>(nComp), 0.0);
    std::vector<double> barIm(static_cast<size_t>(nComp), 0.0);
    for (int m = 0; m < nComp; m++) {
      barRe[static_cast<size_t>(m)] = qlm[static_cast<size_t>(2 * (row + m))];
      barIm[static_cast<size_t>(m)] = qlm[static_cast<size_t>(2 * (row + m) + 1)];
    }
    int nContrib = 1;
    for (int p = j0; p < j1; p++) {
      const int jatom = g.cols[static_cast<size_t>(p)];
      const int jRow = jatom * nComp;
      for (int m = 0; m < nComp; m++) {
        barRe[static_cast<size_t>(m)] +=
            qlm[static_cast<size_t>(2 * (jRow + m))];
        barIm[static_cast<size_t>(m)] +=
            qlm[static_cast<size_t>(2 * (jRow + m) + 1)];
      }
      nContrib++;
    }
    const double inv = 1.0 / static_cast<double>(nContrib);
    double sumBar = 0.0;
    for (int m = 0; m < nComp; m++) {
      const double re = barRe[static_cast<size_t>(m)] * inv;
      const double im = barIm[static_cast<size_t>(m)] * inv;
      sumBar += re * re + im * im;
    }
    qlBar[static_cast<size_t>(i)] = std::sqrt(prefactor * sumBar);
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
  // Device Ylm has no l=12; force the host path (sphericart or Legendre).
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
  const bool hostPass2 = orderL == 12;

#ifdef SEAMS_HAS_MPI
  if (initialized && nranks > 1) {
    // Lechner-Dellago qlBar averages q_lm over neighbours that may live
    // on another rank, so the first pass is gathered before the second.
    allgathervDoubles(qlm, graph.nop, 2 * nComp, rank, nranks);
  }
#endif
  if (hostPass2) {
    runPass2Host(graph, orderL, begin, end, qlm, result.ql, result.qlBar);
  } else {
    runPass2(graph, orderL, begin, end, qlm, result.ql, result.qlBar);
  }
#ifdef SEAMS_HAS_MPI
  if (initialized && nranks > 1) {
    allgathervDoubles(result.ql, graph.nop, 1, rank, nranks);
    allgathervDoubles(result.qlBar, graph.nop, 1, rank, nranks);
  }
#endif

  return result;
}

} // namespace chill
