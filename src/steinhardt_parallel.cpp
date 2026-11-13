//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) 2018--present d-SEAMS core team
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <chrono>

#include <bop.hpp>
#include <domain.hpp>
#include <generic.hpp>
#include <neighbours.hpp>
#include <ozaki.hpp>
#include <sphericart_ylm.hpp>
#include <steinhardt_device.hpp>
#include <thread_team.hpp>
#include <wigner_wl.hpp>
#include <ylm_batch.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <numeric>
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
  // Frame IDs stay large on a rank that holds only its block.
  const bool dense = maxId > 0 && maxId < 4000000;
  std::vector<int> idToIdx;
  if (dense) {
    idToIdx.assign(static_cast<size_t>(maxId) + 1, -1);
    for (int i = 0; i < g.nop; i++) {
      idToIdx[static_cast<size_t>(yCloud.pts[static_cast<size_t>(i)].atomID)] = i;
    }
  }

  // The cloud index of a neighbour's atom ID; negative when it has none
  const auto indexOf = [&](int atomId) {
    int jidx = -1;
    if (dense && atomId >= 0 && atomId <= maxId) {
      jidx = idToIdx[static_cast<size_t>(atomId)];
    } else {
      const auto it = yCloud.idIndexMap.find(atomId);
      if (it != yCloud.idIndexMap.end()) {
        jidx = it->second;
      }
    }
    return jidx < g.nop ? jidx : -1;
  };
  const int nRows = std::min(g.nop, static_cast<int>(nList.size()));

  // A row keeps the same neighbours in both passes and writes only its own
  // run of cols and dr, so the rows are independent.
  // nvc++ -mp=gpu SIGSEGVs host parallel-fors in libnvomp.
  // The thread test is outside the pragma. An `if` clause on the pragma
  // is read again by each worker, which races with the master's stack.
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  const bool useThreads = g.nop >= kParallelThreshold;
#else
  const bool useThreads = false;
#endif
  auto countRow = [&](int i) {
    int kept = 0;
    std::vector<int> seen;
    seen.reserve(8);
    for (size_t j = 1; j < nList[static_cast<size_t>(i)].size(); j++) {
      const int jidx = indexOf(nList[static_cast<size_t>(i)][j]);
      if (jidx < 0) {
        continue;
      }
      if (std::find(seen.begin(), seen.end(), jidx) != seen.end()) {
        continue;
      }
      seen.push_back(jidx);
      kept++;
    }
    g.offsets[static_cast<size_t>(i) + 1] = kept;
  };
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  if (useThreads) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < nRows; i++) {
      countRow(i);
    }
  } else
#endif
  for (int i = 0; i < nRows; i++) {
    countRow(i);
  }
  std::partial_sum(g.offsets.begin(), g.offsets.end(), g.offsets.begin());
  const int nnz = g.offsets[static_cast<size_t>(g.nop)];
  g.cols.resize(static_cast<size_t>(nnz));
  g.dr.resize(static_cast<size_t>(nnz) * 3);

  const gen::FracBox frame = gen::makeFracBox(yCloud);
  auto fillRow = [&](int i) {
    int k = g.offsets[static_cast<size_t>(i)];
    std::vector<int> seen;
    seen.reserve(8);
    for (size_t j = 1; j < nList[static_cast<size_t>(i)].size(); j++) {
      const int jidx = indexOf(nList[static_cast<size_t>(i)][j]);
      if (jidx < 0) {
        continue;
      }
      if (std::find(seen.begin(), seen.end(), jidx) != seen.end()) {
        continue;
      }
      seen.push_back(jidx);
      const auto &pj = yCloud.pts[static_cast<size_t>(jidx)];
      const auto d = gen::relDistFromPoint(frame, yCloud, i, pj.x, pj.y, pj.z);
      g.cols[static_cast<size_t>(k)] = jidx;
      g.dr[3 * static_cast<size_t>(k)] = d[0];
      g.dr[3 * static_cast<size_t>(k) + 1] = d[1];
      g.dr[3 * static_cast<size_t>(k) + 2] = d[2];
      k++;
    }
  };
#if defined(SEAMS_HAS_OPENMP) && !defined(SEAMS_HAS_OFFLOAD)
  if (useThreads) {
    const int rowCount = nRows;
#pragma omp parallel for schedule(static) firstprivate(rowCount)
    for (int i = 0; i < rowCount; i++) {
      fillRow(i);
    }
  } else
#endif
  for (int i = 0; i < nRows; i++) {
    fillRow(i);
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
  const bool pairwise =
      seams::product::reduction() != seams::product::Reduction::Left;
  if (!pairwise) {
    for (int k = 0; k < b; k++) {
      const int i = owner[static_cast<size_t>(k)];
      const size_t row = static_cast<size_t>(i) * nComp;
      const double *src = ylm.data() + static_cast<size_t>(k) * nComp * 2;
      for (int m = 0; m < nComp; m++) {
        qlm[2 * (row + static_cast<size_t>(m))] += src[2 * m];
        qlm[2 * (row + static_cast<size_t>(m)) + 1] += src[2 * m + 1];
      }
    }
  } else {
    int k = 0;
    while (k < b) {
      const int i = owner[static_cast<size_t>(k)];
      int k1 = k + 1;
      while (k1 < b && owner[static_cast<size_t>(k1)] == i) {
        ++k1;
      }
      const int cnt = k1 - k;
      std::vector<double> comp(static_cast<size_t>(cnt));
      const size_t row = static_cast<size_t>(i) * nComp;
      for (int m = 0; m < nComp; m++) {
        for (int t = 0; t < cnt; t++) {
          const double *src =
              ylm.data() + static_cast<size_t>(k + t) * nComp * 2;
          comp[static_cast<size_t>(t)] = src[2 * m];
        }
        qlm[2 * (row + static_cast<size_t>(m))] = seams::product::reduce(
            comp.data(), cnt, seams::product::reduction());
        for (int t = 0; t < cnt; t++) {
          const double *src =
              ylm.data() + static_cast<size_t>(k + t) * nComp * 2;
          comp[static_cast<size_t>(t)] = src[2 * m + 1];
        }
        qlm[2 * (row + static_cast<size_t>(m)) + 1] = seams::product::reduce(
            comp.data(), cnt, seams::product::reduction());
      }
      k = k1;
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

void qlmOneAtomDrPairwise(int iatom, int orderL, const double *dr,
                          const int *offsets, double *qlmInterleaved) {
  const int nComp = 2 * orderL + 1;
  const int row = iatom * nComp;
  for (int m = 0; m < nComp; m++) {
    qlmInterleaved[2 * (row + m)] = 0.0;
    qlmInterleaved[2 * (row + m) + 1] = 0.0;
  }
  if (orderL > 8) {
    return;
  }
  const int j0 = offsets[iatom];
  const int j1 = offsets[iatom + 1];
  std::vector<double> bonds;
  bonds.reserve(static_cast<std::size_t>(std::max(0, j1 - j0)) * static_cast<std::size_t>(nComp) * 2);
  int nUsed = 0;
  for (int p = j0; p < j1; p++) {
    double tmp[34] = {};
    int one = 0;
    seams::steinhardt::qlmAddBond(orderL, dr[3 * p], dr[3 * p + 1], dr[3 * p + 2],
                                  tmp, 0, nComp, one);
    if (one == 0) {
      continue;
    }
    bonds.insert(bonds.end(), tmp, tmp + nComp * 2);
    ++nUsed;
  }
  if (nUsed == 0) {
    return;
  }
  std::vector<double> comp(static_cast<std::size_t>(nUsed));
  for (int m = 0; m < nComp; m++) {
    for (int t = 0; t < nUsed; t++) {
      comp[static_cast<std::size_t>(t)] = bonds[static_cast<std::size_t>(t) * static_cast<std::size_t>(nComp) * 2 + static_cast<std::size_t>(2 * m)];
    }
    qlmInterleaved[2 * (row + m)] =
        seams::product::reduce(comp.data(), nUsed, seams::product::reduction());
    for (int t = 0; t < nUsed; t++) {
      comp[static_cast<std::size_t>(t)] =
          bonds[static_cast<std::size_t>(t) * static_cast<std::size_t>(nComp) * 2 + static_cast<std::size_t>(2 * m) + 1];
    }
    qlmInterleaved[2 * (row + m) + 1] =
        seams::product::reduce(comp.data(), nUsed, seams::product::reduction());
  }
  const double inv = 1.0 / static_cast<double>(nUsed);
  for (int m = 0; m < nComp; m++) {
    qlmInterleaved[2 * (row + m)] *= inv;
    qlmInterleaved[2 * (row + m) + 1] *= inv;
  }
}

// Sum one atom's bonds. Groups of eight go through the Cartesian batch;
// the sum is still in neighbour-list order. l is 3, 4 or 6.
void qlmOneAtomBatch(int iatom, int orderL, const double *dr, const int *offsets,
                     double *qlm) {
  const int nComp = 2 * orderL + 1;
  const int row = iatom * nComp;
  for (int m = 0; m < nComp; m++) {
    qlm[2 * (row + m)] = 0.0;
    qlm[2 * (row + m) + 1] = 0.0;
  }
  const int j0 = offsets[iatom];
  const int j1 = offsets[iatom + 1];
  double xs[8];
  double ys[8];
  double zs[8];
  // l <= 6, so one batch holds at most 8 * 13 components.
  double re[8 * 13];
  double im[8 * 13];
  int filled = 0;
  int nUsed = 0;
  const auto flush = [&]() {
    if (filled == 0) {
      return;
    }
    seams::ylm::ylmBatch(orderL, xs, ys, zs, filled, re, im);
    for (int b = 0; b < filled; b++) {
      const double *yr = re + static_cast<std::size_t>(b) * static_cast<std::size_t>(nComp);
      const double *yi = im + static_cast<std::size_t>(b) * static_cast<std::size_t>(nComp);
      for (int m = 0; m < nComp; m++) {
        qlm[2 * (row + m)] += yr[m];
        qlm[2 * (row + m) + 1] += yi[m];
      }
    }
    filled = 0;
  };
  for (int p = j0; p < j1; p++) {
    const double dx = dr[3 * p];
    const double dy = dr[3 * p + 1];
    const double dz = dr[3 * p + 2];
    if (dx * dx + dy * dy + dz * dz == 0.0) {
      continue;
    }
    xs[filled] = dx;
    ys[filled] = dy;
    zs[filled] = dz;
    ++filled;
    ++nUsed;
    if (filled == 8) {
      flush();
    }
  }
  flush();
  if (nUsed == 0) {
    return;
  }
  const double inv = 1.0 / static_cast<double>(nUsed);
  for (int m = 0; m < nComp; m++) {
    qlm[2 * (row + m)] *= inv;
    qlm[2 * (row + m) + 1] *= inv;
  }
}

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
  const bool pairwise =
      seams::product::reduction() != seams::product::Reduction::Left && orderL <= 8;
  // Degrees 4 and 6 sum eight bonds through the Cartesian harmonic.
  // Each atom still adds its own bonds in list order. Degree 3 stays
  // on the one-bond kernel: the batch normalises by r^3, which is a
  // square root and a divide, and on a four-neighbour shell that costs
  // more than the seven components it writes.
  const bool batchSolid = (orderL == 4 || orderL == 6) && !pairwise;
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold && omp_get_max_threads() > 1;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = begin; i < end; i++) {
    if (pairwise) {
      qlmOneAtomDrPairwise(i, orderL, g.dr.data(), g.offsets.data(), qlm.data());
    } else if (batchSolid) {
      qlmOneAtomBatch(i, orderL, g.dr.data(), g.offsets.data(), qlm.data());
    } else {
      seams::steinhardt::qlmOneAtomDr(i, orderL, g.dr.data(), g.offsets.data(),
                                      g.cols.data(), qlm.data());
    }
  }
}

void runPass1Mask(const NeighbourCSR &g, int orderL, std::vector<double> &qlm,
                  const std::vector<char> *which) {
  const bool pairwise =
      seams::product::reduction() != seams::product::Reduction::Left && orderL <= 8;
  // Same split as runPass1: degree 3 stays on the one-bond kernel.
  const bool batchSolid = (orderL == 4 || orderL == 6) && !pairwise;
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold && omp_get_max_threads() > 1;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = 0; i < g.nop; i++) {
    if (which != nullptr &&
        (static_cast<int>(which->size()) <= i || (*which)[static_cast<std::size_t>(i)] == 0)) {
      continue;
    }
    if (pairwise) {
      qlmOneAtomDrPairwise(i, orderL, g.dr.data(), g.offsets.data(), qlm.data());
    } else if (batchSolid) {
      qlmOneAtomBatch(i, orderL, g.dr.data(), g.offsets.data(), qlm.data());
    } else {
      seams::steinhardt::qlmOneAtomDr(i, orderL, g.dr.data(), g.offsets.data(),
                                      g.cols.data(), qlm.data());
    }
  }
}

void runPass2Mask(const NeighbourCSR &g, int orderL, const std::vector<double> &qlm,
                  std::vector<double> &ql, std::vector<double> &qlBar,
                  const std::vector<char> *which) {
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold && omp_get_max_threads() > 1;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = 0; i < g.nop; i++) {
    if (which != nullptr &&
        (static_cast<int>(which->size()) <= i || (*which)[static_cast<std::size_t>(i)] == 0)) {
      continue;
    }
    seams::steinhardt::qlOneAtom(i, orderL, qlm.data(), g.offsets.data(), g.cols.data(),
                                 ql.data(), qlBar.data());
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
    // l<=12 so nComp<=25; stack, not a heap vector per atom
    double barRe[25] = {};
    double barIm[25] = {};
    for (int m = 0; m < nComp; m++) {
      barRe[m] = qlm[static_cast<size_t>(2 * (row + m))];
      barIm[m] = qlm[static_cast<size_t>(2 * (row + m) + 1)];
    }
    int nContrib = 1;
    for (int p = j0; p < j1; p++) {
      const int jatom = g.cols[static_cast<size_t>(p)];
      const int jRow = jatom * nComp;
      for (int m = 0; m < nComp; m++) {
        barRe[m] += qlm[static_cast<size_t>(2 * (jRow + m))];
        barIm[m] += qlm[static_cast<size_t>(2 * (jRow + m) + 1)];
      }
      nContrib++;
    }
    const double inv = 1.0 / static_cast<double>(nContrib);
    double sumBar = 0.0;
    for (int m = 0; m < nComp; m++) {
      const double re = barRe[m] * inv;
      const double im = barIm[m] * inv;
      sumBar += re * re + im * im;
    }
    qlBar[static_cast<size_t>(i)] = std::sqrt(prefactor * sumBar);
  }
}

void runPass2(const NeighbourCSR &g, int orderL, int begin, int end,
              const std::vector<double> &qlm, std::vector<double> &ql,
              std::vector<double> &qlBar) {
  const bool stock =
      seams::product::mode() == seams::product::Mode::Fp64 &&
      seams::product::reduction() == seams::product::Reduction::Left;
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold;
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = begin; i < end; i++) {
    if (stock) {
      seams::steinhardt::qlOneAtom(i, orderL, qlm.data(), g.offsets.data(),
                                   g.cols.data(), ql.data(), qlBar.data());
      continue;
    }
    const int nComp = 2 * orderL + 1;
    constexpr double pi = 3.14159265358979323846;
    const double prefactor = 4.0 * pi / (2.0 * orderL + 1.0);
    const int row = i * nComp;
    const auto norm = seams::product::realDot(qlm.data() + 2 * row, qlm.data() + 2 * row,
                                              2 * nComp, seams::product::mode());
    ql[static_cast<size_t>(i)] = std::sqrt(prefactor * norm.value);
    const int j0 = g.offsets[static_cast<size_t>(i)];
    const int j1 = g.offsets[static_cast<size_t>(i) + 1];
    if (j0 == j1) {
      qlBar[static_cast<size_t>(i)] = ql[static_cast<size_t>(i)];
      continue;
    }
    const int nContrib = 1 + (j1 - j0);
    std::vector<double> comp(static_cast<size_t>(nContrib));
    std::vector<double> bar(static_cast<size_t>(2 * nComp));
    for (int m = 0; m < 2 * nComp; m++) {
      comp[0] = qlm[static_cast<size_t>(2 * row + m)];
      for (int p = j0; p < j1; p++) {
        const int jatom = g.cols[static_cast<size_t>(p)];
        comp[static_cast<size_t>(1 + (p - j0))] =
            qlm[static_cast<size_t>(2 * (jatom * nComp) + m)];
      }
      bar[static_cast<size_t>(m)] =
          seams::product::reduce(comp.data(), nContrib, seams::product::reduction()) /
          static_cast<double>(nContrib);
    }
    const auto barNorm = seams::product::realDot(bar.data(), bar.data(), 2 * nComp,
                                                 seams::product::mode());
    qlBar[static_cast<size_t>(i)] = std::sqrt(prefactor * barNorm.value);
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

bool atomSplit = true;
#endif

} // namespace

namespace chill {

void setSteinhardtAtomSplit([[maybe_unused]] bool split) {
#ifdef SEAMS_HAS_MPI
  atomSplit = split;
#endif
}

SteinhardtQl steinhardtQl(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                          const std::vector<std::vector<int>> &nList, int orderL) {
  // OpenMP only. adoptThreads would pause the team and start a Rayon
  // pool this harmonic never uses.
  seams::prepareCompute();
  [[maybe_unused]] const seams::TeamWidth team(yCloud.nop);
  SteinhardtQl result;
  result.ql.assign(static_cast<std::size_t>(yCloud.nop), 0.0);
  result.qlBar.assign(static_cast<std::size_t>(yCloud.nop), 0.0);

  if (orderL != 3 && orderL != 4 && orderL != 6 && orderL != 8 &&
      orderL != 12) {
    return result;
  }
  if (yCloud.nop <= 0) {
    return result;
  }

  const NeighbourCSR graph = flatten(yCloud, nList);
  const int nComp = 2 * orderL + 1;
  std::vector<double> qlm(static_cast<size_t>(graph.nop) * static_cast<std::size_t>(nComp) * 2, 0.0);

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
  // A frame-split driver turns the split off: its ranks hold different
  // frames, and under FUNNELED its frame workers may not call the
  // communicator
  int initialized = 0;
  if (atomSplit) {
    MPI_Initialized(&initialized);
  }
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

namespace {

void localQlOnly(const NeighbourCSR &g, int orderL, const std::vector<double> &qlm,
                 std::vector<double> &ql) {
  const int nComp = 2 * orderL + 1;
  constexpr double pi = 3.14159265358979323846;
  const double prefactor = 4.0 * pi / (2.0 * orderL + 1.0);
#ifdef SEAMS_HAS_OPENMP
  const bool useThreads = g.nop >= kParallelThreshold && !omp_in_parallel();
#pragma omp parallel for schedule(static) if (useThreads)
#endif
  for (int i = 0; i < g.nop; ++i) {
    const int row = i * nComp;
    double sum = 0.0;
    for (int m = 0; m < nComp; ++m) {
      const double re = qlm[static_cast<std::size_t>(2 * (row + m))];
      const double im = qlm[static_cast<std::size_t>(2 * (row + m) + 1)];
      sum += re * re + im * im;
    }
    ql[static_cast<std::size_t>(i)] = std::sqrt(prefactor * sum);
  }
}

bool rowLists(const int *nom, int k, int row, int target) {
  const int *base = nom + static_cast<std::size_t>(row) * static_cast<std::size_t>(k);
  for (int t = 0; t < k; ++t) {
    if (base[t] == target) {
      return true;
    }
  }
  return false;
}

} // namespace

// Cutoff shell, solid harmonic, buffers kept across calls. Ortho boxes
// only. The mutual four-nearest nomination is a different graph and a
// heavier search; on ice the cutoff shell is those four neighbours.
SteinhardtQl cutoffSolid(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int orderL, int typeI, bool average, bool invariants,
    std::vector<double> *q6out = nullptr, bool rebuild = true) {
  const int n = yCloud.nop;
  SteinhardtQl result;
  result.ql.assign(static_cast<std::size_t>(n), 0.0);
  result.qlBar.assign(static_cast<std::size_t>(n), 0.0);
  const int nComp = 2 * orderL + 1;
  const int width = 2 * nComp;
  struct Scratch {
    std::vector<int> offsets, cols;
    std::vector<double> xyz, qlm;
    std::vector<double> box, boxLow;
    int n = -1;
    int orderL = -1;
    int geomN = -1;
    int geomType = 0;
    double geomCut = 0;
  };
  static Scratch scratch;
  // The fast local q4/q6 kernel can skip the row. Every other degree
  // writes q_lm, including a local q3, so the row exists before that
  // loop. A wider buffer left by a higher degree is reused; growing is
  // the only allocation.
  if (scratch.xyz.size() != static_cast<std::size_t>(n) * 3) {
    scratch.xyz.resize(static_cast<std::size_t>(n) * 3);
  }
  const std::size_t qlmNeed =
      static_cast<std::size_t>(n) * static_cast<std::size_t>(width);
  if (scratch.qlm.size() < qlmNeed) {
    scratch.qlm.assign(qlmNeed, 0.0);
  }
  scratch.n = n;
  scratch.orderL = orderL;
  // A PDB with no CRYST1 line has no box. Indexing it is a fault, and a
  // cutoff on a missing cell is an empty q_l.
  if (yCloud.box.size() < 3 || !(yCloud.box[0] > 0.0 && yCloud.box[1] > 0.0 &&
                                 yCloud.box[2] > 0.0)) {
    if (q6out != nullptr) {
      q6out->clear();
    }
    return {};
  }
  const double ox = yCloud.boxLow.size() > 0 ? yCloud.boxLow[0] : 0.0;
  const double oy = yCloud.boxLow.size() > 1 ? yCloud.boxLow[1] : 0.0;
  const double oz = yCloud.boxLow.size() > 2 ? yCloud.boxLow[2] : 0.0;
  const double Lx = yCloud.box[0];
  const double Ly = yCloud.box[1];
  const double Lz = yCloud.box[2];
  const double hx = 0.5 * Lx;
  const double hy = 0.5 * Ly;
  const double hz = 0.5 * Lz;
  const bool tilt = yCloud.box.size() >= 6;
  const gen::FracBox frac = gen::makeFracBox(yCloud);
  bool uniform = true;
  for (int i = 0; i < n; i++) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type != typeI) {
      uniform = false;
      break;
    }
  }
  std::vector<int> subset;
  const std::vector<int> *which = nullptr;
  if (!uniform) {
    subset.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
      if (yCloud.pts[static_cast<std::size_t>(i)].type == typeI) {
        subset.push_back(i);
      }
    }
    which = &subset;
  }
#ifdef SEAMS_HAS_OPENMP
  const int saved = std::max(1, omp_get_max_threads());
  omp_set_num_threads(seams::teamFor(n));
#endif
  const bool timing = std::getenv("SEAMS_QL_TIME") != nullptr;
  const auto t0 = std::chrono::steady_clock::now();
  const bool sameGeom = !rebuild && scratch.geomN == n && scratch.geomType == typeI &&
                        scratch.geomCut == cutoff && scratch.box == yCloud.box &&
                        scratch.boxLow == yCloud.boxLow &&
                        scratch.offsets.size() == static_cast<std::size_t>(n) + 1;
  bool have = true;
  if (!sameGeom) {
    have = nneigh::orthoFullCsr(yCloud, which, cutoff, scratch.offsets, scratch.cols);
    scratch.geomN = -1;
  }
  const auto tCsr = std::chrono::steady_clock::now();
#ifdef SEAMS_HAS_OPENMP
  const bool parallel = seams::teamFor(n) > 1;
#endif
  if (!sameGeom) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (int i = 0; i < n; i++) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(i)];
    if (tilt) {
      // Reduced coordinates. The bond loop then wraps those lengths
      // instead of repeating the fractional conversion on every neighbour.
      double sz = (p.z - frac.oz) / frac.lz;
      sz -= std::floor(sz);
      if (!(sz < 1.0)) sz = 0.0;
      double sy = (p.y - frac.oy - frac.yz * sz) / frac.ly;
      sy -= std::floor(sy);
      if (!(sy < 1.0)) sy = 0.0;
      double sx = (p.x - frac.ox - frac.xy * sy - frac.xz * sz) / frac.lx;
      sx -= std::floor(sx);
      if (!(sx < 1.0)) sx = 0.0;
      scratch.xyz[static_cast<std::size_t>(3 * i)] = sx * frac.lx;
      scratch.xyz[static_cast<std::size_t>(3 * i) + 1] = sy * frac.ly;
      scratch.xyz[static_cast<std::size_t>(3 * i) + 2] = sz * frac.lz;
    } else {
      auto wrap = [](double v, double origin, double L) {
        double t = (v - origin) / L;
        t -= std::floor(t);
        if (!(t < 1.0)) {
          t = 0.0;
        }
        return t * L;
      };
      scratch.xyz[static_cast<std::size_t>(3 * i)] = wrap(p.x, ox, Lx);
      scratch.xyz[static_cast<std::size_t>(3 * i) + 1] = wrap(p.y, oy, Ly);
      scratch.xyz[static_cast<std::size_t>(3 * i) + 2] = wrap(p.z, oz, Lz);
    }
  }
    if (have) {
      scratch.geomN = n;
      scratch.geomType = typeI;
      scratch.geomCut = cutoff;
      scratch.box = yCloud.box;
      scratch.boxLow = yCloud.boxLow;
    }
  }
  if (!have) {
#ifdef SEAMS_HAS_OPENMP
    omp_set_num_threads(saved);
#endif
    if (q6out != nullptr) {
      q6out->clear();
    }
    return {};
  }
  constexpr double pi = 3.14159265358979323846;
  const double prefactor = 4.0 * pi / (2.0 * orderL + 1.0);
  const double axy = tilt ? frac.xy / frac.ly : 0.0;
  const double axz = tilt ? frac.xz / frac.lz : 0.0;
  const double ayz = tilt ? frac.yz / frac.lz : 0.0;
  // Local q4 and q6. Eight atoms share one solid-harmonic sum. The
  // neighbour shell is whatever the cutoff list stored, not a fixed
  // four, and the q_lm row is not written.
  bool wroteLocal = false;
  if (!average && !invariants && uniform && (orderL == 4 || orderL == 6)) {
    wroteLocal = true;
    // The eight-wide sum declines a tile the CPU cannot take. The scalar
    // fallback writes q_lm, so that row exists even when the average does not.
    if (scratch.qlm.size() !=
        static_cast<std::size_t>(n) * static_cast<std::size_t>(width)) {
      scratch.qlm.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(width), 0.0);
    }
    const int nTiles = (n + 7) / 8;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
    for (int tile = 0; tile < nTiles; ++tile) {
      const int begin = tile * 8;
      const int alive = std::min(8, n - begin);
      int nUsed[8] = {};
      int upper = 0;
      for (int a = 0; a < alive; ++a) {
        const int i = begin + a;
        const int span = scratch.offsets[static_cast<std::size_t>(i) + 1] -
                         scratch.offsets[static_cast<std::size_t>(i)];
        if (span > upper) {
          upper = span;
        }
      }
      thread_local std::vector<double> dxs, dys, dzs;
      if (static_cast<int>(dxs.size()) < upper * 8) {
        dxs.assign(static_cast<std::size_t>(upper * 8), 0.0);
        dys.assign(static_cast<std::size_t>(upper * 8), 0.0);
        dzs.assign(static_cast<std::size_t>(upper * 8), 0.0);
      }
      int maxB = 0;
      for (int a = 0; a < alive; ++a) {
        const int i = begin + a;
        const int left = scratch.offsets[static_cast<std::size_t>(i)];
        const int right = scratch.offsets[static_cast<std::size_t>(i) + 1];
        const double xi = scratch.xyz[static_cast<std::size_t>(3 * i)];
        const double yi = scratch.xyz[static_cast<std::size_t>(3 * i) + 1];
        const double zi = scratch.xyz[static_cast<std::size_t>(3 * i) + 2];
        int filled = 0;
        for (int p = left; p < right; ++p) {
          const int j = scratch.cols[static_cast<std::size_t>(p)];
          double dx, dy, dz;
          if (tilt) {
            double dX = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
            double dY = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
            double dZ = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
            dX -= frac.lx * std::round(dX / frac.lx);
            dY -= frac.ly * std::round(dY / frac.ly);
            dZ -= frac.lz * std::round(dZ / frac.lz);
            dx = dX + axy * dY + axz * dZ;
            dy = dY + ayz * dZ;
            dz = dZ;
          } else {
            dx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
            dy = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
            dz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
            if (dx > hx) dx -= Lx;
            else if (dx < -hx) dx += Lx;
            if (dy > hy) dy -= Ly;
            else if (dy < -hy) dy += Ly;
            if (dz > hz) dz -= Lz;
            else if (dz < -hz) dz += Lz;
          }
          if (!(dx * dx + dy * dy + dz * dz > 0.0)) {
            continue;
          }
          dxs[static_cast<std::size_t>(filled * 8 + a)] = dx;
          dys[static_cast<std::size_t>(filled * 8 + a)] = dy;
          dzs[static_cast<std::size_t>(filled * 8 + a)] = dz;
          ++filled;
        }
        nUsed[a] = filled;
        if (filled > maxB) {
          maxB = filled;
        }
      }
      for (int a = 0; a < alive; ++a) {
        for (int bond = nUsed[a]; bond < maxB; ++bond) {
          dxs[static_cast<std::size_t>(bond * 8 + a)] = 0.0;
          dys[static_cast<std::size_t>(bond * 8 + a)] = 0.0;
          dzs[static_cast<std::size_t>(bond * 8 + a)] = 0.0;
        }
      }
      for (int a = alive; a < 8; ++a) {
        nUsed[a] = 0;
      }
      double tileQl[8] = {};
      double tileQ6[8] = {};
      // q4 and q6 share this shell. A second call would rebuild the cutoff.
      const bool also6 = q6out != nullptr && orderL == 4;
      const double pre6 = 4.0 * pi / 13.0;
      if (maxB > 0) {
        seams::ylm::ylmLocalQl(orderL, dxs.data(), dys.data(), dzs.data(), maxB, nUsed,
                               prefactor, tileQl);
        if (also6) {
          seams::ylm::ylmLocalQl(6, dxs.data(), dys.data(), dzs.data(), maxB, nUsed, pre6,
                                 tileQ6);
        }
      }
      for (int a = 0; a < alive; ++a) {
        result.ql[static_cast<std::size_t>(begin + a)] = tileQl[a];
        if (also6) {
          (*q6out)[static_cast<std::size_t>(begin + a)] = tileQ6[a];
        }
      }
    }
  }
  // Pairing two atoms into one AVX-512 reduction was slower than the
  // four-bond kernel. The horizontal sum cost more than the idle lanes.
  bool quad = false && uniform && (orderL == 4 || orderL == 6);
  for (int i = 0; quad && i < n; ++i) {
    if (scratch.offsets[static_cast<std::size_t>(i) + 1] -
            scratch.offsets[static_cast<std::size_t>(i)] !=
        4) {
      quad = false;
    }
  }
#if defined(SEAMS_YLM_AVX512)
  const bool simd8 = quad && __builtin_cpu_supports("avx512f") &&
                     __builtin_cpu_supports("avx512dq");
#else
  const bool simd8 = false;
  (void)quad;
#endif
  if (!wroteLocal && simd8) {
    auto oneAtom = [&](int atom) {
      const int a = scratch.offsets[static_cast<std::size_t>(atom)];
      const double xi = scratch.xyz[static_cast<std::size_t>(3 * atom)];
      const double yi = scratch.xyz[static_cast<std::size_t>(3 * atom) + 1];
      const double zi = scratch.xyz[static_cast<std::size_t>(3 * atom) + 2];
      double xs[4], ys[4], zs[4], re[4 * 13], im[4 * 13];
      int nUsed = 0;
      for (int p = 0; p < 4; ++p) {
        const int j = scratch.cols[static_cast<std::size_t>(a + p)];
        double dx, dy, dz;
        if (tilt) {
          double dX = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
          double dY = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
          double dZ = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
          dX -= frac.lx * std::round(dX / frac.lx);
          dY -= frac.ly * std::round(dY / frac.ly);
          dZ -= frac.lz * std::round(dZ / frac.lz);
          dx = dX + axy * dY + axz * dZ;
          dy = dY + ayz * dZ;
          dz = dZ;
        } else {
          dx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
          dy = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
          dz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
          if (dx > hx) dx -= Lx;
          else if (dx < -hx) dx += Lx;
          if (dy > hy) dy -= Ly;
          else if (dy < -hy) dy += Ly;
          if (dz > hz) dz -= Lz;
          else if (dz < -hz) dz += Lz;
        }
        if (!(dx * dx + dy * dy + dz * dz > 0.0)) {
          continue;
        }
        xs[nUsed] = dx;
        ys[nUsed] = dy;
        zs[nUsed] = dz;
        ++nUsed;
      }
      double *dst = scratch.qlm.data() +
                    static_cast<std::size_t>(atom) * static_cast<std::size_t>(width);
      std::fill_n(dst, width, 0.0);
      if (nUsed == 0) {
        return;
      }
      seams::ylm::ylmBatch(orderL, xs, ys, zs, nUsed, re, im);
      for (int bond = 0; bond < nUsed; ++bond) {
        for (int m = 0; m < nComp; ++m) {
          dst[2 * m] += re[bond * nComp + m];
          dst[2 * m + 1] += im[bond * nComp + m];
        }
      }
      const double inv = 1.0 / static_cast<double>(nUsed);
      double sum = 0.0;
      for (int m = 0; m < nComp; ++m) {
        dst[2 * m] *= inv;
        dst[2 * m + 1] *= inv;
        sum += dst[2 * m] * dst[2 * m] + dst[2 * m + 1] * dst[2 * m + 1];
      }
      result.ql[static_cast<std::size_t>(atom)] = std::sqrt(prefactor * sum);
    };
    auto pairRange = [&](int begin, int end) {
      for (int i = begin; i + 1 < end; i += 2) {
        double xs[8], ys[8], zs[8];
        int slot = 0;
        for (int atom = i; atom < i + 2; ++atom) {
          const int a = scratch.offsets[static_cast<std::size_t>(atom)];
          const double xi = scratch.xyz[static_cast<std::size_t>(3 * atom)];
          const double yi = scratch.xyz[static_cast<std::size_t>(3 * atom) + 1];
          const double zi = scratch.xyz[static_cast<std::size_t>(3 * atom) + 2];
          for (int p = 0; p < 4; ++p) {
            const int j = scratch.cols[static_cast<std::size_t>(a + p)];
            double dx, dy, dz;
            if (tilt) {
              double dX = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
              double dY = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
              double dZ = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
              dX -= frac.lx * std::round(dX / frac.lx);
              dY -= frac.ly * std::round(dY / frac.ly);
              dZ -= frac.lz * std::round(dZ / frac.lz);
              dx = dX + axy * dY + axz * dZ;
              dy = dY + ayz * dZ;
              dz = dZ;
            } else {
              dx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
              dy = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
              dz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
              if (dx > hx) dx -= Lx;
              else if (dx < -hx) dx += Lx;
              if (dy > hy) dy -= Ly;
              else if (dy < -hy) dy += Ly;
              if (dz > hz) dz -= Lz;
              else if (dz < -hz) dz += Lz;
            }
            xs[slot] = dx;
            ys[slot] = dy;
            zs[slot] = dz;
            ++slot;
          }
        }
        double accA[26] = {};
        double accB[26] = {};
        seams::ylm::ylmAdd8(orderL, xs, ys, zs, accA, accB);
        auto commit = [&](int atom, const double *acc) {
          double *dst = scratch.qlm.data() +
                        static_cast<std::size_t>(atom) * static_cast<std::size_t>(width);
          constexpr double inv = 0.25;
          double sum = 0.0;
          for (int m = 0; m < nComp; ++m) {
            dst[2 * m] = acc[2 * m] * inv;
            dst[2 * m + 1] = acc[2 * m + 1] * inv;
            sum += dst[2 * m] * dst[2 * m] + dst[2 * m + 1] * dst[2 * m + 1];
          }
          result.ql[static_cast<std::size_t>(atom)] = std::sqrt(prefactor * sum);
        };
        commit(i, accA);
        commit(i + 1, accB);
      }
    };
#ifdef SEAMS_HAS_OPENMP
    if (parallel) {
#pragma omp parallel
      {
        const int nt = std::max(1, omp_get_num_threads());
        const int tid = omp_get_thread_num();
        const int begin = static_cast<int>((static_cast<long long>(n) * tid) / nt);
        const int end = static_cast<int>((static_cast<long long>(n) * (tid + 1)) / nt);
        const int endEven = begin + ((end - begin) & ~1);
        pairRange(begin, endEven);
        if (endEven < end) {
          oneAtom(endEven);
        }
      }
    } else
#endif
    {
      const int endEven = n & ~1;
      pairRange(0, endEven);
      if (endEven < n) {
        oneAtom(endEven);
      }
    }
  }
  if (!wroteLocal && !simd8) {
  const int nTiles = (n + 7) / 8;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
  for (int tile = 0; tile < nTiles; ++tile) {
    const int begin = tile * 8;
    const int end = std::min(begin + 8, n);
    // Four bonds on eight atoms: the harmonics stay in SoA and one
    // blocked transpose writes each atom row. A per-bond lane
    // transpose was the cost this replaces.
    bool rows = end - begin == 8 && (orderL == 4 || orderL == 6);
#if !defined(SEAMS_YLM_AVX512)
    rows = false;
#endif
    alignas(64) double dxs[32];
    alignas(64) double dys[32];
    alignas(64) double dzs[32];
    for (int a = 0; rows && a < 8; ++a) {
      const int i = begin + a;
      if (!uniform && yCloud.pts[static_cast<std::size_t>(i)].type != typeI) {
        rows = false;
        break;
      }
      const int left = scratch.offsets[static_cast<std::size_t>(i)];
      const int right = scratch.offsets[static_cast<std::size_t>(i) + 1];
      if (right - left != 4) {
        rows = false;
        break;
      }
      const double xi = scratch.xyz[static_cast<std::size_t>(3 * i)];
      const double yi = scratch.xyz[static_cast<std::size_t>(3 * i) + 1];
      const double zi = scratch.xyz[static_cast<std::size_t>(3 * i) + 2];
      for (int p = 0; p < 4; ++p) {
        const int j = scratch.cols[static_cast<std::size_t>(left + p)];
        double dx, dy, dz;
        if (tilt) {
          double dX = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
          double dY = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
          double dZ = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
          dX -= frac.lx * std::round(dX / frac.lx);
          dY -= frac.ly * std::round(dY / frac.ly);
          dZ -= frac.lz * std::round(dZ / frac.lz);
          dx = dX + axy * dY + axz * dZ;
          dy = dY + ayz * dZ;
          dz = dZ;
        } else {
          dx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
          dy = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
          dz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
          if (dx > hx) dx -= Lx;
          else if (dx < -hx) dx += Lx;
          if (dy > hy) dy -= Ly;
          else if (dy < -hy) dy += Ly;
          if (dz > hz) dz -= Lz;
          else if (dz < -hz) dz += Lz;
        }
        if (!(dx * dx + dy * dy + dz * dz > 0.0)) {
          rows = false;
          break;
        }
        dxs[p * 8 + a] = dx;
        dys[p * 8 + a] = dy;
        dzs[p * 8 + a] = dz;
      }
    }
#if defined(SEAMS_YLM_AVX512)
    if (rows && !__builtin_cpu_supports("avx512f")) {
      rows = false;
    }
    if (rows) {
      double *row = scratch.qlm.data() +
                    static_cast<std::size_t>(begin) * static_cast<std::size_t>(width);
      rows = seams::ylm::ylmSumRows4(orderL, dxs, dys, dzs, row, width);
      if (rows) {
        constexpr double inv = 0.25;
        for (int a = 0; a < 8; ++a) {
          double *dst = row + static_cast<std::size_t>(a) * static_cast<std::size_t>(width);
          double sum = 0.0;
          for (int m = 0; m < nComp; ++m) {
            dst[2 * m] *= inv;
            dst[2 * m + 1] *= inv;
            sum += dst[2 * m] * dst[2 * m] + dst[2 * m + 1] * dst[2 * m + 1];
          }
          result.ql[static_cast<std::size_t>(begin + a)] = std::sqrt(prefactor * sum);
        }
      }
    }
#endif
    if (rows) {
      continue;
    }
    for (int i = begin; i < end; ++i) {
    double *dst = scratch.qlm.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(width);
    std::fill_n(dst, width, 0.0);
    const bool live = uniform || yCloud.pts[static_cast<std::size_t>(i)].type == typeI;
    if (!live) {
      continue;
    }
    const int a = scratch.offsets[static_cast<std::size_t>(i)];
    const int b = scratch.offsets[static_cast<std::size_t>(i) + 1];
    const double xi = scratch.xyz[static_cast<std::size_t>(3 * i)];
    const double yi = scratch.xyz[static_cast<std::size_t>(3 * i) + 1];
    const double zi = scratch.xyz[static_cast<std::size_t>(3 * i) + 2];
    double xs[8], ys[8], zs[8], re[8 * 13], im[8 * 13];
    int filled = 0;
    int nUsed = 0;
    const auto flush = [&]() {
      if (filled == 0) {
        return;
      }
      seams::ylm::ylmBatch(orderL, xs, ys, zs, filled, re, im);
      for (int bond = 0; bond < filled; bond++) {
        const double *yr = re + static_cast<std::size_t>(bond) * static_cast<std::size_t>(nComp);
        const double *imag = im + static_cast<std::size_t>(bond) * static_cast<std::size_t>(nComp);
        for (int m = 0; m < nComp; m++) {
          dst[2 * m] += yr[m];
          dst[2 * m + 1] += imag[m];
        }
      }
      filled = 0;
    };
    for (int p = a; p < b; p++) {
      const int j = scratch.cols[static_cast<std::size_t>(p)];
      double dx, dy, dz;
      if (tilt) {
        double dX = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
        double dY = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
        double dZ = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
        dX -= frac.lx * std::round(dX / frac.lx);
        dY -= frac.ly * std::round(dY / frac.ly);
        dZ -= frac.lz * std::round(dZ / frac.lz);
        dx = dX + axy * dY + axz * dZ;
        dy = dY + ayz * dZ;
        dz = dZ;
      } else {
        dx = xi - scratch.xyz[static_cast<std::size_t>(3 * j)];
        dy = yi - scratch.xyz[static_cast<std::size_t>(3 * j) + 1];
        dz = zi - scratch.xyz[static_cast<std::size_t>(3 * j) + 2];
        if (dx > hx) dx -= Lx;
        else if (dx < -hx) dx += Lx;
        if (dy > hy) dy -= Ly;
        else if (dy < -hy) dy += Ly;
        if (dz > hz) dz -= Lz;
        else if (dz < -hz) dz += Lz;
      }
      if (!(dx * dx + dy * dy + dz * dz > 0.0)) {
        continue;
      }
      xs[filled] = dx;
      ys[filled] = dy;
      zs[filled] = dz;
      ++nUsed;
      if (++filled == 8) {
        flush();
      }
    }
    flush();
    if (nUsed == 0) {
      continue;
    }
    const double inv = 1.0 / static_cast<double>(nUsed);
    double sum = 0.0;
    for (int m = 0; m < nComp; m++) {
      dst[2 * m] *= inv;
      dst[2 * m + 1] *= inv;
      sum += dst[2 * m] * dst[2 * m] + dst[2 * m + 1] * dst[2 * m + 1];
    }
    result.ql[static_cast<std::size_t>(i)] = std::sqrt(prefactor * sum);
    }
  }
  }
  if (invariants && (orderL == 4 || orderL == 6)) {
    result.wl.assign(static_cast<std::size_t>(n), 0.0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
    for (int i = 0; i < n; i++) {
      const double *row = scratch.qlm.data() +
                          static_cast<std::size_t>(i) * static_cast<std::size_t>(width);
      result.wl[static_cast<std::size_t>(i)] = seams::wigner::hatW(orderL, row);
    }
  }
  const auto tYlm = std::chrono::steady_clock::now();
  if (average) {
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (parallel)
#endif
    for (int i = 0; i < n; i++) {
      const int a = scratch.offsets[static_cast<std::size_t>(i)];
      const int b = scratch.offsets[static_cast<std::size_t>(i) + 1];
      const double *self = scratch.qlm.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(width);
      double acc[26];
      for (int c = 0; c < width; ++c) acc[c] = self[c];
      for (int p = a; p < b; p++) {
        const int j = scratch.cols[static_cast<std::size_t>(p)];
        const double *other = scratch.qlm.data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(width);
        for (int c = 0; c < width; ++c) acc[c] += other[c];
      }
      const double inv = 1.0 / static_cast<double>(b - a + 1);
      double sum = 0.0;
      for (int c = 0; c < width; ++c) {
        acc[c] *= inv;
        sum += acc[c] * acc[c];
      }
      result.qlBar[static_cast<std::size_t>(i)] = std::sqrt(prefactor * sum);
    }
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
  if (q6out != nullptr && !wroteLocal) {
    q6out->clear();
  }
  if (timing) {
    const auto tEnd = std::chrono::steady_clock::now();
    const auto ms = [](auto a, auto b) {
      return std::chrono::duration<double, std::milli>(b - a).count();
    };
    std::cerr << "ql-csr l=" << orderL << " n=" << n << " L=" << Lx
              << " pairs=" << scratch.cols.size() << " team=" << seams::teamFor(n)
              << " csr_ms=" << ms(t0, tCsr) << " ylm_ms=" << ms(tCsr, tYlm)
              << " bar_ms=" << ms(tYlm, tEnd) << "\n";
  }
  return result;
}

bool localCutoffQ4Q6(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI, std::vector<double> &q4, std::vector<double> &q6) {
  q6.assign(static_cast<std::size_t>(std::max(yCloud.nop, 0)), 0.0);
  SteinhardtQl got = cutoffSolid(yCloud, cutoff, 4, typeI, false, false, &q6);
  if (got.ql.empty() || q6.size() != got.ql.size()) {
    return false;
  }
  q4 = std::move(got.ql);
  return true;
}

bool steinhardtQ4Q6Cutoff(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int typeI, bool average, SteinhardtQl &q4, SteinhardtQl &q6) {
  SteinhardtQl a4 = cutoffSolid(yCloud, cutoff, 4, typeI, average, false, nullptr, true);
  if (a4.ql.empty()) {
    return false;
  }
  SteinhardtQl a6 = cutoffSolid(yCloud, cutoff, 6, typeI, average, false, nullptr, false);
  if (a6.ql.empty()) {
    return false;
  }
  q4 = std::move(a4);
  q6 = std::move(a6);
  return true;
}

SteinhardtQl steinhardtQlCutoff(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double cutoff,
    int orderL, int typeI, bool average, bool invariants) {
  if ((orderL == 3 || orderL == 4 || orderL == 6) &&
      (yCloud.box.size() == 3 || yCloud.box.size() >= 6) && yCloud.nop > 0 &&
      cutoff > 0.0) {
    SteinhardtQl fast = cutoffSolid(yCloud, cutoff, orderL, typeI, average, invariants);
    if (!fast.ql.empty()) {
      return fast;
    }
  }
  SteinhardtQl result;
  result.ql.assign(static_cast<std::size_t>(std::max(yCloud.nop, 0)), 0.0);
  result.qlBar.assign(result.ql.size(), 0.0);
  if (orderL != 3 && orderL != 4 && orderL != 6 && orderL != 8 && orderL != 12) {
    return result;
  }
  if (yCloud.nop <= 0 || !(cutoff > 0.0)) {
    return result;
  }

  std::vector<int> subset;
  subset.reserve(static_cast<std::size_t>(yCloud.nop));
  for (int i = 0; i < yCloud.nop; ++i) {
    if (yCloud.pts[static_cast<std::size_t>(i)].type == typeI) {
      subset.push_back(i);
    }
  }
  std::vector<int> offsets;
  std::vector<int> cols;
  const bool have = nneigh::orthoFullCsr(yCloud, &subset, cutoff, offsets, cols);
  if (!have) {
    const auto rows = nneigh::neighListO(cutoff, yCloud, typeI);
    result = steinhardtQl(yCloud, rows, orderL);
    if (!average) {
      std::fill(result.qlBar.begin(), result.qlBar.end(), 0.0);
    }
    return result;
  }

  NeighbourCSR g;
  g.nop = yCloud.nop;
  g.offsets = std::move(offsets);
  g.cols = std::move(cols);
  g.dr.assign(g.cols.size() * 3, 0.0);
  const gen::FracBox frame = gen::makeFracBox(yCloud);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (g.nop >= kParallelThreshold && !omp_in_parallel())
#endif
  for (int i = 0; i < g.nop; ++i) {
    const int a = g.offsets[static_cast<std::size_t>(i)];
    const int b = g.offsets[static_cast<std::size_t>(i) + 1];
    for (int k = a; k < b; ++k) {
      const int j = g.cols[static_cast<std::size_t>(k)];
      const auto &pj = yCloud.pts[static_cast<std::size_t>(j)];
      const auto d = gen::relDistFromPoint(frame, yCloud, i, pj.x, pj.y, pj.z);
      g.dr[static_cast<std::size_t>(3 * k)] = d[0];
      g.dr[static_cast<std::size_t>(3 * k + 1)] = d[1];
      g.dr[static_cast<std::size_t>(3 * k + 2)] = d[2];
    }
  }

  const int nComp = 2 * orderL + 1;
  std::vector<double> qlm(static_cast<std::size_t>(g.nop) * static_cast<std::size_t>(nComp) * 2,
                          0.0);
#ifdef SEAMS_HAS_OFFLOAD
  if (average && orderL != 12 && wantOffload()) {
    runOffload(g, orderL, qlm, result.ql, result.qlBar);
    return result;
  }
#endif
  runPass1(g, orderL, 0, g.nop, qlm);
  if (!average) {
    localQlOnly(g, orderL, qlm, result.ql);
    return result;
  }
  if (orderL == 12) {
    runPass2Host(g, orderL, 0, g.nop, qlm, result.ql, result.qlBar);
  } else {
    runPass2(g, orderL, 0, g.nop, qlm, result.ql, result.qlBar);
  }
  return result;
}

SteinhardtQl steinhardtQlMutual(
    const molSys::PointCloud<molSys::Point<double>, double> &yCloud, int k,
    double candidateCutoff, int typeI, int orderL) {
  const int n = yCloud.nop;
  if (n <= 0 || k <= 0 || k > 32 ||
      (orderL != 3 && orderL != 4 && orderL != 6)) {
    const auto rows =
        nneigh::kNearestNeighbourList(yCloud, std::max(k, 0), candidateCutoff, typeI, true);
    return steinhardtQl(yCloud, rows, orderL);
  }

  const bool timing = std::getenv("SEAMS_QL_TIME") != nullptr;
  const auto t0 = std::chrono::steady_clock::now();
  // Repeated frames of one size reuse the buffers. The benchmark calls
  // this once per sample; allocating the q_lm row every time was visible.
  struct Scratch {
    std::vector<double> xyz;
    std::vector<int> nom;
    std::vector<int> nkeep;
    std::vector<int> neigh;
    std::vector<double> qlm;
    int n = -1;
    int k = -1;
    int orderL = -1;
  };
  static Scratch scratch;
  if (scratch.n != n || scratch.k != k || scratch.orderL != orderL) {
    scratch.nkeep.assign(static_cast<std::size_t>(n), 0);
    scratch.neigh.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(k), -1);
    scratch.qlm.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(2 * (2 * orderL + 1)),
                       0.0);
    scratch.n = n;
    scratch.k = k;
    scratch.orderL = orderL;
  }
  nneigh::kNearestPacked(yCloud, k, candidateCutoff, typeI, scratch.xyz, scratch.nom);
  std::vector<double> &xyz = scratch.xyz;
  std::vector<int> &nom = scratch.nom;
  std::vector<int> &nkeep = scratch.nkeep;
  std::vector<int> &neigh = scratch.neigh;
  std::vector<double> &qlm = scratch.qlm;
  const auto tNom = std::chrono::steady_clock::now();

  const int nComp = 2 * orderL + 1;
  const int width = 2 * nComp;
  constexpr double pi = 3.14159265358979323846;
  const double prefactor = 4.0 * pi / (2.0 * orderL + 1.0);
  SteinhardtQl result;
  result.ql.assign(static_cast<std::size_t>(n), 0.0);
  result.qlBar.assign(static_cast<std::size_t>(n), 0.0);
  const gen::FracBox frame = gen::makeFracBox(yCloud);
  const bool ortho = frame.ok && !frame.triclinic;
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= kParallelThreshold)
#endif
  for (int i = 0; i < n; i++) {
    int seen[32];
    int count = 0;
    const int *row = nom.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
    const double xi = xyz[static_cast<std::size_t>(3 * i)];
    const double yi = xyz[static_cast<std::size_t>(3 * i) + 1];
    const double zi = xyz[static_cast<std::size_t>(3 * i) + 2];
    double dxs[32];
    double dys[32];
    double dzs[32];
    for (int t = 0; t < k; t++) {
      const int j = row[t];
      if (j < 0 || j >= n || j == i || !rowLists(nom.data(), k, j, i)) {
        continue;
      }
      bool dup = false;
      for (int s = 0; s < count; s++) {
        if (seen[s] == j) {
          dup = true;
          break;
        }
      }
      if (dup) {
        continue;
      }
      const double xj = xyz[static_cast<std::size_t>(3 * j)];
      const double yj = xyz[static_cast<std::size_t>(3 * j) + 1];
      const double zj = xyz[static_cast<std::size_t>(3 * j) + 2];
      double dx, dy, dz;
      if (ortho) {
        dx = xi - xj;
        dy = yi - yj;
        dz = zi - zj;
        dx -= frame.lx * std::floor(dx / frame.lx + 0.5);
        dy -= frame.ly * std::floor(dy / frame.ly + 0.5);
        dz -= frame.lz * std::floor(dz / frame.lz + 0.5);
      } else {
        const auto image = gen::fracDelta(frame, xi, yi, zi, xj, yj, zj);
        dx = image[0];
        dy = image[1];
        dz = image[2];
      }
      if (!(dx * dx + dy * dy + dz * dz > 0.0)) {
        continue;
      }
      seen[count] = j;
      dxs[count] = dx;
      dys[count] = dy;
      dzs[count] = dz;
      ++count;
    }
    nkeep[static_cast<std::size_t>(i)] = count;
    int *dstN = neigh.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
    for (int s = 0; s < count; s++) {
      dstN[s] = seen[s];
    }
    double *dst = qlm.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(width);
    // The row is reused across calls. A leftover harmonic would be added
    // to this frame.
    std::fill_n(dst, width, 0.0);
    if (count == 0) {
      continue;
    }
    double re[32 * 13];
    double im[32 * 13];
    seams::ylm::ylmBatch(orderL, dxs, dys, dzs, count, re, im);
    for (int b = 0; b < count; b++) {
      const double *yr = re + static_cast<std::size_t>(b) * static_cast<std::size_t>(nComp);
      const double *imag = im + static_cast<std::size_t>(b) * static_cast<std::size_t>(nComp);
      for (int m = 0; m < nComp; m++) {
        dst[2 * m] += yr[m];
        dst[2 * m + 1] += imag[m];
      }
    }
    const double inv = 1.0 / static_cast<double>(count);
    double sum = 0.0;
    for (int m = 0; m < nComp; m++) {
      dst[2 * m] *= inv;
      dst[2 * m + 1] *= inv;
      sum += dst[2 * m] * dst[2 * m] + dst[2 * m + 1] * dst[2 * m + 1];
    }
    result.ql[static_cast<std::size_t>(i)] = std::sqrt(prefactor * sum);
  }
  const auto tY = std::chrono::steady_clock::now();
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (n >= kParallelThreshold)
#endif
  for (int i = 0; i < n; i++) {
    const int count = nkeep[static_cast<std::size_t>(i)];
    if (count == 0) {
      result.qlBar[static_cast<std::size_t>(i)] = result.ql[static_cast<std::size_t>(i)];
      continue;
    }
    const double *self = qlm.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(width);
    const int *nb = neigh.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
    double sum = 0.0;
    const double inv = 1.0 / static_cast<double>(count + 1);
    for (int m = 0; m < nComp; m++) {
      double re = self[2 * m];
      double im = self[2 * m + 1];
      for (int s = 0; s < count; s++) {
        const double *other =
            qlm.data() + static_cast<std::size_t>(nb[s]) * static_cast<std::size_t>(width);
        re += other[2 * m];
        im += other[2 * m + 1];
      }
      re *= inv;
      im *= inv;
      sum += re * re + im * im;
    }
    result.qlBar[static_cast<std::size_t>(i)] = std::sqrt(prefactor * sum);
  }
  if (timing) {
    const auto ms = [](auto a, auto b) {
      return std::chrono::duration<double, std::milli>(b - a).count();
    };
    const auto tEnd = std::chrono::steady_clock::now();
    std::cerr << "ql-stages l=" << orderL << " n=" << n << " nom_ms=" << ms(t0, tNom)
              << " ylm_ms=" << ms(tNom, tY) << " bar_ms=" << ms(tY, tEnd) << "\n";
  }
  return result;
}

void writeSteinhardtQlm(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                        const std::vector<std::vector<int>> &nList, int orderL,
                        std::vector<double> &qlm, const std::vector<char> *which) {
  const int width = steinhardtQlmWidth(orderL);
  qlm.assign(static_cast<std::size_t>(std::max(yCloud.nop, 0)) * static_cast<std::size_t>(width),
             0.0);
  if (yCloud.nop <= 0 || (orderL != 3 && orderL != 4 && orderL != 6 && orderL != 8 &&
                          orderL != 12)) {
    return;
  }
  const NeighbourCSR graph = flatten(yCloud, nList);
  if (orderL == 12) {
    runPass1(graph, orderL, 0, graph.nop, qlm);
    if (which != nullptr) {
      for (int i = 0; i < graph.nop; i++) {
        if (static_cast<int>(which->size()) <= i || (*which)[static_cast<std::size_t>(i)] == 0) {
          std::fill_n(qlm.begin() + static_cast<std::ptrdiff_t>(i * width), width, 0.0);
        }
      }
    }
    return;
  }
  runPass1Mask(graph, orderL, qlm, which);
}

void readSteinhardtQl(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                      const std::vector<std::vector<int>> &nList, int orderL,
                      const std::vector<double> &qlm, std::vector<double> &ql,
                      std::vector<double> &qlBar, const std::vector<char> *which) {
  if (static_cast<int>(ql.size()) != yCloud.nop ||
      static_cast<int>(qlBar.size()) != yCloud.nop) {
    ql.assign(static_cast<std::size_t>(std::max(yCloud.nop, 0)), 0.0);
    qlBar.assign(ql.size(), 0.0);
  }
  if (yCloud.nop <= 0) {
    return;
  }
  const NeighbourCSR graph = flatten(yCloud, nList);
  runPass2Mask(graph, orderL, qlm, ql, qlBar, which);
}

#ifdef SEAMS_HAS_MPI
AnalysisPhases gQ6Phases;

namespace {

const std::vector<std::vector<int>> &qRows(
    const molSys::PointCloud<molSys::Point<double>, double> &cloud, double cutoff,
    const std::vector<char> &owned) {
  bool single = true;
  for (const char flag : owned) {
    if (!flag) {
      single = false;
      break;
    }
  }
  static std::vector<std::vector<int>> window;
  if (nneigh::blockWindowRows(cloud, cutoff, single ? nullptr : &owned, window)) {
    return window;
  }
  if (single || cloud.nop < 40000) {
    if (!nneigh::refillNeighListO(cutoff, cloud, 1, window)) {
      window = nneigh::neighListO(cutoff, cloud, 1);
    }
  } else {
    window = nneigh::reusedShellList(cutoff, cloud, &owned);
  }
  return window;
}

struct Q6During {
  double cutoff = 3.5;
  int orderL = 6;
  std::vector<char> interior;
  std::vector<double> qlm;
  std::vector<std::vector<int>> rows;
  NeighbourCSR graph;
  bool ran = false;
};

Q6During gQ6;

void duringQ6Halo(void *) {
  const auto t0 = std::chrono::steady_clock::now();
  auto &work = seams::domain::gatheredWork();
  bool single = true;
  for (const char flag : work.share.owned) {
    if (!flag) {
      single = false;
      break;
    }
  }
  if (!nneigh::blockWindowRows(work.cloud, gQ6.cutoff, single ? nullptr : &work.share.owned,
                               gQ6.rows)) {
    gQ6.rows = qRows(work.cloud, gQ6.cutoff, work.share.owned);
  }
  gQ6Phases.rebuildMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  const auto pass0 = std::chrono::steady_clock::now();
  gQ6.graph = flatten(work.cloud, gQ6.rows);
  const int width = steinhardtQlmWidth(gQ6.orderL);
  gQ6.qlm.assign(static_cast<std::size_t>(gQ6.graph.nop) * static_cast<std::size_t>(width), 0.0);
  runPass1Mask(gQ6.graph, gQ6.orderL, gQ6.qlm, &work.share.owned);
  gQ6Phases.computeMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - pass0).count();
  gQ6.ran = true;
}

} // namespace

SteinhardtQl distributedSteinhardt(
    const molSys::PointCloud<molSys::Point<double>, double> &full, double cutoff, int orderL,
    MPI_Comm comm) {
  setSteinhardtAtomSplit(false);
  SteinhardtQl result;
  result.ql.assign(static_cast<std::size_t>(std::max(full.nop, 0)), 0.0);
  result.qlBar.assign(result.ql.size(), 0.0);
  int nRanks = 1;
  MPI_Comm_size(comm, &nRanks);
  if (nRanks <= 1) {
    const auto t0 = std::chrono::steady_clock::now();
    thread_local std::vector<std::vector<int>> held;
    if (!nneigh::refillNeighListO(cutoff, full, 1, held)) {
      held = nneigh::neighListO(cutoff, full, 1);
    }
    const auto &rows = held;
    const auto t1 = std::chrono::steady_clock::now();
    auto out = steinhardtQl(full, rows, orderL);
    const auto t2 = std::chrono::steady_clock::now();
    gQ6Phases.rebuildMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    gQ6Phases.computeMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
    gQ6Phases.exchangeMs = 0;
    gQ6Phases.reduceMs = 0;
    return out;
  }
  gQ6Phases = {};
  gQ6.cutoff = cutoff;
  gQ6.orderL = orderL;
  gQ6.ran = false;
  seams::domain::setHaloHook(duringQ6Halo, nullptr);
  auto &work = seams::domain::gatherWork(full, cutoff, comm, seams::domain::kWorkQ6);
  seams::domain::setHaloHook(nullptr, nullptr);
  const int width = steinhardtQlmWidth(orderL);
  const std::vector<std::vector<int>> &rows =
      gQ6.ran ? gQ6.rows : qRows(work.cloud, cutoff, work.share.owned);
  if (!gQ6.ran) {
    const auto t0 = std::chrono::steady_clock::now();
    gQ6.graph = flatten(work.cloud, rows);
    gQ6Phases.rebuildMs =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  }
  const NeighbourCSR &graph = gQ6.graph;
  std::vector<double> qlm;
  if (gQ6.ran && static_cast<int>(gQ6.qlm.size()) == graph.nop * width) {
    qlm = std::move(gQ6.qlm);
  } else {
    const auto pass0 = std::chrono::steady_clock::now();
    qlm.assign(static_cast<std::size_t>(graph.nop) * static_cast<std::size_t>(width), 0.0);
    runPass1Mask(graph, orderL, qlm, &work.share.owned);
    gQ6Phases.computeMs =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - pass0).count();
  }
  const auto exchange0 = std::chrono::steady_clock::now();
  seams::domain::postBoundaryValues(qlm.data(), width, comm);
  std::vector<char> interior(static_cast<std::size_t>(work.cloud.nop), 0);
  std::vector<char> boundary(interior.size(), 0);
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (omp_get_max_threads() > 1)
#endif
  for (int i = 0; i < work.cloud.nop; i++) {
    if (!work.share.owned[static_cast<std::size_t>(i)]) {
      continue;
    }
    bool inside = true;
    const int begin = graph.offsets[static_cast<std::size_t>(i)];
    const int end = graph.offsets[static_cast<std::size_t>(i) + 1];
    for (int t = begin; t < end; t++) {
      const int j = graph.cols[static_cast<std::size_t>(t)];
      if (j < 0 || j >= work.cloud.nop || !work.share.owned[static_cast<std::size_t>(j)]) {
        inside = false;
        break;
      }
    }
    (inside ? interior : boundary)[static_cast<std::size_t>(i)] = 1;
  }
  std::vector<double> ql(static_cast<std::size_t>(graph.nop), 0.0);
  std::vector<double> qlBar(ql.size(), 0.0);
  const auto pass0 = std::chrono::steady_clock::now();
  runPass2Mask(graph, orderL, qlm, ql, qlBar, &interior);
  seams::domain::waitBoundaryValues();
  gQ6Phases.exchangeMs = seams::domain::takeValueCommMs();
  runPass2Mask(graph, orderL, qlm, ql, qlBar, &boundary);
  gQ6Phases.computeMs +=
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - pass0)
          .count();
  if (!rows.empty()) {
    std::vector<double> effort(rows.size(), 0.0);
    for (std::size_t k = 0; k < rows.size() && k < work.share.owned.size(); k++) {
      if (work.share.owned[k]) {
        effort[k] = static_cast<double>(rows[k].size());
      }
    }
    seams::domain::captureLocalWork(seams::domain::kWorkQ6, effort.data(),
                                    static_cast<int>(effort.size()), comm);
  }
  const auto red0 = std::chrono::steady_clock::now();
  for (std::size_t k = 0; k < work.share.local.size(); k++) {
    if (!work.share.owned[k]) {
      continue;
    }
    const int frame = work.share.local[k];
    result.ql[static_cast<std::size_t>(frame)] = ql[k];
    result.qlBar[static_cast<std::size_t>(frame)] = qlBar[k];
  }
  gQ6Phases.reduceMs =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - red0).count();
  return result;
}

AnalysisPhases lastSteinhardtPhases() { return gQ6Phases; }
#endif

} // namespace chill
