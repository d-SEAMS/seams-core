//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <neighbours.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

molSys::PointCloud<molSys::Point<double>, double> iceIc(int n, double a) {
  const double b[8][3] = {
      {0, 0, 0},          {0, 0.5, 0.5},      {0.5, 0, 0.5},
      {0.5, 0.5, 0},      {0.25, 0.25, 0.25}, {0.25, 0.75, 0.75},
      {0.75, 0.25, 0.75}, {0.75, 0.75, 0.25}};
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {n * a, n * a, n * a};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  cloud.nop = n * n * n * 8;
  cloud.pts.resize(static_cast<std::size_t>(cloud.nop));
  int t = 0;
  for (int ix = 0; ix < n; ++ix) {
    for (int iy = 0; iy < n; ++iy) {
      for (int iz = 0; iz < n; ++iz) {
        for (int k = 0; k < 8; ++k, ++t) {
          auto &p = cloud.pts[static_cast<std::size_t>(t)];
          p.type = 1;
          p.atomID = t;
          p.molID = t;
          p.x = (ix + b[k][0]) * a;
          p.y = (iy + b[k][1]) * a;
          p.z = (iz + b[k][2]) * a;
          cloud.idIndexMap[t] = t;
        }
      }
    }
  }
  return cloud;
}

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

std::uint64_t checksum(const std::vector<int> &offsets, const std::vector<int> &cols) {
  std::uint64_t h = 1469598103934665603ull;
  for (std::size_t i = 0; i + 1 < offsets.size(); ++i) {
    for (int k = offsets[i]; k < offsets[i + 1]; ++k) {
      h ^= static_cast<std::uint64_t>(i) * 1315423911ull +
           static_cast<std::uint64_t>(cols[static_cast<std::size_t>(k)]);
      h *= 1099511628211ull;
    }
  }
  return h;
}

} // namespace

int main(int argc, char **argv) {
  int molecules = 110592;
  int threads = 1;
  int reps = 8;
  if (argc > 1) {
    molecules = std::atoi(argv[1]);
  }
  if (argc > 2) {
    threads = std::atoi(argv[2]);
  }
  if (argc > 3) {
    reps = std::atoi(argv[3]);
  }
  if (threads < 1) {
    threads = 1;
  }
  if (reps < 1) {
    reps = 1;
  }
  int n = 1;
  while (n * n * n * 8 < molecules) {
    ++n;
  }
  if (n * n * n * 8 != molecules) {
    std::fprintf(stderr, "molecules must be 8 n^3\n");
    return 2;
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(threads);
#endif
  auto cloud = iceIc(n, 6.35);
  const double cutoff = 3.5;
  std::vector<int> offsets, cols;
  nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols);
  const int pairs = offsets.empty() ? 0 : offsets.back();
  const auto sum = checksum(offsets, cols);
  std::vector<double> csr(static_cast<std::size_t>(reps));
  std::vector<double> rows(static_cast<std::size_t>(reps));
  for (int r = 0; r < reps; ++r) {
    auto t0 = std::chrono::steady_clock::now();
    nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols);
    auto t1 = std::chrono::steady_clock::now();
    auto listed = nneigh::neighListO(cutoff, cloud, 1);
    auto t2 = std::chrono::steady_clock::now();
    csr[static_cast<std::size_t>(r)] =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    rows[static_cast<std::size_t>(r)] =
        std::chrono::duration<double, std::milli>(t2 - t1).count();
    if (listed.size() != static_cast<std::size_t>(cloud.nop) || pairs != 4 * cloud.nop) {
      std::fprintf(stderr, "pair count %d want %d\n", pairs, 4 * cloud.nop);
      return 1;
    }
  }
  std::printf("molecules=%d threads=%d pairs=%d checksum=%llu csr_median_ms=%.4f "
              "rows_median_ms=%.4f csr_samples=",
              molecules, threads, pairs, static_cast<unsigned long long>(sum),
              median(csr), median(rows));
  for (int r = 0; r < reps; ++r) {
    std::printf("%s%.4f", r ? "," : "", csr[static_cast<std::size_t>(r)]);
  }
  std::printf(" rows_samples=");
  for (int r = 0; r < reps; ++r) {
    std::printf("%s%.4f", r ? "," : "", rows[static_cast<std::size_t>(r)]);
  }
  std::printf("\n");
  return 0;
}
