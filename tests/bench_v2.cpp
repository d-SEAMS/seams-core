//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

#include <bop.hpp>
#include <diamond_cna.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

molSys::PointCloud<molSys::Point<double>, double> ice(int n, bool tilt) {
  const double a = 6.35;
  const double L = n * a;
  const double b[8][3] = {
      {0, 0, 0},          {0, 0.5, 0.5},      {0.5, 0, 0.5},
      {0.5, 0.5, 0},      {0.25, 0.25, 0.25}, {0.25, 0.75, 0.75},
      {0.75, 0.25, 0.75}, {0.75, 0.75, 0.25}};
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = n * n * n * 8;
  cloud.currentFrame = 1;
  cloud.pts.resize(static_cast<std::size_t>(cloud.nop));
  const double xy = tilt ? 0.25 * L : 0.0;
  const double xz = tilt ? 0.05 * L : 0.0;
  const double yz = tilt ? 0.12 * L : 0.0;
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  if (tilt) {
    cloud.box = {L + xmax - xmin, L + ymax - ymin, L, xy, xz, yz};
    cloud.boxLow = {xmin, ymin, 0.0};
  } else {
    cloud.box = {L, L, L};
    cloud.boxLow = {0.0, 0.0, 0.0};
  }
  int t = 0;
  for (int ix = 0; ix < n; ++ix) {
    for (int iy = 0; iy < n; ++iy) {
      for (int iz = 0; iz < n; ++iz) {
        for (int k = 0; k < 8; ++k, ++t) {
          const double s0 = (ix + b[k][0]) / n;
          const double s1 = (iy + b[k][1]) / n;
          const double s2 = (iz + b[k][2]) / n;
          auto &p = cloud.pts[static_cast<std::size_t>(t)];
          p.type = 1;
          p.atomID = t;
          p.molID = t;
          p.x = xmin + L * s0 + xy * s1 + xz * s2;
          p.y = ymin + L * s1 + yz * s2;
          p.z = s2 * L;
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
  if (argc < 5) {
    std::fprintf(stderr, "usage: bench_v2 tilt|chill|rings|symmetry molecules threads reps\n");
    return 2;
  }
  const std::string mode = argv[1];
  const int molecules = std::atoi(argv[2]);
  int threads = std::atoi(argv[3]);
  int reps = std::atoi(argv[4]);
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
  const bool tilt = mode == "tilt";
  auto cloud = ice(n, tilt);
  const double cutoff = 3.5;
  std::vector<double> samples(static_cast<std::size_t>(reps));
  int pairs = 0;
  int sink = 0;
  std::uint64_t sum = 0;
  for (int r = 0; r < reps + 1; ++r) {
    auto t0 = std::chrono::steady_clock::now();
    std::vector<int> offsets, cols;
    if (!nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols, false)) {
      std::fprintf(stderr, "csr failed\n");
      return 1;
    }
    pairs = offsets.back();
    sum = checksum(offsets, cols);
    if (mode == "chill") {
      sink = chill::chillPlusFromCsr(cloud, offsets, cols);
    } else if (mode == "rings") {
      const auto rings = primitive::ringNetwork(offsets, cols, 6);
      sink = static_cast<int>(rings.size());
    } else if (mode == "symmetry") {
      std::vector<int> kind;
      symm::diamondFromCsr(cloud, offsets, cols, kind);
      int labeled = 0;
      for (int v : kind) {
        labeled += v > 0 ? 1 : 0;
      }
      sink = labeled;
    }
    auto t1 = std::chrono::steady_clock::now();
    if (r == 0) {
      continue;
    }
    samples[static_cast<std::size_t>(r - 1)] =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
  }
  std::printf("mode=%s molecules=%d threads=%d pairs=%d checksum=%llu sink=%d "
              "median_ms=%.4f samples_ms=",
              mode.c_str(), molecules, threads, pairs,
              static_cast<unsigned long long>(sum), sink, median(samples));
  for (int r = 0; r < reps; ++r) {
    std::printf("%s%.4f", r ? "," : "", samples[static_cast<std::size_t>(r)]);
  }
  std::printf("\n");
  return 0;
}
