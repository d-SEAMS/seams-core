//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Time local Steinhardt q_l. The frame file is a little-endian header
// (uint32 n, float64 side) followed by n*3 float64 coordinates.
// One warm-up, then the repeats, one millisecond sample per line.
// The neighbour build is inside every sample. Optional last argument
// writes the warm-up q_l as raw float64.

#include <bop.hpp>
#include <mol_sys.hpp>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

namespace {

molSys::PointCloud<molSys::Point<double>, double> load(const std::string &path, double &side) {
  std::ifstream in(path, std::ios::binary);
  std::uint32_t n = 0;
  in.read(reinterpret_cast<char *>(&n), 4);
  in.read(reinterpret_cast<char *>(&side), 8);
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  in.read(reinterpret_cast<char *>(xyz.data()),
          static_cast<std::streamsize>(xyz.size() * sizeof(double)));
  if (!in) {
    std::fprintf(stderr, "short frame\n");
    std::exit(1);
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  cloud.pts.resize(n);
  for (std::uint32_t i = 0; i < n; i++) {
    auto &p = cloud.pts[i];
    p.x = xyz[3 * i];
    p.y = xyz[3 * i + 1];
    p.z = xyz[3 * i + 2];
    p.type = 1;
    p.atomID = static_cast<int>(i + 1);
    p.molID = p.atomID;
  }
  cloud.nop = static_cast<int>(n);
  return cloud;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: bench_ql frame.bin order repeats [ql.bin]\n");
    return 2;
  }
  double side = 0;
  auto cloud = load(argv[1], side);
  const int order = std::atoi(argv[2]);
  const int repeats = std::atoi(argv[3]);
  auto warm = chill::steinhardtQlCutoff(cloud, 3.5, order, 1, false);
  if (warm.ql.empty()) {
    std::fprintf(stderr, "empty q_l\n");
    return 1;
  }
  double sink = 0;
  for (double v : warm.ql) {
    sink += v;
  }
  for (int r = 0; r < repeats; r++) {
    const auto t0 = std::chrono::steady_clock::now();
    auto got = chill::steinhardtQlCutoff(cloud, 3.5, order, 1, false);
    const double ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    sink += got.ql[0];
    std::printf("%.6f\n", ms);
  }
  if (argc > 4) {
    std::ofstream out(argv[4], std::ios::binary);
    out.write(reinterpret_cast<const char *>(warm.ql.data()),
              static_cast<std::streamsize>(warm.ql.size() * sizeof(double)));
  }
  std::fprintf(stderr, "n=%d l=%d side=%.6f q0=%.12f sum=%.8f\n", cloud.nop, order, side,
               warm.ql[0], sink);
  return 0;
}
