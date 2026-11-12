// Copyright (c) d-SEAMS developers
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

double seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
      .count();
}

double type7(std::vector<double> samples, double p) {
  std::sort(samples.begin(), samples.end());
  const double h = (static_cast<double>(samples.size()) - 1.0) * p;
  const auto lo = static_cast<std::size_t>(std::floor(h));
  const auto hi = static_cast<std::size_t>(std::ceil(h));
  if (lo == hi || hi >= samples.size()) {
    return samples[lo];
  }
  const double w = h - static_cast<double>(lo);
  return samples[lo] * (1.0 - w) + samples[hi] * w;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: bench_chemfiles_read DUMP\n");
    return 2;
  }
  const std::string path = argv[1];
  auto timeNative = [&]() {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    const auto start = std::chrono::steady_clock::now();
    cloud = sinp::readLammpsTrjO(path, 1, cloud, 1);
    if (cloud.nop <= 0) {
      cloud = sinp::readLammpsTrj(path, 1, cloud);
    }
    return seconds(start);
  };
  auto timeChem = [&]() {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    const auto start = std::chrono::steady_clock::now();
    cloud = sinp::readChemfiles(path, 1, cloud, -1);
    return seconds(start);
  };
  const double warmNative = timeNative();
  const double warmChem = timeChem();
  std::vector<double> native;
  std::vector<double> chem;
  for (int i = 0; i < 6; ++i) {
    native.push_back(timeNative());
    chem.push_back(timeChem());
  }
  std::printf("warmup_native %.17g\n", warmNative);
  std::printf("warmup_chemfiles %.17g\n", warmChem);
  std::printf("native_q1 %.17g\nnative_q3 %.17g\n", type7(native, 0.25),
              type7(native, 0.75));
  std::printf("native_ends %.17g %.17g\n", type7(native, 0.0), type7(native, 1.0));
  std::printf("chemfiles_q1 %.17g\nchemfiles_q3 %.17g\n", type7(chem, 0.25),
              type7(chem, 0.75));
  std::printf("native_samples");
  for (double sample : native) {
    std::printf(" %.17g", sample);
  }
  std::printf("\nchemfiles_samples");
  for (double sample : chem) {
    std::printf(" %.17g", sample);
  }
  std::printf("\n");
  return 0;
}
