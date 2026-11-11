// Time classifyIce on ideal ice Ic. Release build of the symmetry kernel.
// One warm-up, six calls. Prints one line of samples and the label census.

#include <local_symmetry.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

Eigen::MatrixXd iceIc(int repetitions, double bond, double &side) {
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double basis[8][3] = {
      {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75},
  };
  const int nCell = repetitions * repetitions * repetitions;
  const int n = nCell * 8;
  Eigen::MatrixXd points(n, 3);
  int row = 0;
  for (int i = 0; i < repetitions; ++i) {
    for (int j = 0; j < repetitions; ++j) {
      for (int k = 0; k < repetitions; ++k) {
        for (const auto &b : basis) {
          points(row, 0) = (static_cast<double>(i) + b[0]) * lattice;
          points(row, 1) = (static_cast<double>(j) + b[1]) * lattice;
          points(row, 2) = (static_cast<double>(k) + b[2]) * lattice;
          ++row;
        }
      }
    }
  }
  side = static_cast<double>(repetitions) * lattice;
  return points;
}

} // namespace

std::string cpusAllowed() {
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    const std::string key = "Cpus_allowed_list:";
    if (line.rfind(key, 0) == 0) {
      const auto pos = line.find(':');
      std::string value = line.substr(pos + 1);
      const auto begin = value.find_first_not_of(" \t");
      return begin == std::string::npos ? value : value.substr(begin);
    }
  }
  return "unknown";
}

int main(int argc, char **argv) {
  const std::string cpus = cpusAllowed();
  const int repetitions = argc > 1 ? std::atoi(argv[1]) : 24;
  const int warmup = argc > 2 ? std::atoi(argv[2]) : 1;
  const int repeats = argc > 3 ? std::atoi(argv[3]) : 6;
  double side = 0;
  const Eigen::MatrixXd points = iceIc(repetitions, 2.75, side);
  Eigen::Matrix3d vectors = Eigen::Matrix3d::Zero();
  vectors(0, 0) = side;
  vectors(1, 1) = side;
  vectors(2, 2) = side;
  std::vector<int> kind;
  for (int i = 0; i < warmup; ++i) {
    symm::classifyIce(points, vectors, kind);
  }
  std::vector<double> samples;
  samples.reserve(static_cast<size_t>(repeats));
  for (int i = 0; i < repeats; ++i) {
    const auto t0 = std::chrono::steady_clock::now();
    symm::classifyIce(points, vectors, kind);
    const auto t1 = std::chrono::steady_clock::now();
    samples.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
  }
  int census[4] = {0, 0, 0, 0};
  for (int label : kind) {
    if (label >= 0 && label <= 2) {
      census[label] += 1;
    } else {
      census[3] += 1;
    }
  }
  std::cout << "cpus=" << cpus << " molecules=" << points.rows() << " side=" << side
            << " Ih=" << census[0] << " Ic=" << census[1] << " Isd=" << census[2]
            << " other=" << census[3] << " samples_ms=";
  for (double sample : samples) {
    std::cout << sample << " ";
  }
  std::cout << "\n";
  return census[1] == static_cast<int>(points.rows()) ? 0 : 2;
}
