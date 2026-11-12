//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Three per-atom decisions on a thermal frame, scored against the locked
// second shell. The lock is a bond-graph residual below 0.45 A. The hard
// count thresholds the bond cosines at 0.75. The soft count picks the
// nearest ideal cosine pattern. The template fit is classifyIce.

#include <local_symmetry.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Frame {
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  Eigen::MatrixXd xyz;
};

bool loadFrame(const std::string &path, Frame &frame) {
  std::ifstream in(path);
  if (!in) {
    return false;
  }
  std::string tag;
  std::string name;
  std::string role;
  double temperature = 0;
  int n = 0;
  in >> tag >> name >> tag >> temperature >> tag >> role >> tag >> n >> tag;
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      in >> frame.box(r, c);
    }
  }
  in >> tag;
  frame.xyz.resize(n, 3);
  for (int i = 0; i < n; i++) {
    int layer = 0;
    int core = 0;
    in >> frame.xyz(i, 0) >> frame.xyz(i, 1) >> frame.xyz(i, 2) >> layer >> core;
  }
  return static_cast<bool>(in);
}

int hardCount(const double c[4]) {
  int n = 0;
  int ecl = 0;
  for (int i = 0; i < 4; i++) {
    if (c[i] < -1.0) {
      continue;
    }
    n++;
    if (c[i] > 0.75) {
      ecl++;
    }
  }
  if (n < 4 || ecl > 2) {
    return -1;
  }
  if (ecl == 0) {
    return 1;
  }
  if (ecl == 1) {
    return 0;
  }
  return 2;
}

int softCount(const double c[4]) {
  const double proto[3][4] = {{1.0, 0.5, 0.5, 0.5}, {0.5, 0.5, 0.5, 0.5}, {1.0, 1.0, 0.5, 0.5}};
  const int kindOf[3] = {0, 1, 2};
  double best = 1e300;
  int bestK = -1;
  for (int k = 0; k < 3; k++) {
    double sum = 0.0;
    for (int i = 0; i < 4; i++) {
      if (c[i] < -1.0) {
        return -1;
      }
      const double d = c[i] - proto[k][i];
      sum += d * d;
    }
    if (sum < best) {
      best = sum;
      bestK = kindOf[k];
    }
  }
  if (!(best < 0.20)) {
    return -1;
  }
  return bestK;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: bench_robust frame...\n";
    return 2;
  }
  for (int a = 1; a < argc; a++) {
    Frame frame;
    if (!loadFrame(argv[a], frame)) {
      std::cerr << "unreadable " << argv[a] << "\n";
      return 1;
    }
    const int n = static_cast<int>(frame.xyz.rows());
    const double cuts[] = {0.55, 0.75, 0.95, 1.10, 1.30, 1.50, 2.00};
    std::vector<std::vector<int>> byCut(7);
    std::vector<int> timed;
    const auto t0 = std::chrono::steady_clock::now();
    symm::classifyIce(frame.xyz, frame.box, timed, 0.55);
    const double ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    (void)timed;
    for (int c = 0; c < 7; c++) {
      symm::classifyIce(frame.xyz, frame.box, byCut[static_cast<size_t>(c)], cuts[c]);
    }
    int cutOk[7] = {};
    int falsePos[7] = {};
    int conf[3][3] = {};
    int locked = 0;
    std::ofstream gtOut(std::string(argv[a]) + ".gt");
    int hardOk = 0;
    int softOk = 0;
    int fitOk = 0;
    int gtIh = 0;
    int gtIc = 0;
    int gtIsd = 0;
    for (int i = 0; i < n; i++) {
      const Eigen::MatrixXd shell = symm::iceShell(frame.xyz, frame.box, i);
      if (shell.rows() != 17) {
        continue;
      }
      const double r[3] = {symm::iceTemplateRmsd(shell, symm::Kind::iceIh),
                           symm::iceTemplateRmsd(shell, symm::Kind::iceIc),
                           symm::iceTemplateRmsd(shell, symm::Kind::iceIsd)};
      int gt = -1;
      double best = 1e300;
      for (int k = 0; k < 3; k++) {
        if (r[k] >= 0.0 && r[k] < best) {
          best = r[k];
          gt = k;
        }
      }
      gtOut << i << " " << (best < 0.45 ? gt : -1) << "\n";
      if (!(best < 0.45)) {
        for (int c = 0; c < 7; c++) {
          if (byCut[static_cast<size_t>(c)][static_cast<size_t>(i)] >= 0) {
            falsePos[c]++;
          }
        }
        continue;
      }
      locked++;
      if (gt == 0) {
        gtIh++;
      } else if (gt == 1) {
        gtIc++;
      } else {
        gtIsd++;
      }
      double cosine[4];
      symm::bondCosines(shell, cosine);
      if (hardCount(cosine) == gt) {
        hardOk++;
      }
      if (softCount(cosine) == gt) {
        softOk++;
      }
      for (int c = 0; c < 7; c++) {
        const int pred = byCut[static_cast<size_t>(c)][static_cast<size_t>(i)];
        if (pred == gt) {
          cutOk[c]++;
        }
        if (c == 0 && pred >= 0 && pred < 3) {
          conf[gt][pred]++;
        }
      }
    }
    std::cout << argv[a] << " n=" << n << " locked=" << locked << " gt Ih=" << gtIh
              << " Ic=" << gtIc << " Isd=" << gtIsd << " soft=" << softOk
              << " classify_ms=" << ms;
    for (int c = 0; c < 7; c++) {
      std::cout << " c" << cuts[c] << "=" << cutOk[c] << "/" << falsePos[c];
    }
    std::cout << " confusion";
    for (int g = 0; g < 3; g++) {
      std::cout << " [" << conf[g][0] << " " << conf[g][1] << " " << conf[g][2] << "]";
    }
    std::cout << "\n";
    (void)hardOk;
    (void)fitOk;
  }
  return 0;
}
