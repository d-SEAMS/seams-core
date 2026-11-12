//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Precision and recall of the scale-free shell label against the
// steepest-descent inherent structure of the same frame.

#include <local_symmetry.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool load(const std::string &path, Eigen::Matrix3d &box, Eigen::MatrixXd &xyz) {
  std::ifstream in(path);
  if (!in) {
    return false;
  }
  std::string tag;
  int n = 0;
  double temperature = 0.0;
  in >> tag >> tag >> tag >> temperature >> tag >> tag >> tag >> n >> tag;
  box.setZero();
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      in >> box(r, c);
    }
  }
  in >> tag;
  xyz.resize(n, 3);
  int bilayer = 0;
  int core = 0;
  for (int i = 0; i < n; i++) {
    in >> xyz(i, 0) >> xyz(i, 1) >> xyz(i, 2) >> bilayer >> core;
  }
  return static_cast<bool>(in);
}

void score(const std::vector<int> &pred, const std::vector<int> &truth, int &tp, int &fp, int &fn) {
  tp = 0;
  fp = 0;
  fn = 0;
  const int n = static_cast<int>(pred.size());
  for (int i = 0; i < n; i++) {
    const bool ice = truth[static_cast<size_t>(i)] >= 0;
    const bool hit = pred[static_cast<size_t>(i)] == truth[static_cast<size_t>(i)];
    if (ice && hit) {
      tp++;
    } else if (pred[static_cast<size_t>(i)] >= 0) {
      fp++;
    } else if (ice) {
      fn++;
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: bench_roc thermal.frame quench.frame\n";
    return 2;
  }
  Eigen::Matrix3d box;
  Eigen::Matrix3d qbox;
  Eigen::MatrixXd xyz;
  Eigen::MatrixXd qxyz;
  if (!load(argv[1], box, xyz) || !load(argv[2], qbox, qxyz) || xyz.rows() != qxyz.rows()) {
    std::cerr << "unreadable pair\n";
    return 1;
  }
  // The quench sits on the lattice. A tight ceiling is the ground truth.
  std::vector<int> truth;
  symm::classifyIce(qxyz, qbox, truth, 0.06);
  const double gates[] = {0.06, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22};
  std::cout << argv[1] << " n=" << xyz.rows() << "\n";
  int gtIce = 0;
  for (int label : truth) {
    if (label >= 0) {
      gtIce++;
    }
  }
  std::cout << "  quench_ice " << gtIce << "\n";
  for (double gate : gates) {
    std::vector<int> pred;
    symm::classifyIce(xyz, box, pred, gate);
    int tp = 0;
    int fp = 0;
    int fn = 0;
    score(pred, truth, tp, fp, fn);
    const double precision = (tp + fp) > 0 ? static_cast<double>(tp) / static_cast<double>(tp + fp) : 1.0;
    const double recall = (tp + fn) > 0 ? static_cast<double>(tp) / static_cast<double>(tp + fn) : 1.0;
    std::cout << "  g=" << gate << " tp=" << tp << " fp=" << fp << " fn=" << fn
              << " precision=" << precision << " recall=" << recall << "\n";
  }
  return 0;
}
