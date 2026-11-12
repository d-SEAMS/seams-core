//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Ring topology plus a tolerant shape fit, for a thermal hydrate and the
// cages of its steepest-descent inherent structure.

#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <local_symmetry.hpp>

#include <cmath>
#include <cstdlib>
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

std::vector<std::vector<int>> bonds(const Eigen::MatrixXd &xyz, const Eigen::Matrix3d &box,
                                    double cutoff) {
  const int n = static_cast<int>(xyz.rows());
  std::vector<std::vector<int>> list(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    list[static_cast<size_t>(i)].push_back(i);
  }
  const double cut2 = cutoff * cutoff;
  const double lx = box(0, 0);
  const double ly = box(1, 1);
  const double lz = box(2, 2);
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double dx = xyz(j, 0) - xyz(i, 0);
      double dy = xyz(j, 1) - xyz(i, 1);
      double dz = xyz(j, 2) - xyz(i, 2);
      dx -= lx * std::nearbyint(dx / lx);
      dy -= ly * std::nearbyint(dy / ly);
      dz -= lz * std::nearbyint(dz / lz);
      if (dx * dx + dy * dy + dz * dz <= cut2) {
        list[static_cast<size_t>(i)].push_back(j);
        list[static_cast<size_t>(j)].push_back(i);
      }
    }
  }
  return list;
}

struct CageHit {
  int kind = 0;
  double cx = 0;
  double cy = 0;
  double cz = 0;
  double rmsd = 0;
};

std::vector<CageHit> cagesOf(const Eigen::MatrixXd &xyz, const Eigen::Matrix3d &box, double cutoff,
                             double shapeCeiling) {
  const auto list = bonds(xyz, box, cutoff);
  const auto rings = primitive::ringNetwork(list, 6);
  const cage::Signature sigs[] = {cage::Signature::parse("5:12"), cage::Signature::parse("5:12,6:2"),
                                  cage::Signature::parse("5:12,6:4")};
  const symm::Kind kinds[] = {symm::Kind::cage512, symm::Kind::cage51262, symm::Kind::cage51264};
  std::vector<CageHit> hits;
  for (int s = 0; s < 3; s++) {
    const auto found = cage::findBySignature(rings, sigs[s]);
    for (const auto &cage : found) {
      if (cage.vertices.empty()) {
        continue;
      }
      Eigen::MatrixXd cluster(static_cast<Eigen::Index>(cage.vertices.size()), 3);
      double cx = 0;
      double cy = 0;
      double cz = 0;
      for (size_t i = 0; i < cage.vertices.size(); i++) {
        const int v = cage.vertices[i];
        cluster.row(static_cast<Eigen::Index>(i)) = xyz.row(v);
        cx += xyz(v, 0);
        cy += xyz(v, 1);
        cz += xyz(v, 2);
      }
      const double inv = 1.0 / static_cast<double>(cage.vertices.size());
      const auto hit = symm::matchKernel(cluster);
      if (hit.rmsd > shapeCeiling) {
        continue;
      }
      CageHit row;
      row.kind = static_cast<int>(kinds[s]);
      row.cx = cx * inv;
      row.cy = cy * inv;
      row.cz = cz * inv;
      row.rmsd = hit.rmsd;
      hits.push_back(row);
    }
  }
  return hits;
}

int matched(const std::vector<CageHit> &thermal, const std::vector<CageHit> &inherent, int kind,
            double box) {
  int n = 0;
  std::vector<int> used(inherent.size(), 0);
  for (const auto &cage : thermal) {
    if (cage.kind != kind) {
      continue;
    }
    double best = 1e300;
    int at = -1;
    for (size_t i = 0; i < inherent.size(); i++) {
      if (used[i] || inherent[i].kind != kind) {
        continue;
      }
      double dx = cage.cx - inherent[i].cx;
      double dy = cage.cy - inherent[i].cy;
      double dz = cage.cz - inherent[i].cz;
      dx -= box * std::nearbyint(dx / box);
      dy -= box * std::nearbyint(dy / box);
      dz -= box * std::nearbyint(dz / box);
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 < best) {
        best = d2;
        at = static_cast<int>(i);
      }
    }
    if (at >= 0 && best < 2.0 * 2.0) {
      used[static_cast<size_t>(at)] = 1;
      n++;
    }
  }
  return n;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: report_cages thermal.frame quench.frame [cutoff]\n";
    return 2;
  }
  const double cutoff = argc > 3 ? std::atof(argv[3]) : 3.5;
  Eigen::Matrix3d box;
  Eigen::Matrix3d qbox;
  Eigen::MatrixXd xyz;
  Eigen::MatrixXd qxyz;
  if (!load(argv[1], box, xyz) || !load(argv[2], qbox, qxyz)) {
    std::cerr << "unreadable\n";
    return 1;
  }
  const double shape = 1.5;
  const auto thermal = cagesOf(xyz, box, cutoff, shape);
  const auto inherent = cagesOf(qxyz, qbox, 3.3, shape);
  const char *names[] = {"512", "51262", "51264"};
  const int kinds[] = {static_cast<int>(symm::Kind::cage512), static_cast<int>(symm::Kind::cage51262),
                       static_cast<int>(symm::Kind::cage51264)};
  std::cout << argv[1] << " cutoff=" << cutoff << " shape<" << shape << "\n";
  for (int s = 0; s < 3; s++) {
    int nT = 0;
    int nQ = 0;
    for (const auto &cage : thermal) {
      if (cage.kind == kinds[s]) {
        nT++;
      }
    }
    for (const auto &cage : inherent) {
      if (cage.kind == kinds[s]) {
        nQ++;
      }
    }
    const int agree = matched(thermal, inherent, kinds[s], box(0, 0));
    std::cout << "  " << names[s] << " thermal=" << nT << " inherent=" << nQ << " matched=" << agree
              << "\n";
  }
  return 0;
}
