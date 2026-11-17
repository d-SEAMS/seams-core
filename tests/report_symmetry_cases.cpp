//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

// Counts for the cells PTM does not name: the stacking fault, and the
// sI and sII oxygen lattices. Cage totals use the primitive rings.

#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <local_symmetry.hpp>
#include <seams_c_api.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "hydrate_sites.inc"

namespace {

std::vector<std::vector<int>> bonds(const Eigen::MatrixXd &xyz, double box, double cutoff) {
  const int n = static_cast<int>(xyz.rows());
  std::vector<std::vector<int>> list(static_cast<size_t>(n));
  for (int i = 0; i < n; i++) {
    // ringNetwork treats the first entry of each row as the vertex itself.
    list[static_cast<size_t>(i)].push_back(i);
  }
  const double cut2 = cutoff * cutoff;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double d2 = 0.0;
      for (int k = 0; k < 3; k++) {
        double d = xyz(j, k) - xyz(i, k);
        d -= box * std::nearbyint(d / box);
        d2 += d * d;
      }
      if (d2 <= cut2) {
        list[static_cast<size_t>(i)].push_back(j);
        list[static_cast<size_t>(j)].push_back(i);
      }
    }
  }
  return list;
}

void chillHist(const Eigen::MatrixXd &xyz, double box, const char *tag) {
  const int n = static_cast<int>(xyz.rows());
  std::vector<double> raw(static_cast<size_t>(3 * n));
  for (int i = 0; i < n; i++) {
    raw[static_cast<size_t>(3 * i)] = xyz(i, 0);
    raw[static_cast<size_t>(3 * i + 1)] = xyz(i, 1);
    raw[static_cast<size_t>(3 * i + 2)] = xyz(i, 2);
  }
  const double cell[6] = {0, 0, 0, box, box, box};
  std::vector<int> labels(static_cast<size_t>(n), -1);
  const char *names[] = {"cubic", "hexagonal", "water", "interfacial", "clathrate",
                         "interClathrate", "unclassified", "reCubic", "reHex"};
  if (seams_chill_plus_flat(raw.data(), n, cell, labels.data()) != 0) {
    std::cout << tag << " chill failed\n";
    return;
  }
  std::map<std::string, int> hist;
  for (int v : labels) {
    const char *name = (v >= 0 && v <= 8) ? names[v] : "bad";
    hist[name] += 1;
  }
  std::cout << tag << " chill";
  for (const auto &item : hist) {
    std::cout << " " << item.first << "=" << item.second;
  }
  std::cout << "\n";
}

void iceHist(const Eigen::MatrixXd &xyz, double box, const char *tag) {
  Eigen::Matrix3d vectors = Eigen::Matrix3d::Identity() * box;
  std::vector<int> kind;
  symm::classifyIce(xyz, vectors, kind);
  int nIc = 0;
  int nIh = 0;
  int nIsd = 0;
  int nNo = 0;
  for (int k : kind) {
    if (k == static_cast<int>(symm::Kind::iceIc)) {
      nIc++;
    } else if (k == static_cast<int>(symm::Kind::iceIh)) {
      nIh++;
    } else if (k == static_cast<int>(symm::Kind::iceIsd)) {
      nIsd++;
    } else {
      nNo++;
    }
  }
  std::cout << tag << " kernel ice-Ic=" << nIc << " ice-Ih=" << nIh << " ice-Isd=" << nIsd
            << " other=" << nNo << "\n";
}

void referenceShells(const Eigen::MatrixXd &xyz, const Eigen::Matrix3d &vectors, const char *tag) {
  int nIc = 0;
  int nIh = 0;
  int nIsd = 0;
  int nNo = 0;
  for (int i = 0; i < xyz.rows(); i++) {
    const auto shell = symm::iceShell(xyz, vectors, i);
    if (shell.rows() != 17) {
      nNo++;
      continue;
    }
    const auto hit = symm::best(shell, "ira");
    if (hit.status != 0) {
      nNo++;
    } else if (hit.kind == symm::Kind::iceIc) {
      nIc++;
    } else if (hit.kind == symm::Kind::iceIh) {
      nIh++;
    } else if (hit.kind == symm::Kind::iceIsd) {
      nIsd++;
    } else {
      nNo++;
    }
  }
  std::cout << tag << " ira-shell ice-Ic=" << nIc << " ice-Ih=" << nIh << " ice-Isd=" << nIsd
            << " other=" << nNo << "\n";
}

void cagesOf(const Eigen::MatrixXd &xyz, double box, const char *tag) {
  const auto list = bonds(xyz, box, 3.3);
  const auto rings = primitive::ringNetwork(list, 6);
  int n5 = 0;
  int n6 = 0;
  for (const auto &ring : rings) {
    if (ring.size() == 5) {
      n5++;
    }
    if (ring.size() == 6) {
      n6++;
    }
  }
  const auto small = cage::findBySignature(rings, cage::Signature::parse("5:12"));
  const auto sI = cage::findBySignature(rings, cage::Signature::parse("5:12,6:2"));
  const auto sII = cage::findBySignature(rings, cage::Signature::parse("5:12,6:4"));
  std::cout << tag << " rings5=" << n5 << " rings6=" << n6 << " cage512=" << small.size()
            << " cage51262=" << sI.size() << " cage51264=" << sII.size() << "\n";
  auto nameOne = [&](const std::vector<int> &verts, const char *expect) {
    if (verts.empty()) {
      return;
    }
    Eigen::MatrixXd cluster(static_cast<Eigen::Index>(verts.size()), 3);
    for (size_t i = 0; i < verts.size(); i++) {
      cluster.row(static_cast<Eigen::Index>(i)) = xyz.row(verts[i]);
    }
    const auto hit = symm::matchKernel(cluster);
    std::cout << tag << " template " << expect << " -> " << hit.name << " rmsd " << hit.rmsd
              << "\n";
  };
  if (!small.empty()) {
    nameOne(small.front().vertices, "5^12");
  }
  if (!sI.empty()) {
    nameOne(sI.front().vertices, "5^12 6^2");
  }
  if (!sII.empty()) {
    nameOne(sII.front().vertices, "5^12 6^4");
  }
}

Eigen::MatrixXd fromTable(const double table[][3], int n) {
  Eigen::MatrixXd xyz(n, 3);
  for (int i = 0; i < n; i++) {
    xyz(i, 0) = table[i][0];
    xyz(i, 1) = table[i][1];
    xyz(i, 2) = table[i][2];
  }
  return xyz;
}

} // namespace

int main() {
  const auto fault = symm::crystal(symm::Crystal::fault);
  std::cout << "fault atoms " << fault.cartesian.rows() << "\n";
  referenceShells(fault.cartesian, fault.vectors, "fault");
  // The fault cell is not orthorhombic. CHILL+ below uses the bounding cube
  // only as a large box around a shifted copy, so the periodic images do not
  // re-enter the first shell.
  Eigen::MatrixXd shifted = fault.cartesian;
  for (int i = 0; i < shifted.rows(); i++) {
    shifted.row(i) += Eigen::RowVector3d(30, 30, 30);
  }
  chillHist(shifted, 80.0, "fault-cluster");

  const auto sI = fromTable(kSIxyz, kSI);
  const auto sII = fromTable(kSIIxyz, kSII);
  iceHist(sI, kSIcell, "sI");
  chillHist(sI, kSIcell, "sI");
  cagesOf(sI, kSIcell, "sI");
  iceHist(sII, kSIIcell, "sII");
  chillHist(sII, kSIIcell, "sII");
  cagesOf(sII, kSIIcell, "sII");

  // A small ice Ic crystallite in a box of uncorrelated points.
  const auto ic = symm::crystal(symm::Crystal::iceIc);
  const int nrep = 2;
  const int base = static_cast<int>(ic.cartesian.rows());
  const int nCrystal = base * nrep * nrep * nrep;
  const int nLiquid = 128;
  Eigen::MatrixXd mix(nCrystal + nLiquid, 3);
  int row = 0;
  for (int ix = 0; ix < nrep; ix++) {
    for (int iy = 0; iy < nrep; iy++) {
      for (int iz = 0; iz < nrep; iz++) {
        const Eigen::Vector3d shift =
            static_cast<double>(ix) * ic.vectors.row(0).transpose() +
            static_cast<double>(iy) * ic.vectors.row(1).transpose() +
            static_cast<double>(iz) * ic.vectors.row(2).transpose();
        for (int i = 0; i < base; i++) {
          mix.row(row++) = ic.cartesian.row(i) + shift.transpose();
        }
      }
    }
  }
  for (int i = 0; i < nLiquid; i++) {
    const double u = 0.5 + 0.37 * std::sin(12.9 * i + 0.2);
    const double v = 0.5 + 0.37 * std::cos(7.1 * i + 1.3);
    const double w = 0.5 + 0.37 * std::sin(3.3 * i + 2.2);
    mix.row(row++) = Eigen::RowVector3d(40.0 + 20.0 * u, 5.0 + 20.0 * v, 40.0 + 20.0 * w);
  }
  Eigen::Matrix3d big = Eigen::Matrix3d::Identity() * 80.0;
  std::vector<int> kind;
  symm::classifyIce(mix, big, kind);
  int coreIc = 0;
  int tail = 0;
  for (int i = 0; i < nCrystal; i++) {
    if (kind[static_cast<size_t>(i)] == static_cast<int>(symm::Kind::iceIc)) {
      coreIc++;
    }
  }
  for (int i = nCrystal; i < mix.rows(); i++) {
    if (kind[static_cast<size_t>(i)] < 0) {
      tail++;
    }
  }
  std::cout << "nucleus crystal_atoms=" << nCrystal << " labelled_ice_Ic=" << coreIc
            << " liquid_atoms=" << nLiquid << " labelled_other=" << tail << "\n";
  chillHist(mix, 80.0, "nucleus");
  return 0;
}
