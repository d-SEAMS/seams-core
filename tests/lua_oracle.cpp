/* Copyright (c) d-SEAMS developers
 * SPDX-License-Identifier: MIT
 *
 * C++ census for the Lua fixture check. The calls are the ones
 * seams_chill_plus, seams_q4q6 and seams_rings_find make on a cutoff list:
 * orthoFullCsr, chillPlusFromCsr, steinhardtQlCutoff, ringNetwork.
 */

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

void iceIh(int nx, int ny, int nz, std::vector<double> &xyz, double box9[9]) {
  const double d = 2.75;
  const double z = 1.0 / 16.0;
  const double c = d / 0.375;
  const double a = c / std::sqrt(8.0 / 3.0);
  const double basis[4][3] = {{1.0 / 3.0, 2.0 / 3.0, z},
                              {2.0 / 3.0, 1.0 / 3.0, 0.5 + z},
                              {2.0 / 3.0, 1.0 / 3.0, 1.0 - z},
                              {1.0 / 3.0, 2.0 / 3.0, 0.5 - z}};
  const double lx = nx * a;
  const double ly = ny * a * std::sqrt(3.0) / 2.0;
  const double lz = nz * c;
  const double xy = ny * (-0.5 * a);
  box9[0] = box9[1] = box9[2] = 0;
  box9[3] = lx;
  box9[4] = ly;
  box9[5] = lz;
  box9[6] = xy;
  box9[7] = box9[8] = 0;
  xyz.clear();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        for (int b = 0; b < 4; b++) {
          const double sx = (i + basis[b][0]) / nx;
          const double sy = (j + basis[b][1]) / ny;
          const double sz = (k + basis[b][2]) / nz;
          xyz.push_back(sx * lx + sy * xy);
          xyz.push_back(sy * ly);
          xyz.push_back(sz * lz);
        }
      }
    }
  }
}

void writeFrame(const std::string &path, const std::vector<double> &xyz, const double box9[9]) {
  std::ofstream out(path);
  out.setf(std::ios::scientific);
  out.precision(17);
  const int n = static_cast<int>(xyz.size() / 3);
  out << "# dseams-frame\n" << n << "\n";
  for (int i = 0; i < 9; i++) {
    if (i) {
      out << " ";
    }
    out << box9[i];
  }
  out << "\n";
  for (int i = 0; i < n; i++) {
    out << xyz[3 * i] << " " << xyz[3 * i + 1] << " " << xyz[3 * i + 2] << "\n";
  }
}

Cloud cloudFrom(const std::vector<double> &xyz, const double box9[9]) {
  Cloud cloud;
  if (!nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, box9[0], box9[1], box9[2],
                                    box9[3], box9[4], box9[5], box9[6], box9[7], box9[8])) {
    throw std::runtime_error("box edge lengths must be positive");
  }
  const int n = static_cast<int>(xyz.size() / 3);
  cloud.nop = n;
  cloud.currentFrame = 1;
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.x = xyz[3 * i];
    point.y = xyz[3 * i + 1];
    point.z = xyz[3 * i + 2];
    point.type = 1;
    point.atomID = i + 1;
    point.molID = i + 1;
    cloud.idIndexMap[point.atomID] = i;
  }
  return cloud;
}

struct Report {
  int n = 0;
  int chill[9] = {};
  int rings[13] = {};
  double q6 = 0;
  int nnz = 0;
};

Report analyze(Cloud cloud, int depth) {
  Report report;
  report.n = cloud.nop;
  std::vector<int> offsets;
  std::vector<int> cols;
  if (!nneigh::orthoFullCsr(cloud, nullptr, 3.5, offsets, cols, false)) {
    throw std::runtime_error("cutoff list failed");
  }
  report.nnz = static_cast<int>(cols.size());
  chill::chillPlusFromCsr(cloud, offsets, cols);
  for (int i = 0; i < cloud.nop; i++) {
    int label = static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
    if (label < 0 || label > 8) {
      label = 6;
    }
    report.chill[label] += 1;
  }
  const auto found = primitive::ringNetwork(offsets, cols, depth);
  for (const auto &ring : found) {
    const int size = static_cast<int>(ring.size());
    if (size >= 0 && size <= depth && size < 13) {
      report.rings[size] += 1;
    }
  }
  const int typ = cloud.nop > 0 ? cloud.pts[0].type : 1;
  const auto q = chill::steinhardtQlCutoff(cloud, 3.5, 6, typ, false);
  if (static_cast<int>(q.ql.size()) < cloud.nop) {
    throw std::runtime_error("q6 length is short");
  }
  double total = 0;
  for (int i = 0; i < cloud.nop; i++) {
    total += q.ql[static_cast<std::size_t>(i)];
  }
  report.q6 = cloud.nop > 0 ? total / cloud.nop : 0;
  return report;
}

void printReport(const char *name, const Report &report, int depth) {
  std::cout << name << " n " << report.n << "\n" << name << " chill";
  for (int i = 0; i < 9; i++) {
    std::cout << " " << report.chill[i];
  }
  std::cout << "\n" << name << " rings";
  for (int i = 0; i <= depth; i++) {
    std::cout << " " << report.rings[i];
  }
  std::cout.setf(std::ios::scientific);
  std::cout.precision(12);
  std::cout << "\n" << name << " q6 " << report.q6 << "\n";
  std::cout.unsetf(std::ios::scientific);
  std::cout << name << " nnz " << report.nnz << "\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: lua_oracle frame-path source-root\n";
    return 2;
  }
#ifdef _OPENMP
  omp_set_num_threads(1);
#endif
  const std::string frame = argv[1];
  const std::string root = argv[2];
  try {
    std::vector<double> xyz;
    double box9[9];
    iceIh(4, 4, 2, xyz, box9);
    writeFrame(frame, xyz, box9);
    printReport("ih", analyze(cloudFrom(xyz, box9), 6), 6);
    Cloud guest;
    guest = sinp::readLammpsTrjO(root + "/input/traj/genice_sI.lammpstrj", 1, guest, 1);
    printReport("si", analyze(guest, 6), 6);
    Cloud face;
    face = sinp::readLammpsTrjO(root + "/input/traj/interface_f0.lammpstrj", 1, face, 1);
    printReport("interface", analyze(face, 6), 6);
    Cloud tip;
    tip = sinp::readLammpsTrjO(root + "/input/traj/exampleTraj.lammpstrj", 1, tip, 2);
    printReport("tip", analyze(tip, 8), 8);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 1;
  }
  return 0;
}
