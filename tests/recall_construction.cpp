//-----------------------------------------------------------------------------------
// d-SEAMS
// SPDX-License-Identifier: MIT
//
// Precision and recall of the symmetry kernel and of published CHILL+
// against labels that come from the construction of the frame. The ice
// block is an ideal cubic lattice. The liquid block is a lattice whose
// spacing is outside the 3.5 A cutoff. Neither label is a vote of the
// two classifiers.
//-----------------------------------------------------------------------------------

#include <chill_plus_scheme.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud diamond(int nrep) {
  const double oo = 2.75;
  const double a = oo * 4.0 / std::sqrt(3.0);
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  Cloud cloud;
  cloud.box = {nrep * a, nrep * a, nrep * a};
  cloud.boxLow = {0, 0, 0};
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (const auto &shift : fcc) {
          for (const auto &site : basis) {
            molSys::Point<double> p;
            p.x = (i + shift[0] + site[0]) * a;
            p.y = (j + shift[1] + site[1]) * a;
            p.z = (k + shift[2] + site[2]) * a;
            p.type = 1;
            p.atomID = cloud.nop + 1;
            p.molID = p.atomID;
            cloud.idIndexMap[p.atomID] = cloud.nop;
            cloud.pts.push_back(p);
            cloud.nop += 1;
          }
        }
      }
    }
  }
  return cloud;
}

Cloud spacedLiquid(int n) {
  Cloud cloud;
  const double step = 8.0;
  cloud.box = {n * step, step, step};
  cloud.boxLow = {0, 0, 0};
  for (int i = 0; i < n; ++i) {
    molSys::Point<double> p;
    p.x = (i + 0.5) * step;
    p.y = 0.5 * step;
    p.z = 0.5 * step;
    p.type = 1;
    p.atomID = cloud.nop + 1;
    p.molID = p.atomID;
    cloud.idIndexMap[p.atomID] = cloud.nop;
    cloud.pts.push_back(p);
    cloud.nop += 1;
  }
  return cloud;
}

void append(Cloud &dst, const Cloud &src) {
  for (const auto &point : src.pts) {
    auto p = point;
    p.atomID = dst.nop + 1;
    p.molID = p.atomID;
    dst.idIndexMap[p.atomID] = dst.nop;
    dst.pts.push_back(p);
    dst.nop += 1;
  }
}

struct Score {
  int tp = 0;
  int fp = 0;
  int fn = 0;
  double precision = 0;
  double recall = 0;
};

Score score(const std::vector<char> &truth, const std::vector<char> &pred) {
  Score out;
  const int n = static_cast<int>(truth.size());
  for (int i = 0; i < n; ++i) {
    const bool ice = truth[static_cast<std::size_t>(i)] != 0;
    const bool hit = pred[static_cast<std::size_t>(i)] != 0;
    if (ice && hit) {
      ++out.tp;
    } else if (hit) {
      ++out.fp;
    } else if (ice) {
      ++out.fn;
    }
  }
  out.precision = (out.tp + out.fp) > 0 ? static_cast<double>(out.tp) / (out.tp + out.fp) : 1.0;
  out.recall = (out.tp + out.fn) > 0 ? static_cast<double>(out.tp) / (out.tp + out.fn) : 1.0;
  return out;
}

void printRow(const char *method, const Score &row) {
  std::cout << method << "," << row.tp << "," << row.fp << "," << row.fn << ","
            << row.precision << "," << row.recall << "\n";
}

} // namespace

int main() {
  Cloud ice = diamond(4);
  Cloud liquid = spacedLiquid(ice.nop);
  // Shift the liquid off the ice block so the two constructions do not overlap.
  const double shift = ice.box[0] + 20.0;
  for (auto &point : liquid.pts) {
    point.x += shift;
  }
  std::vector<char> truth(static_cast<std::size_t>(ice.nop), 1);
  truth.resize(static_cast<std::size_t>(ice.nop + liquid.nop), 0);
  Cloud cloud = ice;
  cloud.box[0] = shift + liquid.box[0];
  append(cloud, liquid);

  auto rows = nneigh::neighListO(3.5, cloud, 1);
  chill::applyChillPlus(cloud, rows, chill::ChillPlusRule::Published, false);
  std::vector<char> chill(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; ++i) {
    const auto kind = cloud.pts[static_cast<std::size_t>(i)].iceType;
    chill[static_cast<std::size_t>(i)] =
        (kind == molSys::atom_state_type::cubic || kind == molSys::atom_state_type::hexagonal) ? 1
                                                                                                : 0;
  }
  std::vector<int> shell;
  symm::classifyIceCloud(cloud, shell);
  std::vector<char> kernel(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; ++i) {
    kernel[static_cast<std::size_t>(i)] = shell[static_cast<std::size_t>(i)] >= 0 ? 1 : 0;
  }

  std::cout << "method,tp,fp,fn,precision,recall\n";
  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  printRow("d-seams", score(truth, kernel));
  printRow("chill+", score(truth, chill));
  return 0;
}
