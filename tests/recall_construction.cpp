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

#include <algorithm>
#include <cmath>
#include <fstream>
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
  int n = 0;
  double precision = 0;
  double recall = 0;
};

// Cubic diamond, fixed before either classifier runs. The first shell is
// 4 neighbours at 2.75 A. The second shell is 12 neighbours at a/sqrt(2).
// A site is interior only when the minimum-image cell holds all 16. A
// short shell is an open face of this block, and neither classifier is
// scored on it.
constexpr double kFirst = 2.75;
constexpr double kShellTol = 0.20;

double secondShell(double lattice) { return lattice / std::sqrt(2.0); }

double minImage(double delta, double length) {
  if (!(length > 0.0)) {
    return delta;
  }
  return delta - length * std::round(delta / length);
}

int latticeNeighbours(const Cloud &cloud, int i, const double box[3], double second) {
  const auto &pi = cloud.pts[static_cast<std::size_t>(i)];
  const double limit = second + kShellTol;
  int n = 0;
  for (int j = 0; j < cloud.nop; ++j) {
    if (j == i) {
      continue;
    }
    const auto &pj = cloud.pts[static_cast<std::size_t>(j)];
    const double dx = minImage(pi.x - pj.x, box[0]);
    const double dy = minImage(pi.y - pj.y, box[1]);
    const double dz = minImage(pi.z - pj.z, box[2]);
    if (dx * dx + dy * dy + dz * dz <= limit * limit) {
      ++n;
    }
  }
  return n;
}

// Distance to the nearest face that the periodic cell does not seal.
// A face is open when the cell is longer than the block plus one second
// shell, so the crystal does not tile that axis.
double surfaceDistance(const Cloud &block, int i, const double box[3], double second) {
  double lo[3] = {block.pts[0].x, block.pts[0].y, block.pts[0].z};
  double hi[3] = {lo[0], lo[1], lo[2]};
  for (const auto &point : block.pts) {
    const double c[3] = {point.x, point.y, point.z};
    for (int a = 0; a < 3; ++a) {
      lo[a] = std::min(lo[a], c[a]);
      hi[a] = std::max(hi[a], c[a]);
    }
  }
  const double c[3] = {block.pts[static_cast<std::size_t>(i)].x,
                       block.pts[static_cast<std::size_t>(i)].y,
                       block.pts[static_cast<std::size_t>(i)].z};
  double best = 1.0e9;
  bool any = false;
  for (int a = 0; a < 3; ++a) {
    if (box[a] <= (hi[a] - lo[a]) + second) {
      continue;
    }
    any = true;
    best = std::min(best, std::min(c[a] - lo[a], hi[a] - c[a]));
  }
  return any ? best : 0.0;
}

Score score(const std::vector<char> &truth, const std::vector<char> &pred,
            const std::vector<char> *keep) {
  Score out;
  const int n = static_cast<int>(truth.size());
  for (int i = 0; i < n; ++i) {
    if (keep != nullptr && (*keep)[static_cast<std::size_t>(i)] == 0) {
      continue;
    }
    ++out.n;
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

void printRow(const char *scope, const char *method, const Score &row) {
  std::cout << scope << "," << method << "," << row.tp << "," << row.fp << "," << row.fn << ","
            << row.precision << "," << row.recall << "\n";
}

} // namespace

int main(int argc, char **argv) {
  const std::string sitesPath = argc > 1 ? argv[1] : "";
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double second = secondShell(a);
  Cloud ice = diamond(4);
  // The mask is the lattice and the cell. Classifiers have not run.
  const double iceBox[3] = {ice.box[0], ice.box[1], ice.box[2]};
  int periodicShort = 0;
  for (int i = 0; i < ice.nop; ++i) {
    if (latticeNeighbours(ice, i, iceBox, second) != 16) {
      ++periodicShort;
    }
  }

  Cloud liquid = spacedLiquid(ice.nop);
  const double shift = ice.box[0] + 20.0;
  for (auto &point : liquid.pts) {
    point.x += shift;
  }
  std::vector<char> truth(static_cast<std::size_t>(ice.nop), 1);
  truth.resize(static_cast<std::size_t>(ice.nop + liquid.nop), 0);
  Cloud cloud = ice;
  cloud.box[0] = shift + liquid.box[0];
  append(cloud, liquid);
  const double box[3] = {cloud.box[0], cloud.box[1], cloud.box[2]};

  std::vector<char> interior(static_cast<std::size_t>(cloud.nop), 0);
  std::vector<int> nLattice(static_cast<std::size_t>(cloud.nop), 0);
  std::vector<double> surface(static_cast<std::size_t>(cloud.nop), 0);
  int nInterior = 0;
  for (int i = 0; i < ice.nop; ++i) {
    nLattice[static_cast<std::size_t>(i)] = latticeNeighbours(ice, i, box, second);
    surface[static_cast<std::size_t>(i)] = surfaceDistance(ice, i, box, second);
    if (nLattice[static_cast<std::size_t>(i)] == 16) {
      interior[static_cast<std::size_t>(i)] = 1;
      ++nInterior;
    }
  }
  for (int i = ice.nop; i < cloud.nop; ++i) {
    interior[static_cast<std::size_t>(i)] = 1;
  }

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

  if (!sitesPath.empty()) {
    std::ofstream out(sitesPath);
    if (!out) {
      std::cerr << "cannot write " << sitesPath << "\n";
      return 2;
    }
    out << "site,x,y,z,n_lattice,surface_A,interior,dseams,dseams_kind,chill\n";
    out.setf(std::ios::fixed);
    out.precision(6);
    for (int i = 0; i < ice.nop; ++i) {
      const auto &point = ice.pts[static_cast<std::size_t>(i)];
      out << i << "," << point.x << "," << point.y << "," << point.z << ","
          << nLattice[static_cast<std::size_t>(i)] << "," << surface[static_cast<std::size_t>(i)]
          << "," << int(interior[static_cast<std::size_t>(i)]) << ","
          << int(kernel[static_cast<std::size_t>(i)]) << "," << shell[static_cast<std::size_t>(i)]
          << "," << int(chill[static_cast<std::size_t>(i)]) << "\n";
    }
  }

  std::cout << "scope,method,tp,fp,fn,precision,recall\n";
  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  printRow("all", "d-seams", score(truth, kernel, nullptr));
  printRow("all", "chill+", score(truth, chill, nullptr));
  printRow("interior", "d-seams", score(truth, kernel, &interior));
  printRow("interior", "chill+", score(truth, chill, &interior));
  std::cout << "periodic_complete " << (periodicShort == 0 ? "yes" : "no") << " short "
            << periodicShort << " interior " << nInterior << " of " << ice.nop << "\n";
  const Score inside = score(truth, kernel, &interior);
  const Score everyone = score(truth, kernel, nullptr);
  if (periodicShort != 0 || inside.fn != 0 || everyone.precision < 1.0 - 1e-9 ||
      everyone.recall < 0.875 - 1e-9) {
    return 3;
  }
  return 0;
}
