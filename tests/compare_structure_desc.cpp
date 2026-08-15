/*
** Print d-SEAMS Steinhardt, Voronoi, and template counts on FCC and BCC
** lattices. Pair with scripts/compare_ql_literature.py on the same geometry.
** Run: bcpp/tests/compare_structure_desc
*/

#include <bop.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <structure_desc.hpp>
#include <voronoi_qlm.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud lattice(const std::vector<std::array<double, 3>> &basis, int reps,
              double a) {
  Cloud cloud;
  int id = 1;
  for (int i = 0; i < reps; i++) {
    for (int j = 0; j < reps; j++) {
      for (int k = 0; k < reps; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = (i + b[0]) * a;
          p.y = (j + b[1]) * a;
          p.z = (k + b[2]) * a;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  const double L = reps * a;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

Cloud fcc() {
  return lattice(
      {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}}, 3,
      4.0);
}

Cloud bcc() {
  return lattice({{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, 3, 4.0);
}

double mean(const std::vector<double> &v) {
  double s = 0.0;
  for (double x : v) {
    s += x;
  }
  return v.empty() ? 0.0 : s / static_cast<double>(v.size());
}

void report(const char *name, const Cloud &cloud, double cutQl, double cutVor,
            int kNeigh) {
  auto nList = nneigh::neighListO(cutQl, cloud, 1);
  const auto q4 = chill::steinhardtQl(cloud, nList, 4);
  const auto q6 = chill::steinhardtQl(cloud, nList, 6);
  const auto q8 = chill::steinhardtQl(cloud, nList, 8);
  const auto v4 = chill::steinhardtQlVoronoi(cloud, cutVor, 4);
  const auto v6 = chill::steinhardtQlVoronoi(cloud, cutVor, 6);
  const auto v8 = chill::steinhardtQlVoronoi(cloud, cutVor, 8);
  const auto hits = chill::classifyTemplates(cloud, nList, kNeigh);
  int fcc = 0, hcp = 0, bccN = 0, sc = 0, other = 0;
  for (const auto &h : hits) {
    switch (h.kind) {
    case chill::CrystalKind::fcc:
      fcc++;
      break;
    case chill::CrystalKind::hcp:
      hcp++;
      break;
    case chill::CrystalKind::bcc:
      bccN++;
      break;
    case chill::CrystalKind::sc:
      sc++;
      break;
    default:
      other++;
      break;
    }
  }
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "tool=dseams lattice=" << name << " n=" << cloud.nop
            << " L=" << cloud.box[0] << " cut_ql=" << cutQl
            << " cut_vor=" << cutVor << "\n";
  std::cout << "  ql4=" << mean(q4.ql) << " ql6=" << mean(q6.ql)
            << " ql8=" << mean(q8.ql) << "\n";
  std::cout << "  vor_q4=" << mean(v4.ql) << " vor_q6=" << mean(v6.ql)
            << " vor_q8=" << mean(v8.ql) << "\n";
  std::cout << "  templates fcc=" << fcc << " hcp=" << hcp << " bcc=" << bccN
            << " sc=" << sc << " other=" << other << "\n";
}

} // namespace

int main() {
  report("fcc", fcc(), 3.2, 4.8, 12);
  report("bcc", bcc(), 4.0, 5.6, 8);
  return 0;
}
