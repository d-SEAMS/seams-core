#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ylm_cartesian.hpp>

#include <array>
#include <cmath>
#include <random>
#include <vector>

namespace {

molSys::PointCloud<molSys::Point<double>, double> iceIc(int cells) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  int id = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (int b = 0; b < 8; b++) {
          molSys::Point<double> p;
          const int f = b % 4;
          const double shift = b < 4 ? 0.0 : 0.25;
          p.x = std::fmod((i + fcc[f][0] + shift) * lattice, side);
          p.y = std::fmod((j + fcc[f][1] + shift) * lattice, side);
          p.z = std::fmod((k + fcc[f][2] + shift) * lattice, side);
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          cloud.idIndexMap[id] = cloud.nop;
          cloud.pts.push_back(p);
          ++cloud.nop;
          ++id;
        }
      }
    }
  }
  return cloud;
}

molSys::PointCloud<molSys::Point<double>, double> iceIhCell() {
  const auto cell = symm::crystal(symm::Crystal::iceIh);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double a = cell.vectors(0, 0);
  const double xy = cell.vectors(1, 0);
  const double ly = cell.vectors(1, 1);
  const double lz = cell.vectors(2, 2);
  const double xmin = std::min(0.0, xy);
  const double xmax = std::max(0.0, xy);
  cloud.box = {a + xmax - xmin, ly, lz, xy, 0.0, 0.0};
  cloud.boxLow = {xmin, 0.0, 0.0};
  cloud.currentFrame = 1;
  int id = 1;
  for (Eigen::Index i = 0; i < cell.cartesian.rows(); ++i) {
    molSys::Point<double> p;
    p.x = cell.cartesian(i, 0);
    p.y = cell.cartesian(i, 1);
    p.z = cell.cartesian(i, 2);
    p.type = 1;
    p.atomID = id;
    p.molID = id;
    cloud.idIndexMap[id] = cloud.nop;
    cloud.pts.push_back(p);
    ++cloud.nop;
    ++id;
  }
  return cloud;
}

molSys::PointCloud<molSys::Point<double>, double> liquidBox(int n, unsigned seed) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(0.0, 20.0);
  for (int i = 0; i < n; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    cloud.idIndexMap[p.atomID] = cloud.nop;
    cloud.pts.push_back(p);
    ++cloud.nop;
  }
  return cloud;
}

void expectSameLabels(const molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  auto ref = cloud;
  const auto rows = nneigh::neighListO(3.5, ref, 1);
  chill::getCorrelPlus(ref, rows, false, 4);
  chill::getIceTypePlusNoPrint(ref, rows, false, true);
  auto got = cloud;
  const auto again = nneigh::neighListO(3.5, got, 1);
  chill::applyChillPlus(got, again, chill::ChillPlusRule::Published, false);
  REQUIRE(got.nop == ref.nop);
  for (int i = 0; i < ref.nop; i++) {
    REQUIRE(got.pts[static_cast<std::size_t>(i)].molID ==
            ref.pts[static_cast<std::size_t>(i)].molID);
    REQUIRE(got.pts[static_cast<std::size_t>(i)].iceType ==
            ref.pts[static_cast<std::size_t>(i)].iceType);
  }
}

double maxAbs(const std::vector<double> &a, const std::vector<double> &b) {
  double worst = 0.0;
  for (std::size_t i = 0; i < a.size(); i++) {
    worst = std::max(worst, std::fabs(a[i] - b[i]));
  }
  return worst;
}

} // namespace

TEST_CASE("solid q4 and q6 match the neighbour-list path", "[ql][solid]") {
  auto cloud = iceIc(2);
  const auto rows = nneigh::kNearestNeighbourList(cloud, 4, 5.0, 1, true);
  for (int order : {3, 4, 6}) {
    const auto ref = chill::steinhardtQl(cloud, rows, order);
    const auto got = chill::steinhardtQlMutual(cloud, 4, 5.0, 1, order);
    REQUIRE(maxAbs(ref.ql, got.ql) < 1e-12);
    REQUIRE(maxAbs(ref.qlBar, got.qlBar) < 1e-12);
  }
}

TEST_CASE("solid q4 and q6 match on a disordered cloud", "[ql][solid]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {18.0, 18.0, 18.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  std::mt19937 rng(91001);
  std::uniform_real_distribution<double> u(0.0, 18.0);
  for (int i = 0; i < 80; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    cloud.idIndexMap[p.atomID] = cloud.nop;
    cloud.pts.push_back(p);
    ++cloud.nop;
  }
  const auto rows = nneigh::kNearestNeighbourList(cloud, 4, 6.0, 1, true);
  for (int order : {4, 6}) {
    const auto ref = chill::steinhardtQl(cloud, rows, order);
    const auto got = chill::steinhardtQlMutual(cloud, 4, 6.0, 1, order);
    const auto again = chill::steinhardtQlMutual(cloud, 4, 6.0, 1, order);
    REQUIRE(maxAbs(ref.ql, got.ql) < 1e-12);
    REQUIRE(maxAbs(ref.qlBar, got.qlBar) < 1e-12);
    REQUIRE(maxAbs(got.ql, again.ql) < 1e-12);
    REQUIRE(maxAbs(got.qlBar, again.qlBar) < 1e-12);
  }
}

TEST_CASE("published CHILL+ CSR labels match the row rule per molecule", "[chill][solid]") {
  expectSameLabels(iceIc(2));
  expectSameLabels(iceIhCell());
  expectSameLabels(liquidBox(80, 91001));
  auto mixed = iceIc(2);
  const double side = mixed.box[0];
  std::mt19937 rng(93003);
  std::uniform_real_distribution<double> u(0.0, side);
  const int base = mixed.nop;
  for (int i = 0; i < 24; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = base + i + 1;
    p.molID = p.atomID;
    mixed.idIndexMap[p.atomID] = mixed.nop;
    mixed.pts.push_back(p);
    ++mixed.nop;
  }
  expectSameLabels(mixed);
}

TEST_CASE("symbolic cubic pre-filter falls back outside its certificate", "[chill][prefilter]") {
  auto ice = iceIc(2);
  auto viaFilter = ice;
  REQUIRE(chill::chillPlusSymbolicPrefilter(viaFilter, 3.5, 1));
  REQUIRE(chill::chillPlusPublishedCutoff(ice, 3.5, 1));
  for (int i = 0; i < ice.nop; i++) {
    REQUIRE(viaFilter.pts[static_cast<std::size_t>(i)].iceType ==
            molSys::atom_state_type::cubic);
    REQUIRE(ice.pts[static_cast<std::size_t>(i)].iceType ==
            molSys::atom_state_type::cubic);
    const auto &a = viaFilter.pts[static_cast<std::size_t>(i)].c_ij;
    const auto &b = ice.pts[static_cast<std::size_t>(i)].c_ij;
    REQUIRE(a.size() == b.size());
    for (std::size_t k = 0; k < a.size(); k++) {
      REQUIRE(std::fabs(a[k].c_value - b[k].c_value) < 1e-12);
    }
  }
  molSys::PointCloud<molSys::Point<double>, double> liquid;
  liquid.box = {20.0, 20.0, 20.0};
  liquid.boxLow = {0.0, 0.0, 0.0};
  liquid.currentFrame = 1;
  std::mt19937 rng(91001);
  std::uniform_real_distribution<double> u(0.0, 20.0);
  for (int i = 0; i < 70; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    liquid.idIndexMap[p.atomID] = liquid.nop;
    liquid.pts.push_back(p);
    ++liquid.nop;
  }
  auto untouched = liquid;
  REQUIRE_FALSE(chill::chillPlusSymbolicPrefilter(liquid, 3.5, 1));
  for (int i = 0; i < liquid.nop; i++) {
    REQUIRE(liquid.pts[static_cast<std::size_t>(i)].iceType ==
            untouched.pts[static_cast<std::size_t>(i)].iceType);
    REQUIRE(liquid.pts[static_cast<std::size_t>(i)].c_ij.size() ==
            untouched.pts[static_cast<std::size_t>(i)].c_ij.size());
  }
}

TEST_CASE("published cutoff CHILL+ matches the neighbour-list rule", "[chill][solid]") {
  auto check = [](molSys::PointCloud<molSys::Point<double>, double> cloud) {
    auto ref = cloud;
    const auto rows = nneigh::neighListO(3.5, ref, 1);
    chill::applyChillPlus(ref, rows, chill::ChillPlusRule::Published, false);
    REQUIRE(chill::chillPlusPublishedCutoff(cloud, 3.5, 1));
    for (int i = 0; i < cloud.nop; i++) {
      REQUIRE(cloud.pts[static_cast<std::size_t>(i)].iceType ==
              ref.pts[static_cast<std::size_t>(i)].iceType);
      const auto &a = cloud.pts[static_cast<std::size_t>(i)].c_ij;
      const auto &b = ref.pts[static_cast<std::size_t>(i)].c_ij;
      REQUIRE(a.size() == b.size());
      for (std::size_t k = 0; k < a.size(); k++) {
        REQUIRE(std::fabs(a[k].c_value - b[k].c_value) < 1e-12);
      }
    }
  };
  check(iceIc(2));
  molSys::PointCloud<molSys::Point<double>, double> liquid;
  liquid.box = {20.0, 20.0, 20.0};
  liquid.boxLow = {0.0, 0.0, 0.0};
  liquid.currentFrame = 1;
  std::mt19937 rng(91001);
  std::uniform_real_distribution<double> u(0.0, 20.0);
  for (int i = 0; i < 70; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    liquid.idIndexMap[p.atomID] = liquid.nop;
    liquid.pts.push_back(p);
    ++liquid.nop;
  }
  check(liquid);
}

TEST_CASE("Cartesian l=4 matches the trigonometric harmonic", "[ql][solid]") {
  const std::vector<std::array<double, 3>> dirs = {
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 2.5}, {0.0, 0.0, -1.2},
      {0.2, -0.4, 0.7}, {3.0, 0.1, 0.05}};
  for (const auto &d : dirs) {
    const auto angles = sph::radialCoord({d[0], d[1], d[2]});
    const auto host = sph::spheriHarmo(4, angles);
    double re[9];
    double im[9];
    REQUIRE(seams::ylm::ylmReal(4, d[0], d[1], d[2], re, im));
    for (int m = 0; m < 9; m++) {
      REQUIRE_THAT(re[m], Catch::Matchers::WithinAbs(host[static_cast<std::size_t>(m)].real(), 1e-12));
      REQUIRE_THAT(im[m], Catch::Matchers::WithinAbs(host[static_cast<std::size_t>(m)].imag(), 1e-12));
    }
  }
}

TEST_CASE("fused q4 q6 and q-bar match the cutoff Steinhardt", "[ql][chill]") {
  auto check = [](molSys::PointCloud<molSys::Point<double>, double> cloud) {
    auto labeled = cloud;
    const auto q4 = chill::steinhardtQlCutoff(cloud, 3.5, 4, 1, true);
    const auto q6 = chill::steinhardtQlCutoff(cloud, 3.5, 6, 1, true);
    REQUIRE(chill::chillPlusPublishedCutoff(labeled, 3.5, 1));
    std::vector<double> fq4, fq6, fq4bar, fq6bar;
    auto fused = cloud;
    REQUIRE(chill::q4q6ChillPlus(fused, 3.5, 1, fq4, fq6, fq4bar, fq6bar));
    REQUIRE(maxAbs(q4.ql, fq4) < 1e-12);
    REQUIRE(maxAbs(q6.ql, fq6) < 1e-12);
    REQUIRE(maxAbs(q4.qlBar, fq4bar) < 1e-12);
    REQUIRE(maxAbs(q6.qlBar, fq6bar) < 1e-12);
    for (int i = 0; i < cloud.nop; i++) {
      REQUIRE(fused.pts[static_cast<std::size_t>(i)].iceType ==
              labeled.pts[static_cast<std::size_t>(i)].iceType);
    }
  };
  check(iceIc(2));
  molSys::PointCloud<molSys::Point<double>, double> liquid;
  liquid.box = {20.0, 20.0, 20.0};
  liquid.boxLow = {0.0, 0.0, 0.0};
  liquid.currentFrame = 1;
  std::mt19937 rng(91001);
  std::uniform_real_distribution<double> u(0.0, 20.0);
  for (int i = 0; i < 70; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    liquid.idIndexMap[p.atomID] = liquid.nop;
    liquid.pts.push_back(p);
    ++liquid.nop;
  }
  check(liquid);
}

TEST_CASE("fused q6 on a triclinic cell matches the cutoff call", "[ql][chill]") {
  auto cloud = iceIc(4);
  const double side = cloud.box[0];
  const double xy = 1.2;
  const double xz = 0.4;
  const double yz = 0.3;
  const double xmin = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
  const double xmax = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
  const double ymin = std::min(0.0, yz);
  const double ymax = std::max(0.0, yz);
  for (auto &p : cloud.pts) {
    const double px = p.x;
    const double py = p.y;
    const double pz = p.z;
    p.x = px + (xy / side) * py + (xz / side) * pz;
    p.y = py + (yz / side) * pz;
  }
  cloud.box = {side + xmax - xmin, side + ymax - ymin, side, xy, xz, yz};
  auto labeled = cloud;
  const auto q6 = chill::steinhardtQlCutoff(cloud, 3.5, 6, 1, true);
  REQUIRE(chill::chillPlusPublishedCutoff(labeled, 3.5, 1));
  std::vector<double> fq4, fq6, fq4bar, fq6bar;
  auto fused = cloud;
  REQUIRE(chill::q4q6ChillPlus(fused, 3.5, 1, fq4, fq6, fq4bar, fq6bar));
  REQUIRE(maxAbs(q6.ql, fq6) < 1e-12);
  REQUIRE(maxAbs(q6.qlBar, fq6bar) < 1e-12);
  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE(fused.pts[static_cast<std::size_t>(i)].iceType ==
            molSys::atom_state_type::cubic);
    REQUIRE(labeled.pts[static_cast<std::size_t>(i)].iceType ==
            molSys::atom_state_type::cubic);
  }
}

TEST_CASE("local q4 and q6 do not need the averaged row", "[ql][solid]") {
  auto ice = iceIc(3);
  molSys::PointCloud<molSys::Point<double>, double> liquid;
  liquid.box = {18.0, 18.0, 18.0};
  liquid.boxLow = {0.0, 0.0, 0.0};
  liquid.currentFrame = 1;
  std::mt19937 rng(91001);
  std::uniform_real_distribution<double> u(0.0, 18.0);
  for (int i = 0; i < 400; i++) {
    molSys::Point<double> p;
    p.x = u(rng);
    p.y = u(rng);
    p.z = u(rng);
    p.type = 1;
    p.atomID = i + 1;
    p.molID = p.atomID;
    liquid.pts.push_back(p);
  }
  liquid.nop = static_cast<int>(liquid.pts.size());
  for (const auto *cloud : {&ice, &liquid}) {
    for (int order : {4, 6}) {
      const auto full = chill::steinhardtQlCutoff(*cloud, 3.5, order, 1, true);
      const auto local = chill::steinhardtQlCutoff(*cloud, 3.5, order, 1, false);
      REQUIRE(maxAbs(full.ql, local.ql) < 1e-12);
    }
  }
}
