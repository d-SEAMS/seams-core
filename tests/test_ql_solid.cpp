#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
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
    REQUIRE(maxAbs(ref.ql, got.ql) < 1e-12);
    REQUIRE(maxAbs(ref.qlBar, got.qlBar) < 1e-12);
  }
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
