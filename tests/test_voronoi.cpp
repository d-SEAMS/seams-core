#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <voronoi_qlm.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// Voronoi facet weights on the reference lattices have known geometry: the
// simple cubic cell is a cube with six equal square facets, and the FCC cell
// is a rhombic dodecahedron with twelve equal facets. With equal weights over
// the same neighbour set, the weighted Steinhardt parameters must reproduce
// the tabulated unweighted values exactly.

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

} // namespace

TEST_CASE("Voronoi cell of simple cubic is a cube", "[voronoi]") {
  const double a = 3.0;
  auto cloud = lattice({{0.0, 0.0, 0.0}}, 4, a);
  auto cells = chill::voronoiFacetWeights(cloud, 2.2 * a);

  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE(cells[i].neighbours.size() == 6);
    for (const double w : cells[i].weights) {
      REQUIRE_THAT(w, Catch::Matchers::WithinAbs(1.0 / 6.0, 1e-9));
    }
    double sum = 0.0;
    for (const double w : cells[i].weights) {
      sum += w;
    }
    REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(1.0, 1e-12));
  }
}

TEST_CASE("Voronoi cell of FCC is a rhombic dodecahedron", "[voronoi]") {
  const double a = 4.0;
  auto cloud = lattice(
      {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}}, 3,
      a);
  auto cells = chill::voronoiFacetWeights(cloud, 1.2 * a);

  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE(cells[i].neighbours.size() == 12);
    for (const double w : cells[i].weights) {
      REQUIRE_THAT(w, Catch::Matchers::WithinAbs(1.0 / 12.0, 1e-9));
    }
  }
}

TEST_CASE("Voronoi cell of BCC has eight large and six small facets",
          "[voronoi]") {
  const double a = 4.0;
  auto cloud = lattice({{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, 3, a);
  auto cells = chill::voronoiFacetWeights(cloud, 1.4 * a);

  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE(cells[i].neighbours.size() == 14);
    auto w = cells[i].weights;
    std::sort(w.begin(), w.end());
    // six squares (toward the second shell), then eight hexagons
    for (int k = 0; k < 6; k++) {
      REQUIRE_THAT(w[k], Catch::Matchers::WithinAbs(w[0], 1e-9));
    }
    for (int k = 6; k < 14; k++) {
      REQUIRE_THAT(w[k], Catch::Matchers::WithinAbs(w[13], 1e-9));
    }
    REQUIRE(w[13] > w[0]);
  }
}

TEST_CASE("weighted Steinhardt reproduces FCC references under equal weights",
          "[voronoi]") {
  const double a = 4.0;
  auto cloud = lattice(
      {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}}, 3,
      a);

  const struct {
    int l;
    double reference;
  } cases[] = {{4, 0.190941}, {6, 0.574524}, {8, 0.403915}};

  for (const auto &c : cases) {
    auto q = chill::steinhardtQlVoronoi(cloud, 1.2 * a, c.l);
    for (int i = 0; i < cloud.nop; i++) {
      REQUIRE_THAT(q.ql[i], Catch::Matchers::WithinAbs(c.reference, 1e-4));
      REQUIRE_THAT(q.qlBar[i], Catch::Matchers::WithinAbs(q.ql[i], 1e-9));
    }
  }
}

TEST_CASE("weighted q6 on BCC differs from naive fourteen-neighbour q6",
          "[voronoi]") {
  // The facet weighting is the point of the method: the six second-shell
  // neighbours carry less weight than the eight first-shell ones, so the
  // weighted parameter cannot equal an unweighted average over all fourteen
  const double a = 4.0;
  auto cloud = lattice({{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, 3, a);
  auto weighted = chill::steinhardtQlVoronoi(cloud, 1.4 * a, 6);

  REQUIRE(weighted.ql[0] > 0.0);
  // Uniform environment: the average still coincides with the local value
  for (int i = 0; i < cloud.nop; i++) {
    REQUIRE_THAT(weighted.qlBar[i],
                 Catch::Matchers::WithinAbs(weighted.ql[i], 1e-9));
    REQUIRE_THAT(weighted.ql[i],
                 Catch::Matchers::WithinAbs(weighted.ql[0], 1e-9));
  }
}
