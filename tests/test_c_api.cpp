#include <catch2/catch_test_macros.hpp>
#include <seams_c_api.h>

#include <cmath>
#include <vector>

TEST_CASE("seams_chill_plus labels a four-atom tetrahedron without crashing",
          "[c_api]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  int labels[4] = {-1, -1, -1, -1};
  REQUIRE(seams_chill_plus(xyz, 4, box, labels) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(labels[i] >= 0);
    REQUIRE(labels[i] <= 8);
  }
}

TEST_CASE("seams_chill_plus is seams_chill_plus_graph on the 3.5 A cutoff list",
          "[c_api]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  int byDefault[4] = {-1, -1, -1, -1};
  int byCutoff[4] = {-2, -2, -2, -2};
  int byNull[4] = {-3, -3, -3, -3};
  int byMutual[4] = {-4, -4, -4, -4};
  REQUIRE(seams_chill_plus(xyz, 4, box, byDefault) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, "cutoff", byCutoff) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, nullptr, byNull) == 0);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, "mutual", byMutual) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(byDefault[i] == byCutoff[i]);
    REQUIRE(byNull[i] == byCutoff[i]);
    REQUIRE(byMutual[i] >= 0);
    REQUIRE(byMutual[i] <= 8);
  }
  int unused[4] = {};
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 3.5, "knn-union", unused) == 1);
  REQUIRE(seams_chill_plus_graph(xyz, 4, box, 0.0, "cutoff", unused) == 1);
  const double zero[6] = {0, 0, 0, 0, 0, 0};
  REQUIRE(seams_chill_plus_graph(xyz, 4, zero, 3.5, "cutoff", unused) == 1);
}

TEST_CASE("seams_chill_plus_graph_lammps matches the ortho call and accepts tilt",
          "[c_api]") {
  const double xyz[6] = {0, 0, 0, 2.0, 0, 0};
  const double ortho[6] = {0, 0, 0, 20, 20, 20};
  const double box9[9] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
  int a[2] = {-1, -1};
  int b[2] = {-1, -1};
  REQUIRE(seams_chill_plus_graph(xyz, 2, ortho, 3.5, "cutoff", a) == 0);
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 2, box9, 3.5, "cutoff", b) == 0);
  REQUIRE(a[0] == b[0]);
  REQUIRE(a[1] == b[1]);
  const double tilted[9] = {0, 0, 0, 10, 8.660254037844386, 10, 5, 0, 0};
  REQUIRE(seams_chill_plus_graph_lammps(xyz, 2, tilted, 3.5, "cutoff", b) == 0);
}

TEST_CASE("seams_q4_q6 returns finite cutoff values and rejects a bad box",
          "[c_api]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  const double box9[9] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
  double q4[4] = {};
  double q6[4] = {};
  double q4bar[4] = {};
  double q6bar[4] = {};
  double q4b[4] = {};
  double q6b[4] = {};
  REQUIRE(seams_q4_q6(xyz, 4, box, 3.5, q4, q6, q4bar, q6bar) == 0);
  REQUIRE(seams_q4_q6_lammps(xyz, 4, box9, 3.5, q4b, q6b, nullptr, nullptr) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(std::isfinite(q4[i]));
    REQUIRE(std::isfinite(q6[i]));
    REQUIRE(std::isfinite(q4bar[i]));
    REQUIRE(std::isfinite(q6bar[i]));
    REQUIRE(q4[i] == q4b[i]);
    REQUIRE(q6[i] == q6b[i]);
  }
  REQUIRE(seams_q4_q6(nullptr, 4, box, 3.5, q4, q6, nullptr, nullptr) == 1);
  REQUIRE(seams_q4_q6(xyz, 4, box, 0.0, q4, q6, nullptr, nullptr) == 1);
  const double zero[6] = {0, 0, 0, 0, 0, 0};
  REQUIRE(seams_q4_q6(xyz, 4, zero, 3.5, q4, q6, nullptr, nullptr) == 1);
  REQUIRE(seams_q4_q6_lammps(xyz, 4, nullptr, 3.5, q4, q6, nullptr, nullptr) == 1);
}
