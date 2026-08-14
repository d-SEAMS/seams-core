#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ira_sofi.hpp>

#include <Eigen/Dense>
#include <cmath>

TEST_CASE("ira::available matches the compile-time flag", "[ira]") {
#ifdef SEAMS_HAS_IRA
  REQUIRE(ira::available());
#else
  REQUIRE_FALSE(ira::available());
  ira::Match m;
  Eigen::MatrixXd pts(3, 3);
  pts.setZero();
  REQUIRE(ira::match(pts, pts, m) == 1);
  ira::PointGroup pg;
  REQUIRE(ira::pointGroup(pts, pg) == 1);
#endif
}

#ifdef SEAMS_HAS_IRA

static Eigen::MatrixXd squareXY() {
  Eigen::MatrixXd p(4, 3);
  p << 1.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, -1.0, 0.0, 1.0, -1.0, 0.0;
  return p;
}

TEST_CASE("IRA overlays a rotated and permuted square", "[ira]") {
  const Eigen::MatrixXd ref = squareXY();
  Eigen::MatrixXd tgt(4, 3);
  // 90 deg about z, then swap first two atoms.
  tgt.row(0) = Eigen::RowVector3d(-1.0, -1.0, 0.0);
  tgt.row(1) = Eigen::RowVector3d(-1.0, 1.0, 0.0);
  tgt.row(2) = Eigen::RowVector3d(1.0, 1.0, 0.0);
  tgt.row(3) = Eigen::RowVector3d(1.0, -1.0, 0.0);
  ira::Match m;
  REQUIRE(ira::match(ref, tgt, m) == 0);
  REQUIRE(m.assignment.size() == 4);
  REQUIRE_THAT(m.rmsd, Catch::Matchers::WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(m.hausdorff, Catch::Matchers::WithinAbs(0.0, 1e-6));
}

TEST_CASE("SOFI reports a non-C1 group for a square", "[ira]") {
  ira::PointGroup pg;
  REQUIRE(ira::pointGroup(squareXY(), pg) == 0);
  REQUIRE(pg.nOperations >= 2);
  REQUIRE(pg.symbol != "C1");
  REQUIRE(pg.symbol != "C1 ");
}

#endif
