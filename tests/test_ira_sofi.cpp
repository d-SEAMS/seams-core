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

TEST_CASE("IRA residual is R times ref plus t versus assigned target",
          "[ira]") {
  // 90 deg about z is not an involution, so R*ref+t and ref-(R*target+t)
  // cannot both be ~0. The C API convention used here is the former.
  const Eigen::MatrixXd ref = squareXY();
  Eigen::MatrixXd tgt(4, 3);
  tgt.row(0) = Eigen::RowVector3d(-1.0, 1.0, 0.0);
  tgt.row(1) = Eigen::RowVector3d(-1.0, -1.0, 0.0);
  tgt.row(2) = Eigen::RowVector3d(1.0, -1.0, 0.0);
  tgt.row(3) = Eigen::RowVector3d(1.0, 1.0, 0.0);
  ira::Match m;
  REQUIRE(ira::match(ref, tgt, m) == 0);
  REQUIRE_THAT(m.rmsd, Catch::Matchers::WithinAbs(0.0, 1e-6));
  double alt = 0.0;
  int used = 0;
  const int n = 4;
  for (int i = 0; i < n; i++) {
    int j = (i < static_cast<int>(m.assignment.size()))
                ? static_cast<int>(m.assignment[static_cast<size_t>(i)])
                : -1;
    if (j < 0 || j >= n) {
      continue;
    }
    const Eigen::Vector3d a(ref(i, 0), ref(i, 1), ref(i, 2));
    const Eigen::Vector3d b(tgt(j, 0), tgt(j, 1), tgt(j, 2));
    alt += (a - (m.rotation * b + m.translation)).squaredNorm();
    used++;
  }
  if (used > 0) {
    alt = std::sqrt(alt / static_cast<double>(used));
  }
  REQUIRE(alt > 0.5);
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
