#include <catch2/catch_test_macros.hpp>

#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <vector>

TEST_CASE("a complete cutoff list past 32-bit offsets is refused", "[neighbours]") {
  // n*(n-1) is past INT_MAX. The cutoff is past the longest minimum image
  // of a 10 A cube, so every pair would be stored.
  constexpr int n = 46342;
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.pts.resize(static_cast<std::size_t>(n));
  cloud.nop = n;
  for (int i = 0; i < n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.type = 1;
    point.atomID = i + 1;
    point.molID = i + 1;
    point.x = 0.0;
    point.y = 0.0;
    point.z = 0.0;
  }
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  std::vector<int> offsets;
  std::vector<int> cols;
  REQUIRE_FALSE(nneigh::orthoFullCsr(cloud, nullptr, 20.0, offsets, cols, false));
}
