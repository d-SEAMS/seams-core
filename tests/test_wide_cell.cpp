#include <catch2/catch_test_macros.hpp>

#include <local_symmetry.hpp>
#include <mol_sys.hpp>

#include <vector>

namespace {

molSys::PointCloud<molSys::Point<double>, double> bondedPair(double side) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  for (int i = 0; i < 2; ++i) {
    molSys::Point<double> point;
    point.x = static_cast<double>(i) * 2.75;
    point.y = 0.0;
    point.z = 0.0;
    point.type = 1;
    point.atomID = i + 1;
    point.molID = i + 1;
    cloud.idIndexMap[point.atomID] = i;
    cloud.pts.push_back(point);
  }
  cloud.nop = 2;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

} // namespace

TEST_CASE("a bonded pair in a wide cell matches the small-cell labels", "[shell]") {
  // ceil(4516 / 3.5) is 1291, and 1291^3 does not fit in a signed 32-bit cell index.
  std::vector<int> small;
  std::vector<int> wide;
  symm::classifyIceCloud(bondedPair(20.0), small);
  symm::classifyIceCloud(bondedPair(4516.0), wide);
  REQUIRE(small.size() == 2);
  REQUIRE(wide == small);
}
