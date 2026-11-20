#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <topo_one_dim.hpp>

#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>

std::vector<std::vector<int>> keepAxialRingsOnly(
    const std::vector<std::vector<int>> &rings,
    molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  const int axialDim = nneigh::dumpAxialDim(cloud);
  auto axialCoordinate = [axialDim](const molSys::Point<double> &point) {
    if (axialDim == 0) {
      return point.x;
    }
    if (axialDim == 1) {
      return point.y;
    }
    return point.z;
  };

  std::vector<std::vector<int>> axialRings;
  for (const auto &ring : rings) {
    const double plane = axialCoordinate(cloud.pts[static_cast<std::size_t>(ring.front())]);
    bool perpendicularToAxis = true;
    for (const int atom : ring) {
      if (std::abs(axialCoordinate(cloud.pts[static_cast<std::size_t>(atom)]) - plane) >
          1e-12) {
        perpendicularToAxis = false;
        break;
      }
    }
    if (perpendicularToAxis) {
      axialRings.push_back(ring);
    }
  }
  return axialRings;
}

using namespace ring;

TEST_CASE("consumer axial-ring filter is not shadowed by the internal header",
          "[keep_axial_name]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 50.0};
  cloud.boxLow = {0.0, 0.0, 0.0};

  const double coordinates[8][3] = {
      {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {2.0, 2.0, 0.0}, {0.0, 2.0, 0.0},
      {0.0, 0.0, 3.0}, {2.0, 0.0, 3.0}, {2.0, 2.0, 3.0}, {0.0, 2.0, 3.0},
  };
  for (int atom = 0; atom < 8; ++atom) {
    molSys::Point<double> point;
    point.atomID = atom;
    point.x = coordinates[atom][0];
    point.y = coordinates[atom][1];
    point.z = coordinates[atom][2];
    cloud.pts.push_back(point);
    cloud.idIndexMap[atom] = atom;
  }
  cloud.nop = static_cast<int>(cloud.pts.size());

  const std::vector<std::vector<int>> rings = {
      {0, 1, 2, 3},
      {4, 5, 6, 7},
      {0, 1, 5, 4},
  };
  const std::vector<std::vector<int>> expected = {
      {0, 1, 2, 3},
      {4, 5, 6, 7},
  };

  REQUIRE(keepAxialRingsOnly(rings, cloud) == expected);
}
