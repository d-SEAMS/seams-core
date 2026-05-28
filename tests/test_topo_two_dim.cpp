#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <topo_two_dim.hpp>

#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

// Helper: build a cloud with atoms forming a simple ring structure
// Creates 6 atoms in a hexagonal arrangement
static molSys::PointCloud<molSys::Point<double>, double>
makeHexCloud(double radius = 2.0) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  for (int i = 0; i < 6; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i;
    double angle = 2.0 * M_PI * i / 6.0;
    pt.x = 10.0 + radius * std::cos(angle);
    pt.y = 10.0 + radius * std::sin(angle);
    pt.z = 10.0;
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 6;
  return cloud;
}

TEST_CASE("polygonRingAnalysis runs with empty rings", "[topo_two_dim]") {
  auto cloud = makeHexCloud();
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 2.5);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topo2d/").string();
  fs::create_directories(tmpPath);

  std::vector<std::vector<int>> rings; // empty
  int ret = ring::polygonRingAnalysis(tmpPath, rings, nList, cloud, 6, 100.0, 1);

  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("polygonRingAnalysis processes known ring set", "[topo_two_dim]") {
  auto cloud = makeHexCloud();
  auto nList = nneigh::getNewNeighbourListByIndex(cloud, 2.5);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_topo2d_rings/").string();
  fs::create_directories(tmpPath);

  // A single 6-membered ring
  std::vector<std::vector<int>> rings = {{0, 1, 2, 3, 4, 5}};

  int ret = ring::polygonRingAnalysis(tmpPath, rings, nList, cloud, 6, 100.0, 1);

  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}
