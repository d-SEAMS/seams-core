#include <catch2/catch_test_macros.hpp>

#include <mol_sys.hpp>
#include <site.hpp>

#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

// A 3x3x3 simple-cubic lattice of "water" at spacing 3.0 in a periodic box,
// with two lattice sites turned into ions (type 3).
Cloud lattice(std::vector<int> ionSites) {
  Cloud cloud;
  const double a = 3.0;
  int id = 0;
  for (int x = 0; x < 3; x++) {
    for (int y = 0; y < 3; y++) {
      for (int z = 0; z < 3; z++) {
        molSys::Point<double> p;
        p.x = x * a;
        p.y = y * a;
        p.z = z * a;
        p.type = 1;
        p.atomID = id + 1;
        p.molID = id + 1;
        for (int s : ionSites) {
          if (s == id) {
            p.type = 3;
          }
        }
        cloud.pts.push_back(p);
        cloud.idIndexMap[id + 1] = id;
        ++id;
      }
    }
  }
  cloud.nop = id;
  cloud.box = {3 * a, 3 * a, 3 * a};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

} // namespace

TEST_CASE("an ion whose shell is all ice is in ice", "[ions]") {
  const auto cloud = lattice({13});
  std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), true);
  const auto env = site::ionEnvironment(cloud, ice, {13}, 1, 3.5);
  REQUIRE(env.ion.size() == 1);
  REQUIRE(env.shell[0] == 6);
  REQUIRE(env.iceFraction[0] == 1.0);
  REQUIRE(env.state[0] == site::IonState::ice);
  REQUIRE(env.nIce == 1);
  REQUIRE(env.nFront == 0);
  REQUIRE(env.nLiquid == 0);
}

TEST_CASE("an ion with no labelled neighbour is in liquid, a mixed shell is at the front",
          "[ions]") {
  const auto cloud = lattice({0, 13});
  std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), false);
  // label half the shell of atom 13 (its +x and -x neighbours are 22 and 4)
  ice[22] = true;
  ice[4] = true;
  const auto env = site::ionEnvironment(cloud, ice, {0, 13}, 1, 3.5);
  REQUIRE(env.ion == std::vector<int>{0, 13});
  REQUIRE(env.state[0] == site::IonState::liquid);
  REQUIRE(env.state[1] == site::IonState::front);
  REQUIRE(env.iceFraction[1] == 2.0 / 6.0);
  REQUIRE(env.nLiquid == 1);
  REQUIRE(env.nFront == 1);
}

TEST_CASE("ice cluster census assigns an ion to the nearest ice component",
          "[ions]") {
  const auto cloud = lattice({13});
  std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), false);
  // 4 and 22 are the +/-x shell of the ion at 13; 0 is a far ice atom.
  ice[4] = true;
  ice[22] = true;
  ice[0] = true;
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(cloud.nop));
  const int n = cloud.nop;
  auto add = [&](int a, int b) {
    rows[static_cast<std::size_t>(a)].push_back(b);
    rows[static_cast<std::size_t>(b)].push_back(a);
  };
  for (int i = 0; i < n; i++) {
    rows[static_cast<std::size_t>(i)].push_back(i);
  }
  add(4, 22);
  const auto census = site::iceClusterIonCensus(cloud, ice, rows, {13}, 3.5);
  REQUIRE(census.nClusters == 2);
  REQUIRE(census.clusterOf[4] == census.clusterOf[22]);
  REQUIRE(census.clusterOf[4] != census.clusterOf[0]);
  int assigned = 0;
  for (int nIon : census.ionsInCluster) {
    assigned += nIon;
  }
  REQUIRE(assigned == 1);
  REQUIRE(census.clusterOfIon[0] == census.clusterOf[4]);
}

TEST_CASE("ions do not count in each other's shells and the periodic image is seen",
          "[ions]") {
  // 0 at the origin and 2 at (0, 0, 6): under the 9.0 box they sit 3.0 apart
  // through the boundary, so the ion at 0 has five water neighbours, not six
  const auto cloud = lattice({0, 2});
  std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), true);
  const auto env = site::ionEnvironment(cloud, ice, {0, 2}, 1, 3.5);
  REQUIRE(env.shell[0] == 5);
  REQUIRE(env.shell[1] == 5);
  REQUIRE(env.state[0] == site::IonState::ice);
}

TEST_CASE("waterType zero accepts every non-ion atom", "[ions]") {
  auto cloud = lattice({13});
  cloud.pts[4].type = 2;  // a second water species in the shell
  std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), true);
  REQUIRE(site::ionEnvironment(cloud, ice, {13}, 1, 3.5).shell[0] == 5);
  REQUIRE(site::ionEnvironment(cloud, ice, {13}, 0, 3.5).shell[0] == 6);
}
