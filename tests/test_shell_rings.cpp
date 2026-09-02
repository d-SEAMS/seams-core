#include <catch2/catch_test_macros.hpp>

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <site.hpp>

#include <numeric>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

// A 5x5x5 simple-cubic lattice at spacing 3.0, periodic, one site an ion.
// With the ring depth capped at four the only rings are the square faces:
// 125 * 12 / 4 = 375 of them, twelve through any one site.
Cloud lattice(int ionSite) {
  Cloud cloud;
  const double a = 3.0;
  int id = 0;
  for (int x = 0; x < 5; x++) {
    for (int y = 0; y < 5; y++) {
      for (int z = 0; z < 5; z++) {
        molSys::Point<double> p;
        p.x = x * a;
        p.y = y * a;
        p.z = z * a;
        p.type = id == ionSite ? 3 : 1;
        p.atomID = id + 1;
        p.molID = id + 1;
        cloud.pts.push_back(p);
        cloud.idIndexMap[id + 1] = id;
        ++id;
      }
    }
  }
  cloud.nop = id;
  cloud.box = {5 * a, 5 * a, 5 * a};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

} // namespace

TEST_CASE("the rings through an ion's shell are the squares that survive it", "[ions][rings]") {
  const int ion = 62;  // (2, 2, 2), the centre
  const auto cloud = lattice(ion);
  const auto rows = nneigh::neighbourListByIndex(cloud, nneigh::neighListO(3.5, cloud, 1));
  const auto rings = primitive::ringNetwork(rows, 4);
  // the ion is not a vertex: its twelve squares are gone
  REQUIRE(rings.size() == 375 - 12);
  const std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), true);
  const auto env = site::ionEnvironment(cloud, ice, {ion}, 1, 3.5);
  REQUIRE(env.members.size() == 1);
  REQUIRE(env.members[0].size() == 6);
  REQUIRE(env.shell[0] == 6);
  const auto census = site::shellRingCensus(rings, env.members[0], 4);
  REQUIRE(census.size() == 5);
  // each shell molecule keeps 12 - 4 squares; no surviving square holds two
  // shell molecules, since the square through any two of them ran through the ion
  REQUIRE(census[4] == 6 * 8);
  REQUIRE(census[3] == 0);
  // the whole census over every molecule is every ring
  std::vector<int> all(static_cast<std::size_t>(cloud.nop));
  std::iota(all.begin(), all.end(), 0);
  REQUIRE(site::shellRingCensus(rings, all, 4)[4] == static_cast<int>(rings.size()));
  // a cap below the ring size counts nothing
  REQUIRE(site::shellRingCensus(rings, env.members[0], 3) == std::vector<int>{0, 0, 0, 0});
}
