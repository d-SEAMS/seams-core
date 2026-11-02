#include <catch2/catch_test_macros.hpp>

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <site.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

// Rows by index with a leading self entry, from an edge list.
std::vector<std::vector<int>> rowsOf(int n, const std::vector<std::pair<int, int>> &edges) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    rows[static_cast<std::size_t>(i)].push_back(i);
  }
  for (const auto &[a, b] : edges) {
    rows[static_cast<std::size_t>(a)].push_back(b);
    rows[static_cast<std::size_t>(b)].push_back(a);
  }
  return rows;
}

std::vector<std::vector<int>> sortedRings(std::vector<std::vector<int>> rings) {
  for (auto &r : rings) {
    std::sort(r.begin(), r.end());
  }
  std::sort(rings.begin(), rings.end());
  return rings;
}

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

TEST_CASE("capped and broken rings on a hand-built graph", "[ions][rings]") {
  // Hexagon 0-1-2-3-4-5, square 5-6-7-8 sharing vertex 5, pendant 8-9.
  // Shell: every vertex but 7 and 8.
  const auto rows = rowsOf(10, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
                                {5, 6}, {6, 7}, {7, 8}, {8, 5}, {8, 9}});
  const std::vector<int> shell = {0, 1, 2, 3, 4, 5, 6, 9};
  const auto rings = primitive::ringNetwork(rows, 6);
  REQUIRE(sortedRings(rings) ==
          std::vector<std::vector<int>>{{0, 1, 2, 3, 4, 5}, {5, 6, 7, 8}});

  // capped: the hexagon lies in the shell; the square has 7 and 8 outside it
  REQUIRE(sortedRings(site::cappedRings(rings, shell, 6)) ==
          std::vector<std::vector<int>>{{0, 1, 2, 3, 4, 5}});
  REQUIRE(site::cappedRings(rings, shell, 5).empty());
  // the shell census counts both rings, capped or not
  const auto census = site::shellRingCensus(rings, shell, 6);
  REQUIRE(census[4] == 1);
  REQUIRE(census[6] == 1);

  // broken, with the virtual ion as vertex 10 bonded to the eight shell
  // waters. A primitive ring through the ion has no other shell vertex but
  // its two ring neighbours (each would be a chord to the ion), so:
  //   seven triangles, one per shell-shell bond (six hexagon bonds and 5-6);
  //   the square 10-5-8-9 through the non-shell 8;
  //   the pentagon 10-6-7-8-9 (10-6-7-8-5 has the chord 5-6).
  const int ion = 10;
  const auto broken = site::brokenRings(rows, shell, 6);
  for (const auto &r : broken) {
    REQUIRE(std::count(r.begin(), r.end(), ion) == 1);
  }
  REQUIRE(sortedRings(broken) ==
          std::vector<std::vector<int>>{{0, 1, 10}, {0, 5, 10}, {1, 2, 10},
                                        {2, 3, 10}, {3, 4, 10}, {4, 5, 10},
                                        {5, 6, 10}, {5, 8, 9, 10}, {6, 7, 8, 9, 10}});
  std::map<std::size_t, int> bySize;
  for (const auto &r : broken) {
    ++bySize[r.size()];
  }
  REQUIRE(bySize == std::map<std::size_t, int>{{3, 7}, {4, 1}, {5, 1}});
  // the ring-size limit applies to the augmented graph too
  REQUIRE(site::brokenRings(rows, shell, 4).size() == 8);
  REQUIRE(site::brokenRings(rows, shell, 3).size() == 7);
  REQUIRE(site::brokenRings(rows, {}, 6).empty());
}

TEST_CASE("the broken rings of a lattice ion are the rings its removal lost",
          "[ions][rings]") {
  const int ion = 62;
  const auto cloud = lattice(ion);
  const auto rows = nneigh::neighbourListByIndex(cloud, nneigh::neighListO(3.5, cloud, 1));
  const auto rings = primitive::ringNetwork(rows, 4);
  const std::vector<bool> ice(static_cast<std::size_t>(cloud.nop), true);
  const auto env = site::ionEnvironment(cloud, ice, {ion}, 1, 3.5);
  // no two shell waters are bonded, so no ring lies wholly in the shell
  REQUIRE(site::cappedRings(rings, env.members[0], 4).empty());
  // the virtual ion restores the twelve squares through the lattice site
  const auto broken = site::brokenRings(rows, env.members[0], 4);
  REQUIRE(broken.size() == 375 - rings.size());
  REQUIRE(broken.size() == 12);
  const int virtualIon = static_cast<int>(rows.size());
  for (const auto &r : broken) {
    REQUIRE(r.size() == 4);
    REQUIRE(std::count(r.begin(), r.end(), virtualIon) == 1);
  }
}
