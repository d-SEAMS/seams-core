// Internal
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <seams_output.hpp>
#include <topo_one_dim.hpp>

// Standard
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <random>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <rang.hpp>

namespace fs = std::filesystem;

// Helper: build the 8-atom prism geometry (reusable)
static void buildPrismCloud(molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud.box = {10, 10, 50};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
  molSys::Point<double> iPoint;
  iPoint.type = 1;
  double coords[8][3] = {
      {0, 0, 0}, {2.75, 0, 0}, {2.75, 2.75, 0}, {0, 2.75, 0},
      {0, 0, 2.75}, {2.75, 0, 2.75}, {2.75, 2.75, 2.75}, {0, 2.75, 2.75}};
  for (int i = 0; i < 8; i++) {
    iPoint.atomID = i;
    iPoint.x = coords[i][0]; iPoint.y = coords[i][1]; iPoint.z = coords[i][2];
    yCloud.pts.push_back(iPoint);
    yCloud.idIndexMap[i] = i;
  }
  yCloud.nop = 8;
}

SCENARIO("Test the prism identification scheme for a single tetragonal prism.",
         "[topo]") {
  GIVEN("A pointCloud") {
    // Hard-coded example of a single tetragonal prism
    // Here, we manually input a tetragonal prism's coordinates
    molSys::PointCloud<molSys::Point<double>, double> yCloud; // pointCloud
    molSys::Point<double> iPoint;                             // A single point
    std::vector<double> box = {10, 10, 50};                   // Box lengths
    yCloud.box = box; // Update the box lengths. z dim is axial
    std::vector<std::vector<int>> nList; // Neighbour list
    std::vector<std::vector<int>> rings; // Rings
    std::vector<ring::strucType>
        ringType; // This vector will have a value for each ring inside
    std::vector<int> listPrism; // Contains indices of rings in each prism
    // --------------------
    // Building the prism block (Assuming an ideal tetragonal prism)
    // Each side length is 2.75
    //
    iPoint.type = 1; // Same for all the points here
    // Element {0}
    iPoint.atomID = 0; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 0;      // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {1}
    iPoint.atomID = 1; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 0;      // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {2}
    iPoint.atomID = 2; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 2.75;   // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {3}
    iPoint.atomID = 3; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 2.75;   // y
    iPoint.z = 0;      // z
    yCloud.pts.push_back(iPoint);
    // Element {4}
    iPoint.atomID = 4; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 0;      // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Element {5}
    iPoint.atomID = 5; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 0;      // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Element {6}
    iPoint.atomID = 6; // iatom
    iPoint.x = 2.75;   // x
    iPoint.y = 2.75;   // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Element {7}
    iPoint.atomID = 7; // iatom
    iPoint.x = 0;      // x
    iPoint.y = 2.75;   // y
    iPoint.z = 2.75;   // z
    yCloud.pts.push_back(iPoint);
    // Update nop
    yCloud.nop = yCloud.pts.size();
    // Update the unordered map
    for (int iatom = 0; iatom < yCloud.nop; iatom++) {
      yCloud.idIndexMap[iatom] = iatom;
    } // end of filling the map
    // --------------------
    WHEN("Given a pointCloud, and a neighbour list") {
      // Calculate a neighbour list
      nList = nneigh::neighListO(3.5, yCloud, 1);
      // Find the vector of vector of rings
      rings = primitive::ringNetwork(nList, 5);
      THEN("There should be exactly one tetragonal prism.") {
        int nPrisms = 0; // The number of perfect prisms
        int nDeformed = 0;
        // Qualifier for the type of atom it is:
        std::vector<double> rmsdPerAtom; // Not used here
        ringType.resize(
            rings.size()); // Has a value for each ring. init to zero.
        // Find the number of tetragonal prisms
        listPrism = ring::findPrisms(rings, ringType, nPrisms, nDeformed,
                                     nList, yCloud, rmsdPerAtom, false);
        // Assert the number of prism blocks
        REQUIRE(nPrisms ==
                1); // Evaluate condition for a single tetragonal prism
      }
    } // End of getting the neighbour list
  }   // End of given
} // End of scenario

// -- Additional tests for wrapper functions --

TEST_CASE("prismAnalysis integration with tetragonal prism", "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildPrismCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto rings = primitive::ringNetwork(nList, 5);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_prismanalysis/").string();
  int atomID = 0;

  int ret = ring::prismAnalysis(tmpPath, rings, nList, yCloud, 5, atomID, 1, 1,
                                 false);
  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}


TEST_CASE("relaxedPrismConditions detects prism with relaxed criteria",
          "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildPrismCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 5);

  // Find two 4-membered rings (basal rings of the tetragonal prism)
  // The prism should have basal rings
  std::vector<std::vector<int>> fourRings;
  for (const auto &r : rings) {
    if (r.size() == 4) fourRings.push_back(r);
  }

  if (fourRings.size() >= 2) {
    bool relaxed =
        ring::relaxedPrismConditions(nList, fourRings[0], fourRings[1]);
    // At least one bond between basal rings should exist
    REQUIRE(relaxed == true);
  }
}

TEST_CASE("prismAnalysis with shape matching on tetragonal prism",
          "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildPrismCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto rings = primitive::ringNetwork(nList, 5);

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_prism_sm/").string();
  int atomID = 0;

  int ret = ring::prismAnalysis(tmpPath, rings, nList, yCloud, 5, atomID, 1, 1,
                                 true);
  REQUIRE(ret == 0);

  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}

TEST_CASE("assignPrismType assigns atom types for known prism list",
          "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildPrismCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto rings = primitive::ringNetwork(nList, 5);

  std::vector<ring::strucType> ringType(rings.size());
  int nPrisms = 0, nDeformed = 0;
  std::vector<double> rmsdPerAtom;
  auto listPrism = ring::findPrisms(rings, ringType, nPrisms, nDeformed,
                                     nList, yCloud, rmsdPerAtom, false);

  // Now assign prism types
  std::vector<int> atomTypes(yCloud.nop, 1); // 1 = dummy
  std::vector<ring::strucType> atomState(yCloud.nop);

  int ret = ring::assignPrismType(rings, listPrism, 4, ringType, atomTypes,
                                   atomState);
  REQUIRE(ret == 0);

  // At least some atoms should be classified as prism type (not dummy)
  bool hasNonDummy = false;
  for (int i = 0; i < yCloud.nop; i++) {
    if (atomTypes[i] != 1) hasNonDummy = true;
  }
  REQUIRE(hasNonDummy);
}

TEST_CASE("rmAxialTranslations wraps through dump H not bound span",
          "[topo_one_dim]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {25.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  molSys::Point<double> a;
  a.x = 1.0;
  a.y = 0.5;
  a.z = 1.0;
  a.atomID = 1;
  molSys::Point<double> b;
  b.x = 0.2;
  b.y = 0.5;
  b.z = 1.0;
  b.atomID = 2;
  cloud.pts.push_back(a);
  cloud.pts.push_back(b);
  cloud.idIndexMap[1] = 0;
  cloud.idIndexMap[2] = 1;
  int atomID = 0;
  REQUIRE(ring::rmAxialTranslations(cloud, atomID, 1, 1) == 0);
  REQUIRE(atomID == 1);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(0.2886751345948129, 1e-10));
  REQUIRE_THAT(cloud.pts[0].y, Catch::Matchers::WithinAbs(0.5, 1e-12));
  REQUIRE_THAT(cloud.pts[1].x, Catch::Matchers::WithinAbs(19.488675134594813, 1e-10));
  REQUIRE_THAT(cloud.pts[1].y, Catch::Matchers::WithinAbs(0.5, 1e-12));
  REQUIRE_THAT(cloud.pts[1].z, Catch::Matchers::WithinAbs(1.0, 1e-12));
}

// findPrisms without shape matching tests only the rings that hold a
// neighbour of each first basal atom. The all-pairs loop below is the
// definition; both must give the same counts, ring types and prism list.
namespace {

molSys::PointCloud<molSys::Point<double>, double> jitteredTubes(double jitter, unsigned seed) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double side = 2.75;
  const int layers = 10;
  std::mt19937 rng(seed);
  std::normal_distribution<double> noise(0.0, jitter);
  double x0 = 0.0;
  int id = 0;
  for (int n : {4, 5, 6, 3}) {
    const double r = side / (2.0 * std::sin(M_PI / n));
    const double cx = x0 + r + 4.0;
    for (int k = 0; k < layers; ++k) {
      for (int j = 0; j < n; ++j) {
        const double a = 2.0 * M_PI * j / n;
        molSys::Point<double> p;
        p.type = 1;
        p.atomID = id;
        p.x = cx + r * std::cos(a) + noise(rng);
        p.y = 10.0 + r * std::sin(a) + noise(rng);
        p.z = std::fmod(k * side + noise(rng) + layers * side, layers * side);
        cloud.pts.push_back(p);
        cloud.idIndexMap[id] = id;
        ++id;
      }
    }
    x0 = cx + r + 4.0;
  }
  cloud.nop = id;
  cloud.box = {x0 + 2.0, 20.0, layers * side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  return cloud;
}

struct PrismResult {
  std::vector<int> list;
  std::vector<ring::strucType> types;
  int perfect = 0;
};

PrismResult allPairs(const std::vector<std::vector<int>> &rings,
                     const std::vector<std::vector<int>> &nList,
                     molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  PrismResult out;
  out.types.assign(rings.size(), ring::strucType::unclassified);
  for (std::size_t i = 0; i + 1 < rings.size(); ++i) {
    for (std::size_t j = i + 1; j < rings.size(); ++j) {
      std::vector<int> a = rings[i], b = rings[j];
      if (a.size() == 4 && !ring::discardExtraTetragonBlocks(a, b, cloud)) continue;
      if (ring::hasCommonElements(a, b)) continue;
      if (!ring::basalPrismConditions(nList, a, b)) continue;
      out.perfect += 1;
      for (std::size_t k : {i, j}) {
        if (out.types[k] == ring::strucType::unclassified) {
          out.types[k] = ring::strucType::Prism;
          out.list.push_back(static_cast<int>(k));
        }
      }
    }
  }
  std::sort(out.list.begin(), out.list.end());
  out.list.erase(std::unique(out.list.begin(), out.list.end()), out.list.end());
  return out;
}

}  // namespace

TEST_CASE("findPrisms matches the all-pairs prism search on jittered nanotubes", "[topo]") {
  int prisms = 0;
  for (double jitter : {0.0, 0.1, 0.25, 0.4}) {
    auto cloud = jitteredTubes(jitter, 11u);
    auto withSelf = nneigh::neighbourListByIndex(cloud, nneigh::neighListO(3.5, cloud, 1));
    // The C API passes rows without the atom in front.
    std::vector<std::vector<int>> noSelf;
    for (const auto &row : withSelf) noSelf.emplace_back(row.begin() + 1, row.end());
    const auto all = primitive::ringNetwork(withSelf, 6);
    for (const auto *nList : {&withSelf, &noSelf}) {
      for (std::size_t size = 3; size <= 6; ++size) {
        std::vector<std::vector<int>> rings;
        for (const auto &r : all) {
          if (r.size() == size) rings.push_back(r);
        }
        if (rings.empty()) continue;
        const PrismResult want = allPairs(rings, *nList, cloud);
        std::vector<ring::strucType> types(rings.size(), ring::strucType::unclassified);
        std::vector<double> rmsd(cloud.nop, -1.0);
        int perfect = -1, imperfect = -1;
        const auto list =
            ring::findPrisms(rings, types, perfect, imperfect, *nList, cloud, rmsd, false);
        REQUIRE(list == want.list);
        REQUIRE(types == want.types);
        REQUIRE(perfect == want.perfect);
        REQUIRE(imperfect == 0);
        prisms += perfect;
      }
    }
  }
  REQUIRE(prisms > 0);
}
