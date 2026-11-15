#include <cage.hpp>
#include <cage_canon.hpp>
#include <cage_enum.hpp>
#include <ira_sofi.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

namespace {

std::set<std::vector<int>>
vertexSetsOf(const std::vector<cage::Cage> &cages,
             cage::cageType type, const std::vector<std::vector<int>> &rings) {
  std::set<std::vector<int>> out;
  for (const auto &c : cages) {
    if (c.type != type) {
      continue;
    }
    std::set<int> verts;
    for (const int ri : c.rings) {
      if (ri < 0 || ri >= static_cast<int>(rings.size())) {
        continue;
      }
      for (const int a : rings[static_cast<size_t>(ri)]) {
        verts.insert(a);
      }
    }
    out.emplace(verts.begin(), verts.end());
  }
  return out;
}

std::set<std::vector<int>>
vertexSetsOf(const std::vector<cage::FoundCage> &cages) {
  std::set<std::vector<int>> out;
  for (const auto &c : cages) {
    out.insert(c.vertices);
  }
  return out;
}

std::vector<std::vector<int>> hexPrismRings() {
  return {{0, 1, 2, 3, 4, 5},
          {6, 7, 8, 9, 10, 11},
          {0, 1, 7, 6},
          {1, 2, 8, 7},
          {2, 3, 9, 8},
          {3, 4, 10, 9},
          {4, 5, 11, 10},
          {5, 0, 6, 11}};
}

void buildHCCloud(molSys::PointCloud<molSys::Point<double>, double> &yCloud) {
  yCloud.box = {53.9690018, 54.5289994, 51.257};
  yCloud.boxLow = {0.0, 0.0, 0.0};
  yCloud.currentFrame = 1;
  molSys::Point<double> iPoint;
  iPoint.type = 1;
  const double coords[12][3] = {
      {8.995, 10.3859997, 15.0939999}, {6.7459998, 14.2810001, 15.0939999},
      {4.4970002, 10.3859997, 15.0939999}, {8.995, 12.9829998, 14.1949997},
      {6.7459998, 9.0880003, 14.1949997}, {4.4970002, 12.9829998, 14.1949997},
      {8.995, 12.9829998, 11.4329996}, {6.7459998, 9.0880003, 11.4329996},
      {4.4970002, 12.9829998, 11.4329996}, {8.995, 10.3859997, 10.5340004},
      {6.7459998, 14.2810001, 10.5340004}, {4.4970002, 10.3859997, 10.5340004}};
  for (int i = 0; i < 12; i++) {
    iPoint.atomID = i;
    iPoint.x = coords[i][0];
    iPoint.y = coords[i][1];
    iPoint.z = coords[i][2];
    yCloud.pts.push_back(iPoint);
    yCloud.idIndexMap[i] = i;
  }
  yCloud.nop = 12;
}

std::vector<std::vector<int>> sixOf(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return six;
}

} // namespace

TEST_CASE("Signature parse accepts lists and the named table", "[cage_enum]") {
  const auto sod = cage::Signature::parse("sodalite");
  REQUIRE(sod == cage::Signature::parse("4:6,6:8"));
  REQUIRE(sod.str() == "4:6,6:8");
  REQUIRE(sod.faceCount() == 14);
  REQUIRE(sod.maxRingSize() == 6);
  REQUIRE(sod.containsSize(4));
  REQUIRE(sod.containsSize(6));
  REQUIRE_FALSE(sod.containsSize(8));

  REQUIRE(cage::Signature::parse("alpha") ==
          cage::Signature::parse("4:6,6:8,8:6"));
  REQUIRE(cage::Signature::parse("512") == cage::Signature::parse("5:12"));
  REQUIRE(cage::Signature::parse("51262") ==
          cage::Signature::parse("5:12,6:2"));
  REQUIRE(cage::Signature::parse("51264") ==
          cage::Signature::parse("5:12,6:4"));
  REQUIRE(cage::Signature::parse("51268") ==
          cage::Signature::parse("5:12,6:8"));
  REQUIRE(cage::Signature::parse("sh") ==
          cage::Signature::parse("4:3,5:6,6:3"));
  REQUIRE(cage::Signature::parse("sH") == cage::Signature::parse("sh"));
  REQUIRE(cage::Signature::parse("hc").counts ==
          cage::Signature::parse("4:6,6:2").counts);
  REQUIRE(cage::Signature::parse("hc").kind == cage::Signature::Kind::HexC);
  REQUIRE(cage::Signature::parse("4:6,6:2").kind ==
          cage::Signature::Kind::Census);
  REQUIRE(cage::Signature::parse("hc") != cage::Signature::parse("4:6,6:2"));
  REQUIRE(cage::Signature::parse("DDC").kind ==
          cage::Signature::Kind::DoubleDiaC);
  REQUIRE(cage::Signature::parse("DDC").counts ==
          cage::Signature::parse("6:7").counts);
  REQUIRE(cage::Signature::parse(" 4:6, 6:8 ") == sod);
  REQUIRE(cage::Signature::parse("4:3,4:3") == cage::Signature::parse("4:6"));
}

TEST_CASE("Signature parse rejects empty and non-positive tokens",
          "[cage_enum]") {
  REQUIRE_THROWS_AS(cage::Signature::parse(""), std::invalid_argument);
  REQUIRE_THROWS_AS(cage::Signature::parse("nope"), std::invalid_argument);
  REQUIRE_THROWS_AS(cage::Signature::parse("4:0"), std::invalid_argument);
  REQUIRE_THROWS_AS(cage::Signature::parse("0:6"), std::invalid_argument);
  REQUIRE_THROWS_AS(cage::Signature::parse("4"), std::invalid_argument);
}

TEST_CASE("hexagonal prism rings close and match the hc signature",
          "[cage_enum]") {
  const auto rings = hexPrismRings();
  std::vector<int> all(rings.size());
  for (size_t i = 0; i < rings.size(); ++i) {
    all[i] = static_cast<int>(i);
  }
  REQUIRE(cage::isClosedPolyhedron(rings, all));
  REQUIRE_FALSE(cage::isClosedPolyhedron(rings, {0, 1}));

  const auto found =
      cage::findBySignature(rings, cage::Signature::parse("4:6,6:2"));
  REQUIRE(found.size() == 1);
  REQUIRE(found[0].vertices.size() == 12);
  REQUIRE(found[0].faces.size() == 8);
  REQUIRE(found[0].signature == cage::Signature::parse("4:6,6:2"));
#ifdef SEAMS_HAS_NAUTY
  REQUIRE(cage::isHexagonalPrism(
      {rings[static_cast<size_t>(found[0].faces[0])],
       rings[static_cast<size_t>(found[0].faces[1])],
       rings[static_cast<size_t>(found[0].faces[2])],
       rings[static_cast<size_t>(found[0].faces[3])],
       rings[static_cast<size_t>(found[0].faces[4])],
       rings[static_cast<size_t>(found[0].faces[5])],
       rings[static_cast<size_t>(found[0].faces[6])],
       rings[static_cast<size_t>(found[0].faces[7])]}));
  REQUIRE_FALSE(found[0].certificate.empty());
#endif
}

TEST_CASE("a prism missing one face is an incomplete cage, not a closed one",
          "[cage_enum]") {
  auto rings = hexPrismRings();
  rings.pop_back();
  const auto sig = cage::Signature::parse("4:6,6:2");
  const auto closed = cage::findBySignature(rings, sig);
  REQUIRE(closed.empty());
  const auto cups = cage::findIncompleteBySignature(rings, sig, 6);
  REQUIRE_FALSE(cups.empty());
  REQUIRE_FALSE(cups[0].closed);
  REQUIRE(cups[0].danglingEdges > 0);
  REQUIRE(cups[0].faces.size() >= 6);
  const auto full = hexPrismRings();
  const auto stillClosed = cage::findBySignature(full, sig);
  REQUIRE(stillClosed.size() == 1);
  REQUIRE(stillClosed[0].closed);
  const auto noCups = cage::findIncompleteBySignature(full, sig, 6);
  REQUIRE(noCups.empty());
}

TEST_CASE("a full-census prism with an extra attaching quad is incomplete",
          "[cage_enum]") {
  auto rings = hexPrismRings();
  rings.pop_back();
  rings.push_back({5, 0, 12, 13});
  const auto sig = cage::Signature::parse("4:6,6:2");
  REQUIRE(sig.faceCount() == 8);
  const auto closed = cage::findBySignature(rings, sig);
  REQUIRE(closed.empty());
  const auto cups = cage::findIncompleteBySignature(rings, sig, 6);
  REQUIRE_FALSE(cups.empty());
  bool sawFull = false;
  for (const auto &c : cups) {
    REQUIRE_FALSE(c.closed);
    REQUIRE(c.danglingEdges > 0);
    if (c.faces.size() == 8) {
      sawFull = true;
    }
  }
  REQUIRE(sawFull);
}

TEST_CASE("a half-prism cup is found at the default minFaces floor",
          "[cage_enum]") {
  auto rings = hexPrismRings();
  rings.erase(rings.begin() + 4, rings.end());
  REQUIRE(rings.size() == 4);
  const auto sig = cage::Signature::parse("4:6,6:2");
  REQUIRE(sig.faceCount() / 2 == 4);
  const auto closed = cage::findBySignature(rings, sig);
  REQUIRE(closed.empty());
  const auto cups = cage::findIncompleteBySignature(rings, sig, 0);
  REQUIRE_FALSE(cups.empty());
  REQUIRE_FALSE(cups[0].closed);
  REQUIRE(cups[0].faces.size() >= 4);
  REQUIRE(cups[0].danglingEdges > 0);
}

TEST_CASE("findBySignature(hc) matches findHC vertices on the twelve-atom HC",
          "[cage_enum]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);
  auto six = sixOf(rings);
  std::map<int, int> census;
  for (const auto &r : rings) {
    census[static_cast<int>(r.size())] += 1;
  }
  std::string censusText;
  for (const auto &kv : census) {
    censusText += " size" + std::to_string(kv.first) + "=" +
                  std::to_string(kv.second);
  }
  INFO("ring census" << censusText << " nRings=" << rings.size()
                     << " nSix=" << six.size());

  std::vector<ring::strucType> ringType(six.size());
  std::vector<cage::Cage> cageList;
  ring::findHC(six, ringType, nList, cageList);

  const auto found =
      cage::findBySignature(rings, nList, cage::Signature::parse("hc"));
  INFO("found=" << found.size() << " findHC=" << cageList.size());
  REQUIRE(vertexSetsOf(found) ==
          vertexSetsOf(cageList, cage::cageType::HexC, six));
  REQUIRE(found.size() == 1);
  REQUIRE(found[0].vertices.size() == 12);
}

TEST_CASE("findBySignature(hc) and findBySignature(ddc) match findHC/findDDC "
          "on mW cubic",
          "[cage_enum]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);

  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nList, 7);
  auto six = sixOf(rings);

  std::vector<ring::strucType> ringType(six.size());
  std::vector<cage::Cage> cageList;
  auto listHC = ring::findHC(six, ringType, nList, cageList);
  ring::findDDC(six, ringType, listHC, cageList);

  const auto hc =
      cage::findBySignature(rings, nList, cage::Signature::parse("hc"));
  const auto ddc =
      cage::findBySignature(rings, nList, cage::Signature::parse("ddc"));

  REQUIRE(vertexSetsOf(hc) ==
          vertexSetsOf(cageList, cage::cageType::HexC, six));
  REQUIRE(vertexSetsOf(ddc) ==
          vertexSetsOf(cageList, cage::cageType::DoubleDiaC, six));
  // bench cages_v2 deposited 2040. That was the unpatched conflict set.
  // The maximum packing is one cage for every two molecules.
  REQUIRE(yCloud.nop == 4096);
  REQUIRE(ddc.size() == 2048);
}

struct DumpCages {
  std::vector<std::vector<int>> rings;
  std::vector<cage::FoundCage> found;
};

static DumpCages cagesOnDump(const char *path, const char *spec) {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO(path, 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  const auto sig = cage::Signature::parse(spec);
  DumpCages out;
  out.rings = primitive::ringNetwork(nList, std::max(sig.maxRingSize(), 6));
  out.found = cage::findBySignature(out.rings, nList, sig);
  return out;
}

TEST_CASE("sodalite signature finds 24-vertex cages on GenIce SOD",
          "[cage_enum]") {
  const auto run = cagesOnDump("traj/genice_sod.lammpstrj", "sodalite");
  REQUIRE_FALSE(run.found.empty());
  for (const auto &c : run.found) {
    REQUIRE(c.vertices.size() == 24);
    REQUIRE(c.faces.size() == 14);
    REQUIRE(cage::isClosedPolyhedron(run.rings, c.faces));
  }
}

TEST_CASE("sodalite signature finds 24-vertex cages on GenIce FAU",
          "[cage_enum]") {
  const auto run = cagesOnDump("traj/genice_fau.lammpstrj", "sodalite");
  REQUIRE_FALSE(run.found.empty());
  for (const auto &c : run.found) {
    REQUIRE(c.vertices.size() == 24);
    REQUIRE(c.faces.size() == 14);
    REQUIRE(cage::isClosedPolyhedron(run.rings, c.faces));
  }
}

TEST_CASE("512 and 51262 signatures find hydrate cages on GenIce sI",
          "[cage_enum]") {
  const auto dodeca = cagesOnDump("traj/genice_sI.lammpstrj", "512");
  REQUIRE_FALSE(dodeca.found.empty());
  for (const auto &c : dodeca.found) {
    REQUIRE(c.vertices.size() == 20);
    REQUIRE(c.faces.size() == 12);
    REQUIRE(cage::isClosedPolyhedron(dodeca.rings, c.faces));
  }
  const auto tetra = cagesOnDump("traj/genice_sI.lammpstrj", "51262");
  REQUIRE_FALSE(tetra.found.empty());
  for (const auto &c : tetra.found) {
    REQUIRE(c.vertices.size() == 24);
    REQUIRE(c.faces.size() == 14);
    REQUIRE(cage::isClosedPolyhedron(tetra.rings, c.faces));
  }
}

TEST_CASE("perfect sI keeps closed 512 counts; leftover cups stay open",
          "[cage_enum]") {
  const auto run = cagesOnDump("traj/genice_sI.lammpstrj", "512");
  const auto sig = cage::Signature::parse("512");
  const auto nClosed = run.found.size();
  REQUIRE(cage::findBySignature(run.rings, sig).size() == nClosed);
  REQUIRE_FALSE(run.found.empty());
  const auto cups = cage::findIncompleteBySignature(run.rings, sig, 0);
  for (const auto &c : cups) {
    REQUIRE_FALSE(c.closed);
    REQUIRE(c.danglingEdges > 0);
  }
}

TEST_CASE("a nucleation dump reports incomplete cages", "[cage_enum]") {
  const auto nuc = cagesOnDump("traj/exampleTraj.lammpstrj", "hc");
  const auto sig = cage::Signature::parse("hc");
  const auto cups = cage::findIncompleteBySignature(nuc.rings, sig, 0);
  REQUIRE_FALSE(cups.empty());
  for (const auto &c : cups) {
    REQUIRE_FALSE(c.closed);
    REQUIRE(c.danglingEdges > 0);
  }
  const auto si = cagesOnDump("traj/genice_sI.lammpstrj", "512");
  const auto nClosed = cage::findBySignature(si.rings, cage::Signature::parse("512")).size();
  REQUIRE(nClosed == si.found.size());
  REQUIRE(nClosed > 0);
}

TEST_CASE("512 signature finds 20-vertex cages on GenIce sII", "[cage_enum]") {
  const auto run = cagesOnDump("traj/genice_sII.lammpstrj", "512");
  REQUIRE_FALSE(run.found.empty());
  for (const auto &c : run.found) {
    REQUIRE(c.vertices.size() == 20);
    REQUIRE(c.faces.size() == 12);
    REQUIRE(cage::isClosedPolyhedron(run.rings, c.faces));
  }
}

TEST_CASE("51264 signature finds 28-vertex cages on GenIce sII", "[cage_enum]") {
  const auto run = cagesOnDump("traj/genice_sII.lammpstrj", "51264");
  REQUIRE_FALSE(run.found.empty());
  for (const auto &c : run.found) {
    REQUIRE(c.vertices.size() == 28);
    REQUIRE(c.faces.size() == 16);
    REQUIRE(cage::isClosedPolyhedron(run.rings, c.faces));
  }
}

TEST_CASE("sH and 51268 signatures find cages on GenIce sH", "[cage_enum]") {
  const auto medium = cagesOnDump("traj/genice_sH.lammpstrj", "sH");
  REQUIRE_FALSE(medium.found.empty());
  for (const auto &c : medium.found) {
    REQUIRE(c.vertices.size() == 20);
    REQUIRE(c.faces.size() == 12);
    REQUIRE(c.signature == cage::Signature::parse("4:3,5:6,6:3"));
    REQUIRE(cage::isClosedPolyhedron(medium.rings, c.faces));
  }
  const auto large = cagesOnDump("traj/genice_sH.lammpstrj", "51268");
  REQUIRE_FALSE(large.found.empty());
  for (const auto &c : large.found) {
    REQUIRE(c.vertices.size() == 36);
    REQUIRE(c.faces.size() == 20);
    REQUIRE(c.signature == cage::Signature::parse("5:12,6:8"));
    REQUIRE(cage::isClosedPolyhedron(large.rings, c.faces));
  }
}

TEST_CASE("a species cycle matches under rotation and reversal", "[cage_enum]") {
  REQUIRE(cage::speciesCycleMatches({1, 2, 1, 2}, {2, 1, 2, 1}));
  REQUIRE(cage::speciesCycleMatches({1, 2, 1, 2}, {2, 1, 2, 1}));
  REQUIRE(cage::speciesCycleMatches({1, 2, 3}, {1, 3, 2}));
  REQUIRE_FALSE(cage::speciesCycleMatches({1, 1, 2, 2}, {1, 2, 1, 2}));
}

TEST_CASE("alternating faces close a cube and a homopolar face does not",
          "[cage_enum]") {
  const std::vector<std::vector<int>> faces = {
      {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
      {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  const auto sig = cage::Signature::parse("4:6");
  const std::vector<int> alternating = {1, 2, 1, 2, 2, 1, 2, 1};
  const std::vector<std::vector<int>> pattern = {{1, 2, 1, 2}};
  const auto closed =
      cage::findBySignature(faces, sig, alternating, pattern);
  REQUIRE(closed.size() == 1);
  REQUIRE(closed[0].closed);
  REQUIRE(closed[0].vertices.size() == 8);
  REQUIRE(cage::isClosedPolyhedron(faces, closed[0].faces));

  std::vector<int> homopolar = alternating;
  homopolar[6] = 1;
  const auto rejected =
      cage::findBySignature(faces, sig, homopolar, pattern);
  REQUIRE(rejected.empty());

  const auto uncolored = cage::findBySignature(faces, sig);
  REQUIRE(uncolored.size() == 1);
}

TEST_CASE("a homopolar edge is counted on the atoms it joins", "[cage_enum]") {
  const std::vector<std::vector<int>> nList = {{1, 3}, {0, 2}, {1, 3}, {0, 2}};
  const std::vector<std::vector<int>> rings = {{0, 1, 2, 3}};
  const std::vector<int> species = {1, 1, 2, 2};
  const auto rows = cage::formerRows(nList, rings, species, -1);
  REQUIRE(rows.size() == 4);
  REQUIRE(rows[0].coord == 2);
  REQUIRE(rows[0].homopolar == 1);
  REQUIRE(rows[0].rings.at(4) == 1);
  REQUIRE(rows[2].homopolar == 1);
  const auto formers = cage::formerRows(nList, rings, species, 1);
  REQUIRE(formers.size() == 2);
  REQUIRE(formers[0].index == 0);
  REQUIRE(formers[1].index == 1);
  REQUIRE(cage::sameNetwork(rows, rows));
  std::vector<int> split = {1, 2, 2, 2};
  const auto late = cage::formerRows(nList, rings, split, -1);
  REQUIRE(late[0].homopolar == 0);
  REQUIRE_FALSE(cage::sameNetwork(rows, late));
}

TEST_CASE("coordination skips the leading self index", "[cage_enum]") {
  // neighbourListByIndex rows are [self, neighbours...].
  const std::vector<std::vector<int>> nList = {
      {0, 1, 3}, {1, 0, 2}, {2, 1, 3}, {3, 0, 2}};
  const std::vector<std::vector<int>> rings = {{0, 1, 2, 3}};
  const std::vector<int> species = {1, 1, 2, 2};
  const auto rows = cage::formerRows(nList, rings, species, -1);
  REQUIRE(rows[0].coord == 2);
  REQUIRE(rows[0].homopolar == 1);
  REQUIRE(rows[1].coord == 2);
  REQUIRE(rows[1].homopolar == 1);
  REQUIRE(rows[2].homopolar == 1);
  REQUIRE(rows[3].homopolar == 1);
}

TEST_CASE("IRA and SOFI see the cage vertices and not the frame", "[cage_enum]") {
  std::vector<std::array<double, 3>> all(10);
  for (int i = 0; i < 8; ++i) {
    const int x = i & 1;
    const int y = (i >> 1) & 1;
    const int z = (i >> 2) & 1;
    all[static_cast<std::size_t>(i)] = {static_cast<double>(x),
                                        static_cast<double>(y),
                                        static_cast<double>(z)};
  }
  all[8] = {99.0, 99.0, 99.0};
  all[9] = {-99.0, -99.0, -99.0};
  const std::vector<int> vertices = {0, 1, 2, 3, 4, 5, 6, 7};
  const auto cageXyz = cage::coordsOfVertices(all, vertices);
  REQUIRE(cageXyz.size() == 8);
  for (const auto &p : cageXyz) {
    REQUIRE(p[0] > -2.0);
    REQUIRE(p[0] < 2.0);
    REQUIRE(p[1] > -2.0);
    REQUIRE(p[1] < 2.0);
    REQUIRE(p[2] > -2.0);
    REQUIRE(p[2] < 2.0);
  }
  const auto shape = cage::shapeOfVertices(cageXyz);
  REQUIRE(shape.nVertices == 8);
  const auto mismatch = cage::overlayVertices(all, cageXyz);
  REQUIRE(mismatch.nVertices == 8);
  REQUIRE(mismatch.status == 1);
  if (!ira::available()) {
    REQUIRE(shape.status == 1);
    REQUIRE(shape.pointGroup.empty());
  } else {
    REQUIRE(shape.status == 0);
    REQUIRE_FALSE(shape.pointGroup.empty());
    const auto same = cage::overlayVertices(cageXyz, cageXyz);
    REQUIRE(same.status == 0);
    REQUIRE(same.rmsd < 1e-6);
  }
}

static DumpCages contractibleCagesOnDump(const char *path, const char *spec) {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO(path, 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  const auto sig = cage::Signature::parse(spec);
  DumpCages out;
  out.rings = ring::cageRings(
      primitive::ringNetwork(nList, std::max(sig.maxRingSize(), 6)), yCloud,
      3.5);
  out.found = cage::findBySignature(out.rings, nList, sig);
  return out;
}

static std::size_t cagesOnCloud(
    molSys::PointCloud<molSys::Point<double>, double> cloud, const char *spec,
    double cutoff = 3.5) {
  auto nList = nneigh::neighListO(cutoff, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);
  const auto sig = cage::Signature::parse(spec);
  const auto rings = ring::cageRings(
      primitive::ringNetwork(nList, std::max(sig.maxRingSize(), 6)), cloud,
      cutoff);
  return cage::findBySignature(rings, nList, sig).size();
}

TEST_CASE("rings that wind the sI unit cell bound no cage", "[cage_enum]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, yCloud, 1);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  const auto all = primitive::ringNetwork(nList, 6);
  const auto kept = ring::contractibleRings(all, yCloud);
  REQUIRE(kept.size() < all.size());
  REQUIRE(contractibleCagesOnDump("traj/genice_sI.lammpstrj", "512")
              .found.size() == 2);
  REQUIRE(contractibleCagesOnDump("traj/genice_sI.lammpstrj", "51262")
              .found.size() == 6);
  REQUIRE(contractibleCagesOnDump("traj/genice_sI.lammpstrj", "51264")
              .found.empty());
  REQUIRE(contractibleCagesOnDump("traj/genice_sI.lammpstrj", "51268")
              .found.empty());
  REQUIRE(contractibleCagesOnDump("traj/genice_sI.lammpstrj", "sh")
              .found.empty());
}

TEST_CASE("a 2x2x2 sI supercell has eight times the unit-cell cages",
          "[cage_enum]") {
  molSys::PointCloud<molSys::Point<double>, double> unit;
  unit = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, unit, 1);
  REQUIRE(unit.nop == 46);
  REQUIRE(unit.box.size() >= 3);
  const double lx = unit.box[0];
  const double ly = unit.box[1];
  const double lz = unit.box[2];
  molSys::PointCloud<molSys::Point<double>, double> super;
  super.box = {2.0 * lx, 2.0 * ly, 2.0 * lz};
  super.boxLow = unit.boxLow;
  super.currentFrame = 1;
  int id = 1;
  for (int ix = 0; ix < 2; ++ix) {
    for (int iy = 0; iy < 2; ++iy) {
      for (int iz = 0; iz < 2; ++iz) {
        for (const auto &src : unit.pts) {
          auto p = src;
          p.x += ix * lx;
          p.y += iy * ly;
          p.z += iz * lz;
          p.atomID = id;
          super.idIndexMap[id] = super.nop;
          super.pts.push_back(p);
          ++super.nop;
          ++id;
        }
      }
    }
  }
  REQUIRE(super.nop == 8 * unit.nop);
  const char *specs[] = {"512", "51262", "51264", "51268", "sh"};
  const std::size_t unitCounts[] = {2, 6, 0, 0, 0};
  for (int i = 0; i < 5; ++i) {
    INFO(specs[i]);
    const auto one = cagesOnCloud(unit, specs[i]);
    const auto eight = cagesOnCloud(super, specs[i]);
    REQUIRE(one == unitCounts[i]);
    REQUIRE(eight == 8 * one);
  }
}

TEST_CASE("doubling the thin sH edge doubles each cage kind and adds none",
          "[cage_enum]") {
  molSys::PointCloud<molSys::Point<double>, double> cell;
  cell = sinp::readLammpsTrjO("traj/genice_sH.lammpstrj", 1, cell, 1);
  REQUIRE(cell.box.size() >= 6);
  const double lz = cell.box[2];
  molSys::PointCloud<molSys::Point<double>, double> doubled = cell;
  doubled.pts.clear();
  doubled.idIndexMap.clear();
  doubled.nop = 0;
  doubled.box[2] = 2.0 * lz;
  int id = 1;
  for (int iz = 0; iz < 2; ++iz) {
    for (const auto &src : cell.pts) {
      auto p = src;
      p.z += iz * lz;
      p.atomID = id;
      doubled.idIndexMap[id] = doubled.nop;
      doubled.pts.push_back(p);
      ++doubled.nop;
      ++id;
    }
  }
  const char *specs[] = {"512", "51262", "51264", "51268", "sh"};
  const std::size_t cellCounts[] = {12, 0, 0, 4, 8};
  for (int i = 0; i < 5; ++i) {
    INFO(specs[i]);
    const auto one = cagesOnCloud(cell, specs[i]);
    const auto two = cagesOnCloud(doubled, specs[i]);
    REQUIRE(one == cellCounts[i]);
    REQUIRE(two == 2 * one);
  }
}

TEST_CASE("a cell thinner than twice the cutoff is refused", "[cage_enum]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, yCloud, 1);
  const std::vector<std::vector<int>> none;
  yCloud.box = {6.0, 12.0, 12.0};
  REQUIRE_THROWS_AS(ring::cageRings(none, yCloud, 3.5), std::invalid_argument);
  yCloud.box = {12.0034886, 12.0034886, 12.0034886};
  REQUIRE_NOTHROW(ring::cageRings(none, yCloud, 3.5));
}

TEST_CASE("a cloud without a box keeps every ring", "[cage_enum]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  buildHCCloud(yCloud);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  const auto rings = primitive::ringNetwork(nList, 7);
  REQUIRE(ring::contractibleRings(rings, yCloud) == rings);
  yCloud.box.clear();
  REQUIRE(ring::contractibleRings(rings, yCloud) == rings);
}

TEST_CASE("the triclinic sH cell gives its three cage kinds and no others",
          "[cage_enum]") {
  const char *path = "traj/genice_sH.lammpstrj";
  REQUIRE(contractibleCagesOnDump(path, "512").found.size() == 12);
  REQUIRE(contractibleCagesOnDump(path, "51262").found.empty());
  REQUIRE(contractibleCagesOnDump(path, "51264").found.empty());
  REQUIRE(contractibleCagesOnDump(path, "51268").found.size() == 4);
  REQUIRE(contractibleCagesOnDump(path, "sh").found.size() == 8);
}

TEST_CASE("paired edges alone do not make a sphere", "[cage_enum]") {
  const std::vector<std::vector<int>> pillow = {{0, 1, 2}, {2, 1, 0}};
  REQUIRE(cage::isClosedPolyhedron(pillow, {0, 1}));
  REQUIRE_FALSE(cage::isSpherePolyhedron(pillow, {0, 1}));
  REQUIRE(cage::findBySignature(pillow, cage::Signature::parse("3:2")).empty());

  const std::vector<std::vector<int>> pinched = {
      {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2},
      {0, 4, 5}, {0, 5, 6}, {0, 6, 4}, {4, 6, 5}};
  const std::vector<int> both = {0, 1, 2, 3, 4, 5, 6, 7};
  REQUIRE(cage::isClosedPolyhedron(pinched, both));
  REQUIRE_FALSE(cage::isSpherePolyhedron(pinched, both));
  REQUIRE(cage::isSpherePolyhedron(pinched, {0, 1, 2, 3}));
  REQUIRE(cage::isSpherePolyhedron(pinched, {4, 5, 6, 7}));
  REQUIRE(cage::findBySignature(pinched, cage::Signature::parse("3:8")).empty());
  REQUIRE(cage::findBySignature(pinched, cage::Signature::parse("3:4")).size() ==
          2);

  const auto prism = hexPrismRings();
  REQUIRE(cage::isSpherePolyhedron(prism, {0, 1, 2, 3, 4, 5, 6, 7}));
}

TEST_CASE("liquid mW at the start of a nucleation run has no hydrate cages",
          "[cage_enum]") {
  for (const char *spec : {"512", "51262", "51264", "51268", "sh"}) {
    INFO("signature " << spec);
    const auto run = cagesOnDump("traj/nucleation_f1.lammpstrj", spec);
    REQUIRE(run.found.empty());
  }
}

TEST_CASE("a seven-vertex torus has paired edges but is no cage",
          "[cage_enum]") {
  std::vector<std::vector<int>> torus;
  for (int i = 0; i < 7; ++i) {
    torus.push_back({i, (i + 1) % 7, (i + 3) % 7});
    torus.push_back({i, (i + 3) % 7, (i + 2) % 7});
  }
  std::vector<int> all(torus.size());
  for (size_t i = 0; i < torus.size(); ++i) {
    all[i] = static_cast<int>(i);
  }
  REQUIRE(cage::isClosedPolyhedron(torus, all));
  REQUIRE(cage::findBySignature(torus, cage::Signature::parse("3:14")).empty());
  REQUIRE_FALSE(cage::isSpherePolyhedron(torus, all));
}

TEST_CASE("signatures searched together give the cages of each alone",
          "[cage_enum]") {
  for (const char *path :
       {"traj/genice_sH.lammpstrj", "traj/genice_sII.lammpstrj",
        "traj/nucleation_f1.lammpstrj"}) {
    INFO(path);
    molSys::PointCloud<molSys::Point<double>, double> yCloud;
    yCloud = sinp::readLammpsTrjO(path, 1, yCloud, 1);
    auto nList = nneigh::neighListO(3.5, yCloud, 1);
    nList = nneigh::neighbourListByIndex(yCloud, nList);
    const auto rings = ring::contractibleRings(
        primitive::ringNetwork(nList, 6), yCloud);
    std::vector<cage::Signature> sigs;
    for (const char *spec :
         {"512", "51262", "51264", "51268", "sh", "hc", "512", "4:6"}) {
      sigs.push_back(cage::Signature::parse(spec));
    }
    const auto together = cage::findBySignatures(rings, nList, sigs, {});
    REQUIRE(together.size() == sigs.size());
    for (size_t i = 0; i < sigs.size(); ++i) {
      INFO("signature " << sigs[i].str());
      const auto alone = cage::findBySignature(rings, nList, sigs[i], {});
      REQUIRE(together[i].size() == alone.size());
      for (size_t c = 0; c < alone.size(); ++c) {
        REQUIRE(together[i][c].faces == alone[c].faces);
        REQUIRE(together[i][c].vertices == alone[c].vertices);
        REQUIRE(together[i][c].certificate == alone[c].certificate);
        REQUIRE(together[i][c].signature == alone[c].signature);
      }
    }
  }
}

TEST_CASE("cage edge window parses none and a range", "[cage_enum]") {
  const auto fallback = cage::parseEdgeWindow("");
  REQUIRE(fallback.enabled);
  REQUIRE(fallback.lo == 2.65);
  REQUIRE(fallback.hi == 3.15);
  REQUIRE(cage::edgeWindowSpec(fallback) == "2.65,3.15");
  REQUIRE_FALSE(cage::parseEdgeWindow("none").enabled);
  REQUIRE_FALSE(cage::parseEdgeWindow("OFF").enabled);
  REQUIRE_FALSE(cage::parseEdgeWindow("false").enabled);
  REQUIRE(cage::edgeWindowSpec(cage::parseEdgeWindow("none")) == "none");
  const auto colon = cage::parseEdgeWindow("2.40:3.40");
  REQUIRE(colon.enabled);
  REQUIRE(colon.lo == 2.4);
  REQUIRE(colon.hi == 3.4);
  const auto slash = cage::parseEdgeWindow(" 0.265/0.315 ");
  REQUIRE(std::abs(slash.lo - 0.265) < 1e-12);
  REQUIRE(std::abs(slash.hi - 0.315) < 1e-12);
  REQUIRE_THROWS_AS(cage::parseEdgeWindow("2.65-3.15"), std::invalid_argument);
  REQUIRE_THROWS_AS(cage::parseEdgeWindow("3,2"), std::invalid_argument);
  REQUIRE_THROWS_AS(cage::parseEdgeWindow("-1,3"), std::invalid_argument);
}

namespace {

using CageCloud = molSys::PointCloud<molSys::Point<double>, double>;

CageCloud replicateCloud(const CageCloud &unit, int nx, int ny, int nz) {
  REQUIRE(unit.box.size() >= 3);
  const double lx = unit.box[0];
  const double ly = unit.box[1];
  const double lz = unit.box[2];
  CageCloud super;
  super.box = {static_cast<double>(nx) * lx, static_cast<double>(ny) * ly,
               static_cast<double>(nz) * lz};
  super.boxLow = unit.boxLow;
  super.currentFrame = 1;
  int id = 1;
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        for (const auto &src : unit.pts) {
          auto p = src;
          p.x += ix * lx;
          p.y += iy * ly;
          p.z += iz * lz;
          p.atomID = id;
          super.idIndexMap[id] = super.nop;
          super.pts.push_back(p);
          ++super.nop;
          ++id;
        }
      }
    }
  }
  return super;
}

std::vector<std::vector<std::vector<int>>>
cagesForWindow(const CageCloud &cloud,
               const std::vector<std::vector<int>> &rings,
               const std::vector<std::vector<int>> &nList,
               const cage::EdgeWindow &window) {
  std::vector<std::array<double, 3>> xyz(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
  }
  const std::vector<cage::Signature> sigs = {
      cage::Signature::parse("512"),
      cage::Signature::parse("51262"),
      cage::Signature::parse("51264"),
  };
  const auto found = cage::findBySignatures(rings, nList, sigs, {}, &xyz,
                                            &cloud.box, window, &cloud.boxLow);
  std::vector<std::vector<std::vector<int>>> vertices(found.size());
  for (std::size_t i = 0; i < found.size(); ++i) {
    vertices[i].reserve(found[i].size());
    for (const auto &cage : found[i]) {
      vertices[i].push_back(cage.vertices);
    }
  }
  return vertices;
}

void expectWindowAgrees(const CageCloud &cloud, std::size_t n512,
                        std::size_t n51262, std::size_t n51264) {
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);
  const auto rings =
      ring::cageRings(primitive::ringNetwork(nList, 6), cloud, 3.5);
  const auto on = cagesForWindow(cloud, rings, nList, {});
  const auto off = cagesForWindow(cloud, rings, nList, cage::parseEdgeWindow("none"));
  REQUIRE(on.size() == 3);
  REQUIRE(off.size() == 3);
  const std::size_t want[] = {n512, n51262, n51264};
  for (std::size_t i = 0; i < 3; ++i) {
    REQUIRE(on[i].size() == want[i]);
    REQUIRE(off[i].size() == want[i]);
    REQUIRE(on[i] == off[i]);
  }
}

} // namespace

TEST_CASE("default and disabled edge windows agree on the reference clathrates",
          "[cage_enum]") {
  CageCloud si;
  si = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, si, 1);
  REQUIRE(si.nop == 46);
  expectWindowAgrees(si, 2, 6, 0);
  {
    auto nList = nneigh::neighListO(3.5, si, 1);
    nList = nneigh::neighbourListByIndex(si, nList);
    const auto rings =
        ring::cageRings(primitive::ringNetwork(nList, 6), si, 3.5);
    const auto tight =
        cagesForWindow(si, rings, nList, cage::parseEdgeWindow("1,2"));
    REQUIRE(tight[0].empty());
    REQUIRE(tight[1].empty());
    REQUIRE(tight[2].empty());
  }

  CageCloud sii;
  sii = sinp::readLammpsTrjO("traj/genice_sII.lammpstrj", 1, sii, 1);
  REQUIRE(sii.nop == 136);
  expectWindowAgrees(sii, 16, 0, 8);

  const auto si222 = replicateCloud(si, 2, 2, 2);
  REQUIRE(si222.nop == 8 * si.nop);
  expectWindowAgrees(si222, 16, 48, 0);

  const auto sii10 = replicateCloud(sii, 10, 10, 10);
  REQUIRE(sii10.nop == 1000 * sii.nop);
  expectWindowAgrees(sii10, 16000, 0, 8000);
}

TEST_CASE("default and disabled edge windows agree on triclinic sH",
          "[cage_enum]") {
  // 136 waters: four formula units of structure H. Each unit is
  // 3 * 5^12, 2 * 4^3 5^6 6^3 and 1 * 5^12 6^8.
  // Wrapping the bound spans as orthogonal lengths drops the cages
  // whose edges cross the xy shear.
  CageCloud cell;
  cell = sinp::readLammpsTrjO("traj/genice_sH.lammpstrj", 1, cell, 1);
  REQUIRE(cell.nop == 136);
  REQUIRE(cell.box.size() >= 6);
  auto nList = nneigh::neighListO(3.5, cell, 1);
  nList = nneigh::neighbourListByIndex(cell, nList);
  const auto rings =
      ring::cageRings(primitive::ringNetwork(nList, 6), cell, 3.5);
  std::vector<std::array<double, 3>> xyz(static_cast<std::size_t>(cell.nop));
  for (int i = 0; i < cell.nop; ++i) {
    const auto &p = cell.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
  }
  const std::vector<cage::Signature> sigs = {
      cage::Signature::parse("512"),
      cage::Signature::parse("51262"),
      cage::Signature::parse("51264"),
      cage::Signature::parse("51268"),
      cage::Signature::parse("sh"),
  };
  const std::size_t want[] = {12, 0, 0, 4, 8};
  const auto on = cage::findBySignatures(rings, nList, sigs, {}, &xyz,
                                         &cell.box, {});
  const auto off = cage::findBySignatures(rings, nList, sigs, {}, &xyz,
                                          &cell.box, cage::parseEdgeWindow("none"));
  REQUIRE(on.size() == 5);
  REQUIRE(off.size() == 5);
  for (std::size_t i = 0; i < 5; ++i) {
    INFO(sigs[i].str());
    REQUIRE(on[i].size() == want[i]);
    REQUIRE(off[i].size() == want[i]);
  }
}

TEST_CASE("tilted sH edge window keeps every hexagonal hydrate cage",
          "[cage_enum]") {
  CageCloud sh;
  sh = sinp::readLammpsTrjO("traj/genice_sH.lammpstrj", 1, sh, 1);
  REQUIRE(sh.nop == 136);
  REQUIRE(sh.box.size() >= 6);
  // Two 34-water sH cells: 6*2 dodecahedra, no 5^12 6^2, no 5^12 6^4.
  // The tilt used to make the edge window wrap the padded bounds, which
  // dropped the cages whose contacts cross that face.
  expectWindowAgrees(sh, 12, 0, 0);

  auto nList = nneigh::neighListO(3.5, sh, 1);
  nList = nneigh::neighbourListByIndex(sh, nList);
  const auto rings = ring::cageRings(primitive::ringNetwork(nList, 8), sh, 3.5);
  std::vector<std::array<double, 3>> xyz(static_cast<std::size_t>(sh.nop));
  for (int i = 0; i < sh.nop; ++i) {
    const auto &p = sh.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
  }
  const std::vector<cage::Signature> sigs = {
      cage::Signature::parse("sh"),
      cage::Signature::parse("51268"),
  };
  const auto on = cage::findBySignatures(rings, nList, sigs, {}, &xyz, &sh.box,
                                         {}, &sh.boxLow);
  const auto off = cage::findBySignatures(rings, nList, sigs, {}, &xyz, &sh.box,
                                          cage::parseEdgeWindow("none"),
                                          &sh.boxLow);
  REQUIRE(on[0].size() == 8);
  REQUIRE(on[1].size() == 4);
  REQUIRE(off[0].size() == on[0].size());
  REQUIRE(off[1].size() == on[1].size());
}
