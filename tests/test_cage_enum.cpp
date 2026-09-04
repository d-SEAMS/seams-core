#include <cage.hpp>
#include <cage_canon.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
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
