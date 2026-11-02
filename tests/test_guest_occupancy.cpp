#include <catch2/catch_test_macros.hpp>

#include <cage.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <site.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

// Eight "water" vertices on a cube of edge 4 whose centre sits at (2, 2, 2)
// and a second cube shifted along x by the box length minus its edge, so
// that its vertices straddle the periodic boundary; guests are appended.
Cloud twoCubes(const std::vector<std::array<double, 3>> &guests) {
  Cloud cloud;
  const double L = 20.0;
  int id = 0;
  auto add = [&](double x, double y, double z, int type) {
    molSys::Point<double> p;
    p.x = std::fmod(x + L, L);
    p.y = std::fmod(y + L, L);
    p.z = std::fmod(z + L, L);
    p.type = type;
    p.atomID = id + 1;
    p.molID = id + 1;
    cloud.pts.push_back(p);
    cloud.idIndexMap[id + 1] = id;
    ++id;
  };
  for (int c = 0; c < 2; c++) {
    const double ox = c == 0 ? 0.0 : L - 2.0;  // second cube spans x in [18, 22) -> wraps
    for (int i = 0; i < 8; i++) {
      add(ox + 4.0 * (i & 1), 4.0 * ((i >> 1) & 1), 4.0 * ((i >> 2) & 1), 1);
    }
  }
  for (const auto &g : guests) {
    add(g[0], g[1], g[2], 2);
  }
  cloud.nop = id;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

std::vector<int> range(int lo, int hi) {
  std::vector<int> v;
  for (int i = lo; i < hi; i++) {
    v.push_back(i);
  }
  return v;
}

// Faces of the first cube of twoCubes, vertex i at 4 * (i & 1, i >> 1 & 1, i >> 2 & 1).
const std::vector<std::vector<int>> kCubeRings = {{0, 1, 3, 2}, {4, 5, 7, 6},
                                                   {0, 1, 5, 4}, {2, 3, 7, 6},
                                                   {0, 2, 6, 4}, {1, 3, 7, 5}};
const std::vector<int> kCubeFaces = {0, 1, 2, 3, 4, 5};

// A hexahedron whose top face is a skew quadrilateral: corners 4 and 6 at
// z = 4, corners 5 and 7 at z = 6. The four sides are planar. `topFirst`
// names the corner the top ring is listed from.
struct SkewBox {
  Cloud cloud;
  std::vector<std::vector<int>> rings;
};

SkewBox skewBox(int topFirst) {
  SkewBox box;
  const std::vector<std::array<double, 3>> xyz = {
      {0, 0, 0}, {4, 0, 0}, {4, 4, 0}, {0, 4, 0},
      {0, 0, 4}, {4, 0, 6}, {4, 4, 4}, {0, 4, 6}};
  for (int i = 0; i < 8; i++) {
    molSys::Point<double> p;
    p.x = xyz[static_cast<std::size_t>(i)][0] + 5.0;
    p.y = xyz[static_cast<std::size_t>(i)][1] + 5.0;
    p.z = xyz[static_cast<std::size_t>(i)][2] + 5.0;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    box.cloud.pts.push_back(p);
    box.cloud.idIndexMap[i + 1] = i;
  }
  box.cloud.nop = 8;
  box.cloud.box = {20.0, 20.0, 20.0};
  box.cloud.boxLow = {0.0, 0.0, 0.0};
  std::vector<int> top = {4, 5, 6, 7};
  std::rotate(top.begin(), top.begin() + (topFirst - 4), top.end());
  box.rings = {{0, 1, 2, 3}, top,          {0, 1, 5, 4},
               {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
  return box;
}

} // namespace

TEST_CASE("the periodic centroid of a cage across the boundary is inside it", "[guests]") {
  const auto cloud = twoCubes({});
  const auto c0 = site::periodicCentroid(cloud, range(0, 8));
  REQUIRE(std::fabs(c0[0] - 2.0) < 1e-9);
  REQUIRE(std::fabs(c0[1] - 2.0) < 1e-9);
  REQUIRE(std::fabs(c0[2] - 2.0) < 1e-9);
  const auto c1 = site::periodicCentroid(cloud, range(8, 16));
  // the wrapped cube spans x in [18, 22): centre at 20 == 0 modulo the box
  const double wrapped = std::fmod(c1[0] + 20.0, 20.0);
  REQUIRE((std::fabs(wrapped) < 1e-9 || std::fabs(wrapped - 20.0) < 1e-9));
  REQUIRE(std::fabs(c1[1] - 2.0) < 1e-9);
}

TEST_CASE("guests are assigned to the cage whose centre they are nearest", "[guests]") {
  // one guest at the first centre, one near the wrapped centre seen from
  // the far side of the box, one in open space
  const auto cloud = twoCubes({{2.0, 2.0, 2.0}, {19.5, 2.0, 2.2}, {10.0, 10.0, 10.0}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const auto occ = site::guestOccupancy(cloud, cages, {16, 17, 18}, 3.0);
  REQUIRE(occ.cageOfGuest == std::vector<int>{0, 1, -1});
  REQUIRE(occ.guestsPerCage == std::vector<int>{1, 1});
  REQUIRE(occ.occupied == 2);
  REQUIRE(occ.multiply == 0);
  REQUIRE(occ.free == 1);
  REQUIRE(std::fabs(occ.centreDistance[0]) < 1e-9);
  REQUIRE(std::fabs(occ.centreDistance[1] - std::sqrt(0.25 + 0.04)) < 1e-9);
  REQUIRE(occ.centreDistance[2] == -1.0);
}

TEST_CASE("ray-parity occupancy keeps a centre guest and rejects one outside",
          "[guests]") {
  const auto cloud = twoCubes({{2.0, 2.0, 2.0}, {10.0, 10.0, 10.0}, {2.0, 2.0, -0.5}});
  // Cube 0: vertices 0-7. Faces as 6 rings of 4.
  const std::vector<std::vector<int>> rings = {
      {0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4},
      {2, 3, 7, 6}, {0, 2, 6, 4}, {1, 3, 7, 5},
      {8, 9, 11, 10}, {12, 13, 15, 14}, {8, 9, 13, 12},
      {10, 11, 15, 14}, {8, 10, 14, 12}, {9, 11, 15, 13}};
  const std::vector<std::vector<int>> faces = {
      {0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}};
  const auto occ = site::guestOccupancyInside(cloud, rings, faces, {16, 17, 18});
  REQUIRE(occ.cageOfGuest == std::vector<int>{0, -1, -1});
  REQUIRE(occ.guestsPerCage == std::vector<int>{1, 0});
  REQUIRE(occ.occupied == 1);
  REQUIRE(occ.free == 2);
}

TEST_CASE("a guest on a face of the cage is inside", "[guests][inside]") {
  // on the top face, on the edge the top face shares with the x = 4 face,
  // at a corner, and just above the top face
  const auto cloud =
      twoCubes({{2.3, 1.7, 4.0}, {4.0, 1.3, 4.0}, {4.0, 4.0, 4.0}, {2.3, 1.7, 4.001}});
  for (int g = 16; g < 19; g++) {
    const auto &p = cloud.pts[static_cast<std::size_t>(g)];
    const auto tr = site::pointInFacesTrace(cloud, kCubeRings, kCubeFaces, {p.x, p.y, p.z});
    REQUIRE(tr.onFace);
    REQUIRE(tr.inside);
    REQUIRE(tr.rays == 0);
  }
  const auto &above = cloud.pts[19];
  const auto tr = site::pointInFacesTrace(cloud, kCubeRings, kCubeFaces,
                                          {above.x, above.y, above.z});
  REQUIRE_FALSE(tr.onFace);
  REQUIRE_FALSE(tr.inside);
  const auto occ = site::guestOccupancyInside(cloud, kCubeRings, {kCubeFaces}, {16, 17, 18, 19});
  REQUIRE(occ.cageOfGuest == std::vector<int>{0, 0, 0, -1});
}

TEST_CASE("a skew face is fanned from its centroid, not from its first vertex",
          "[guests][inside]") {
  // Under the top face's centroid (2, 2, 5) but above the diagonal 4-6
  // (z = 4): inside the centroid fan, outside a fan from corner 4.
  const std::array<double, 3> under = {7.0, 7.1, 9.5};
  // Above the centroid fan (z = 5 there) but below the diagonal 5-7 (z = 6):
  // outside the centroid fan, inside a fan from corner 5.
  const std::array<double, 3> over = {7.0, 6.9, 10.5};
  for (int first : {4, 5, 6, 7}) {
    const auto box = skewBox(first);
    REQUIRE(site::pointInFaces(box.cloud, box.rings, kCubeFaces, under));
    REQUIRE_FALSE(site::pointInFaces(box.cloud, box.rings, kCubeFaces, over));
  }
  auto box = skewBox(4);
  for (const auto &g : {under, over}) {
    molSys::Point<double> p;
    p.x = g[0];
    p.y = g[1];
    p.z = g[2];
    p.type = 2;
    p.atomID = box.cloud.nop + 1;
    p.molID = box.cloud.nop + 1;
    box.cloud.pts.push_back(p);
    box.cloud.idIndexMap[p.atomID] = box.cloud.nop;
    ++box.cloud.nop;
  }
  const auto occ = site::guestOccupancyInside(box.cloud, box.rings, {kCubeFaces}, {8, 9});
  REQUIRE(occ.cageOfGuest == std::vector<int>{0, -1});
}

TEST_CASE("a ray through a triangle edge or vertex is retried along the next direction",
          "[guests][inside]") {
  const auto cloud = twoCubes({});
  const auto &dirs = site::insideRayDirections();
  REQUIRE(dirs.size() >= 2);
  const auto d0 = dirs[0];
  REQUIRE(d0[0] > 0.0);
  REQUIRE(d0[1] > 0.0);
  REQUIRE(d0[2] > 0.0);
  auto behind = [&](const std::array<double, 3> &target, double s) {
    return std::array<double, 3>{target[0] - s * d0[0], target[1] - s * d0[1],
                                 target[2] - s * d0[2]};
  };
  // the first ray runs through a cube corner, the centroid of the top face
  // (the fan vertex its four triangles share) and the middle of a cube edge
  for (const auto &target : {std::array<double, 3>{4.0, 4.0, 4.0},
                             std::array<double, 3>{2.0, 2.0, 4.0},
                             std::array<double, 3>{4.0, 4.0, 2.0}}) {
    const auto tr = site::pointInFacesTrace(cloud, kCubeRings, kCubeFaces, behind(target, 1.0));
    REQUIRE_FALSE(tr.onFace);
    REQUIRE(tr.rays >= 2);
    REQUIRE(tr.inside);
    REQUIRE(tr.crossings % 2 == 1);
  }
  // from below the cube, the first ray enters through the fan vertex of the
  // bottom face and leaves again
  const auto out = site::pointInFacesTrace(cloud, kCubeRings, kCubeFaces,
                                           behind({2.0, 2.0, 0.0}, 1.0));
  REQUIRE(out.rays >= 2);
  REQUIRE_FALSE(out.inside);
  REQUIRE(out.crossings % 2 == 0);
  // a clean ray from the centre is decided by the first direction
  const auto centre = site::pointInFacesTrace(cloud, kCubeRings, kCubeFaces, {2.0, 2.0, 2.0});
  REQUIRE(centre.rays == 1);
  REQUIRE(centre.crossings == 1);
  REQUIRE(centre.inside);
}

TEST_CASE("two guests in one cage count as multiple occupancy", "[guests]") {
  const auto cloud = twoCubes({{1.5, 2.0, 2.0}, {2.5, 2.0, 2.0}});
  const std::vector<std::vector<int>> cages = {range(0, 8), range(8, 16)};
  const auto occ = site::guestOccupancy(cloud, cages, {16, 17}, 3.0);
  REQUIRE(occ.guestsPerCage == std::vector<int>{2, 0});
  REQUIRE(occ.occupied == 1);
  REQUIRE(occ.multiply == 1);
  REQUIRE(occ.free == 0);
  REQUIRE(occ.occupancyHistogram.size() >= 3);
  REQUIRE(occ.occupancyHistogram[0] == 1);
  REQUIRE(occ.occupancyHistogram[2] == 1);
  // a radius too small leaves every guest free
  const auto none = site::guestOccupancy(cloud, cages, {16, 17}, 0.1);
  REQUIRE(none.free == 2);
  REQUIRE(none.occupied == 0);
}

TEST_CASE("two H2 in a GenIce sII 51264 cage report integer occupancy",
          "[guests]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO("traj/genice_sII.lammpstrj", 1, yCloud, 1);
  REQUIRE(yCloud.nop > 0);
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  nList = nneigh::neighbourListByIndex(yCloud, nList);
  const auto sig = cage::Signature::parse("51264");
  const auto rings = primitive::ringNetwork(nList, std::max(sig.maxRingSize(), 6));
  const auto found = cage::findBySignature(rings, nList, sig);
  REQUIRE_FALSE(found.empty());
  REQUIRE(found[0].vertices.size() == 28);

  const auto centre = site::periodicCentroid(yCloud, found[0].vertices);
  auto addGuest = [&](double dx, double dy, double dz) {
    molSys::Point<double> p;
    p.x = centre[0] + dx;
    p.y = centre[1] + dy;
    p.z = centre[2] + dz;
    p.type = 2;
    p.atomID = yCloud.nop + 1;
    p.molID = yCloud.nop + 1;
    yCloud.pts.push_back(p);
    yCloud.idIndexMap[p.atomID] = yCloud.nop;
    ++yCloud.nop;
  };
  addGuest(-0.4, 0.0, 0.0);
  addGuest(0.4, 0.0, 0.0);

  std::vector<std::vector<int>> cages;
  cages.reserve(found.size());
  for (const auto &c : found) {
    cages.push_back(c.vertices);
  }
  const auto occ =
      site::guestOccupancy(yCloud, cages, {yCloud.nop - 2, yCloud.nop - 1}, 4.0);
  REQUIRE(occ.occupied == 1);
  REQUIRE(occ.multiply == 1);
  REQUIRE(occ.free == 0);
  REQUIRE(occ.guestsPerCage[0] == 2);
  REQUIRE(occ.occupancyHistogram.size() >= 3);
  REQUIRE(occ.occupancyHistogram[2] == 1);
  const int histSum =
      std::accumulate(occ.occupancyHistogram.begin(), occ.occupancyHistogram.end(), 0);
  REQUIRE(histSum == static_cast<int>(found.size()));
}
