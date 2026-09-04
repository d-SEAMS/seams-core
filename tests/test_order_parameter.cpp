#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <utility>
#include <vector>

// Helper to build a PointCloud from a list of (x,y,z) coordinates
static molSys::PointCloud<molSys::Point<double>, double>
makeCloud(const std::vector<std::array<double, 3>> &coords,
          double boxLen = 100.0) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {boxLen, boxLen, boxLen};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  for (int i = 0; i < static_cast<int>(coords.size()); i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i;
    p.molID = 0;
    p.x = coords[i][0];
    p.y = coords[i][1];
    p.z = coords[i][2];
    cloud.pts.push_back(p);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = static_cast<int>(coords.size());
  return cloud;
}

// Regression test for projAreaSingleRing return order.
// The function should return {areaXY, areaXZ, areaYZ} (index 0=XY, 1=XZ, 2=YZ).
// The caller calcCoverageArea reads [0]=XY, [1]=XZ, [2]=YZ.
TEST_CASE("projAreaSingleRing returns areas in XY, XZ, YZ order",
          "[order_parameter]") {
  // A 2x3 rectangle in the XZ plane (y=5 for all points).
  // Vertices: (0,5,0), (2,5,0), (2,5,3), (0,5,3)
  // XY projected area: all y are the same, so area = 0
  // XZ projected area: 2 * 3 = 6
  // YZ projected area: all y are the same, so area = 0
  auto cloud = makeCloud({{0, 5, 0}, {2, 5, 0}, {2, 5, 3}, {0, 5, 3}});
  std::vector<int> ring = {0, 1, 2, 3};

  auto areas = topoparam::projAreaSingleRing(cloud, ring);

  REQUIRE(areas.size() == 3);
  // [0] = XY area = 0 (flat in XY since all y identical)
  REQUIRE_THAT(areas[0], Catch::Matchers::WithinAbs(0.0, 1e-10));
  // [1] = XZ area = 6.0
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(6.0, 1e-10));
  // [2] = YZ area = 0 (flat in YZ since all y identical)
  REQUIRE_THAT(areas[2], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

// Additional test: rectangle in XY plane
TEST_CASE("projAreaSingleRing XY plane rectangle", "[order_parameter]") {
  // A 4x3 rectangle in the XY plane (z=0 for all).
  // Vertices: (0,0,0), (4,0,0), (4,3,0), (0,3,0)
  // XY area: 4*3 = 12
  // XZ area: 0 (all z identical)
  // YZ area: 0 (all z identical)
  auto cloud = makeCloud({{0, 0, 0}, {4, 0, 0}, {4, 3, 0}, {0, 3, 0}});
  std::vector<int> ring = {0, 1, 2, 3};

  auto areas = topoparam::projAreaSingleRing(cloud, ring);

  REQUIRE_THAT(areas[0], Catch::Matchers::WithinAbs(12.0, 1e-10));
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(0.0, 1e-10));
  REQUIRE_THAT(areas[2], Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("rodgerF4 is cos 3 phi on a known H-O-O-H dihedral",
          "[order_parameter]") {
  // Two waters. Outer hydrogens give a 90 degree H-O-O-H dihedral
  // (phi = pi/2, cos 3phi = 0) or a planar 0 degree one (cos 3phi = 1).
  auto twoWaters = [](double h2y, double h2z) {
    // O1 (0,0,0) mol 1; H_outer (0,1,0); H_inner (1,0,0)
    // O2 (3,0,0) mol 2; H_outer (3,h2y,h2z); H_inner (2,0,0)
    auto cloud = makeCloud({{0, 0, 0}, {0, 1, 0}, {1, 0, 0},
                            {3, 0, 0}, {3, h2y, h2z}, {2, 0, 0}});
    cloud.pts[0].type = 1;
    cloud.pts[0].molID = 1;
    cloud.pts[0].atomID = 1;
    cloud.pts[1].type = 2;
    cloud.pts[1].molID = 1;
    cloud.pts[1].atomID = 2;
    cloud.pts[2].type = 2;
    cloud.pts[2].molID = 1;
    cloud.pts[2].atomID = 3;
    cloud.pts[3].type = 1;
    cloud.pts[3].molID = 2;
    cloud.pts[3].atomID = 4;
    cloud.pts[4].type = 2;
    cloud.pts[4].molID = 2;
    cloud.pts[4].atomID = 5;
    cloud.pts[5].type = 2;
    cloud.pts[5].molID = 2;
    cloud.pts[5].atomID = 6;
    cloud.idIndexMap.clear();
    for (int i = 0; i < 6; i++) {
      cloud.idIndexMap[cloud.pts[static_cast<std::size_t>(i)].atomID] = i;
    }
    return cloud;
  };
  // nList by atom ID, leading self
  const std::vector<std::vector<int>> nList = {
      {1, 4}, {2}, {3}, {4, 1}, {5}, {6}};
  auto planar = twoWaters(1.0, 0.0);
  const auto f0 = topoparam::rodgerF4(planar, nList, 1, 2);
  REQUIRE(std::isfinite(f0[0]));
  REQUIRE_THAT(f0[0], Catch::Matchers::WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(f0[3], Catch::Matchers::WithinAbs(1.0, 1e-9));
  REQUIRE_THAT(topoparam::meanFinite(f0), Catch::Matchers::WithinAbs(1.0, 1e-9));

  auto right = twoWaters(0.0, 1.0);
  const auto f90 = topoparam::rodgerF4(right, nList, 1, 2);
  REQUIRE_THAT(f90[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
  REQUIRE_THAT(f90[3], Catch::Matchers::WithinAbs(0.0, 1e-9));
}

TEST_CASE("rodgerF4 is NaN on mW with no hydrogens", "[order_parameter]") {
  auto cloud = makeCloud({{0, 0, 0}, {3, 0, 0}});
  cloud.pts[0].atomID = 1;
  cloud.pts[1].atomID = 2;
  cloud.idIndexMap = {{1, 0}, {2, 1}};
  const std::vector<std::vector<int>> nList = {{1, 2}, {2, 1}};
  const auto f4 = topoparam::rodgerF4(cloud, nList, 1, 2);
  REQUIRE_FALSE(std::isfinite(f4[0]));
  REQUIRE_FALSE(std::isfinite(f4[1]));
  REQUIRE_FALSE(std::isfinite(topoparam::meanFinite(f4)));
}

namespace {

double wrapLen(double x, double L) {
  x -= L * std::floor(x / L);
  if (x < 0.0) {
    x += L;
  }
  if (x >= L) {
    x = 0.0;
  }
  return x;
}

// Ice-rule orientation: each oxygen donates exactly two hydrogens.
bool assignIceRules(const std::vector<std::vector<int>> &adj,
                    std::vector<std::pair<int, int>> &owned) {
  const int nO = static_cast<int>(adj.size());
  std::vector<std::pair<int, int>> bonds;
  std::vector<std::vector<int>> inc(static_cast<std::size_t>(nO));
  for (int i = 0; i < nO; i++) {
    for (int j : adj[static_cast<std::size_t>(i)]) {
      if (i < j) {
        inc[static_cast<std::size_t>(i)].push_back(
            static_cast<int>(bonds.size()));
        inc[static_cast<std::size_t>(j)].push_back(
            static_cast<int>(bonds.size()));
        bonds.push_back({i, j});
      }
    }
  }
  if (bonds.empty()) {
    return false;
  }
  std::vector<int> owner(bonds.size(), -1);
  std::vector<int> outc(static_cast<std::size_t>(nO), 0);

  auto propagate = [&]() -> bool {
    bool changed = true;
    while (changed) {
      changed = false;
      for (int v = 0; v < nO; v++) {
        int freeB = 0;
        for (int b : inc[static_cast<std::size_t>(v)]) {
          if (owner[static_cast<std::size_t>(b)] < 0) {
            ++freeB;
          }
        }
        if (outc[static_cast<std::size_t>(v)] > 2) {
          return false;
        }
        if (outc[static_cast<std::size_t>(v)] + freeB < 2) {
          return false;
        }
        if (outc[static_cast<std::size_t>(v)] == 2) {
          for (int b : inc[static_cast<std::size_t>(v)]) {
            if (owner[static_cast<std::size_t>(b)] < 0) {
              const int oth = bonds[static_cast<std::size_t>(b)].first == v
                                  ? bonds[static_cast<std::size_t>(b)].second
                                  : bonds[static_cast<std::size_t>(b)].first;
              owner[static_cast<std::size_t>(b)] = oth;
              ++outc[static_cast<std::size_t>(oth)];
              changed = true;
              if (outc[static_cast<std::size_t>(oth)] > 2) {
                return false;
              }
            }
          }
        }
        if (outc[static_cast<std::size_t>(v)] + freeB == 2 && freeB > 0) {
          for (int b : inc[static_cast<std::size_t>(v)]) {
            if (owner[static_cast<std::size_t>(b)] < 0) {
              owner[static_cast<std::size_t>(b)] = v;
              ++outc[static_cast<std::size_t>(v)];
              changed = true;
            }
          }
        }
      }
    }
    return true;
  };

  std::function<bool()> dfs = [&]() -> bool {
    if (!propagate()) {
      return false;
    }
    int undecided = -1;
    for (std::size_t b = 0; b < bonds.size(); b++) {
      if (owner[b] < 0) {
        undecided = static_cast<int>(b);
        break;
      }
    }
    if (undecided < 0) {
      for (int c : outc) {
        if (c != 2) {
          return false;
        }
      }
      return true;
    }
    const auto snapOwner = owner;
    const auto snapOut = outc;
    const int i = bonds[static_cast<std::size_t>(undecided)].first;
    const int j = bonds[static_cast<std::size_t>(undecided)].second;
    for (int cand : {i, j}) {
      owner = snapOwner;
      outc = snapOut;
      owner[static_cast<std::size_t>(undecided)] = cand;
      ++outc[static_cast<std::size_t>(cand)];
      if (dfs()) {
        return true;
      }
    }
    owner = snapOwner;
    outc = snapOut;
    return false;
  };

  if (!dfs()) {
    return false;
  }
  owned.clear();
  owned.reserve(bonds.size());
  for (std::size_t b = 0; b < bonds.size(); b++) {
    const int o = owner[b];
    const int oth =
        bonds[b].first == o ? bonds[b].second : bonds[b].first;
    owned.push_back({o, oth});
  }
  return true;
}

void addIceHydrogens(molSys::PointCloud<molSys::Point<double>, double> &cloud,
                     const std::vector<std::pair<int, int>> &owned) {
  int nextId = 0;
  for (const auto &p : cloud.pts) {
    nextId = std::max(nextId, p.atomID);
  }
  ++nextId;
  for (const auto &bond : owned) {
    const int o = bond.first;
    const int oth = bond.second;
    const auto dr = gen::relDist(cloud, oth, o);
    const double r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
    REQUIRE(r2 > 0.0);
    const double inv = 1.0 / std::sqrt(r2);
    molSys::Point<double> h;
    h.type = 2;
    h.molID = cloud.pts[static_cast<std::size_t>(o)].molID;
    h.atomID = nextId++;
    h.x = wrapLen(cloud.pts[static_cast<std::size_t>(o)].x + dr[0] * inv,
                  cloud.box[0]);
    h.y = wrapLen(cloud.pts[static_cast<std::size_t>(o)].y + dr[1] * inv,
                  cloud.box[1]);
    h.z = wrapLen(cloud.pts[static_cast<std::size_t>(o)].z + dr[2] * inv,
                  cloud.box[2]);
    cloud.pts.push_back(h);
    cloud.idIndexMap[h.atomID] = static_cast<int>(cloud.pts.size()) - 1;
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
}

std::vector<std::vector<int>>
oxygenAdj(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
          const std::vector<std::vector<int>> &nList) {
  std::vector<std::vector<int>> adj(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; i++) {
    if (static_cast<std::size_t>(i) >= nList.size() ||
        nList[static_cast<std::size_t>(i)].size() < 2) {
      continue;
    }
    for (std::size_t k = 1; k < nList[static_cast<std::size_t>(i)].size();
         k++) {
      const auto it =
          cloud.idIndexMap.find(nList[static_cast<std::size_t>(i)][k]);
      if (it == cloud.idIndexMap.end()) {
        continue;
      }
      adj[static_cast<std::size_t>(i)].push_back(it->second);
    }
  }
  return adj;
}

molSys::PointCloud<molSys::Point<double>, double> iceIhOxygens() {
  const double a = 4.5115;
  const double c = 7.3463;
  const double b = a * std::sqrt(3.0);
  const int nx = 2;
  const int ny = 2;
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {nx * a, ny * b, c};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  const double frac[4][3] = {{1.0 / 3.0, 2.0 / 3.0, 1.0 / 16.0},
                             {2.0 / 3.0, 1.0 / 3.0, 9.0 / 16.0},
                             {2.0 / 3.0, 1.0 / 3.0, 15.0 / 16.0},
                             {1.0 / 3.0, 2.0 / 3.0, 7.0 / 16.0}};
  for (int ix = 0; ix < nx; ix++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int s = 0; s < 2; s++) {
        const double dx = (static_cast<double>(ix) + 0.5 * s) * a;
        const double dy = (static_cast<double>(iy) + 0.5 * s) * b;
        for (const auto &f : frac) {
          molSys::Point<double> p;
          p.type = 1;
          p.x = wrapLen(a * f[0] + (-0.5 * a) * f[1] + dx, cloud.box[0]);
          p.y = wrapLen((0.5 * a * std::sqrt(3.0)) * f[1] + dy, cloud.box[1]);
          p.z = wrapLen(c * f[2], cloud.box[2]);
          p.atomID = static_cast<int>(cloud.pts.size()) + 1;
          p.molID = p.atomID;
          cloud.pts.push_back(p);
          cloud.idIndexMap[p.atomID] = static_cast<int>(cloud.pts.size()) - 1;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  return cloud;
}

double f4OnOxygenCloud(
    molSys::PointCloud<molSys::Point<double>, double> oxy) {
  auto nO = nneigh::kNearestNeighbourList(oxy, 4, 3.5, 1, true);
  for (int i = 0; i < oxy.nop; i++) {
    REQUIRE(nO[static_cast<std::size_t>(i)].size() == 5);
  }
  std::vector<std::pair<int, int>> owned;
  REQUIRE(assignIceRules(oxygenAdj(oxy, nO), owned));
  addIceHydrogens(oxy, owned);
  auto nList = nneigh::kNearestNeighbourList(oxy, 4, 3.5, 1, true);
  const auto f4 = topoparam::rodgerF4(oxy, nList, 1, 2);
  const double mean = topoparam::meanFinite(f4);
  REQUIRE(std::isfinite(mean));
  return mean;
}

} // namespace

TEST_CASE("rodgerF4 is near -0.4 on ice Ih when hydrogens exist",
          "[order_parameter]") {
  const double mean = f4OnOxygenCloud(iceIhOxygens());
  UNSCOPED_INFO("F4 ice Ih mean=" << mean);
  REQUIRE(mean > -0.55);
  REQUIRE(mean < -0.25);
}

TEST_CASE("rodgerF4 is finite on exampleTraj when hydrogens are kept",
          "[order_parameter]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj("traj/exampleTraj.lammpstrj", 1, cloud);
  REQUIRE(cloud.nop > 0);
  auto nList = nneigh::kNearestNeighbourList(cloud, 4, 3.5, 2, true);
  const auto f4 = topoparam::rodgerF4(cloud, nList, 2, 1);
  int nFinite = 0;
  for (double v : f4) {
    nFinite += std::isfinite(v) ? 1 : 0;
  }
  REQUIRE(nFinite > 0);
  REQUIRE(std::isfinite(topoparam::meanFinite(f4)));
}

TEST_CASE("rodgerF4 is near 0.7 on filled sI when hydrogens exist",
          "[order_parameter]") {
  molSys::PointCloud<molSys::Point<double>, double> sI;
  sI = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, sI, 1);
  REQUIRE(sI.nop == 46);
  const double mean = f4OnOxygenCloud(sI);
  UNSCOPED_INFO("F4 sI mean=" << mean);
  // Perfect 0 K ice-rule sI sits above the thermal literature 0.7.
  REQUIRE(mean > 0.50);
  REQUIRE(mean < 0.95);
}

TEST_CASE("jumpRotorTau90 is finite on a 90-degree H-H rotation",
          "[order_parameter]") {
  auto a = makeCloud({{0, 0, 0}, {1, 0, 0}, {-1, 0, 0}}, 10.0);
  a.pts[0].type = 1;
  a.pts[1].type = 2;
  a.pts[2].type = 2;
  a.pts[0].molID = 1;
  a.pts[1].molID = 1;
  a.pts[2].molID = 1;
  auto b = a;
  b.pts[1].x = 0.0;
  b.pts[1].y = 1.0;
  b.pts[2].x = 0.0;
  b.pts[2].y = -1.0;
  const double tau = topoparam::jumpRotorTau90(a, b, 1.0, 1, 2);
  REQUIRE(std::isfinite(tau));
  REQUIRE_THAT(tau, Catch::Matchers::WithinAbs(1.0, 1e-12));
}

TEST_CASE("jumpRotorTau90 is NaN on a static ice Ih dump",
          "[order_parameter]") {
  auto ih = iceIhOxygens();
  auto nO = nneigh::kNearestNeighbourList(ih, 4, 3.5, 1, true);
  std::vector<std::pair<int, int>> owned;
  REQUIRE(assignIceRules(oxygenAdj(ih, nO), owned));
  addIceHydrogens(ih, owned);
  const double tau = topoparam::jumpRotorTau90(ih, ih, 1.0, 1, 2);
  REQUIRE_FALSE(std::isfinite(tau));
}

TEST_CASE("CHILL+ layerCubicity reproduces the literature I_sd string",
          "[order_parameter]") {
  // Four layers along z at 0, 3.7, 7.4, 11.1 in a 14.8 box: C H C H
  std::vector<std::array<double, 3>> coords;
  for (int k = 0; k < 4; k++) {
    coords.push_back({1.0, 1.0, k * 3.7});
    coords.push_back({2.0, 2.0, k * 3.7});
  }
  auto cloud = makeCloud(coords, 14.8);
  for (int i = 0; i < 8; i++) {
    cloud.pts[static_cast<std::size_t>(i)].iceType =
        (i / 2) % 2 == 0 ? molSys::atom_state_type::cubic
                         : molSys::atom_state_type::hexagonal;
  }
  const auto st = topoparam::layerCubicity(cloud, 2, 3.7);
  REQUIRE(st.sequence == "CHCH");
  REQUIRE_THAT(st.phiC, Catch::Matchers::WithinAbs(0.5, 1e-12));
  REQUIRE(st.phiC > 0.0);
  REQUIRE(st.phiC < 1.0);
}

TEST_CASE("TUM stacking uses ring planes and stays empty without basal rings",
          "[order_parameter]") {
  // Two disjoint six-rings stacked along z, marked basal. Their centroids
  // fall in two H layers. A five-ring (clathrate face) does not vote.
  auto cloud = makeCloud({{0, 0, 0}, {2, 0, 0}, {3, 1.5, 0},
                          {2, 3, 0}, {0, 3, 0}, {-1, 1.5, 0},
                          {0, 0, 7.4}, {2, 0, 7.4}, {3, 1.5, 7.4},
                          {2, 3, 7.4}, {0, 3, 7.4}, {-1, 1.5, 7.4},
                          {1, 1, 3.7}, {2, 1, 3.7}, {2.5, 2, 3.7},
                          {1.5, 2.8, 3.7}, {0.5, 2, 3.7}},
                         14.8);
  const std::vector<std::vector<int>> rings = {
      {0, 1, 2, 3, 4, 5}, {6, 7, 8, 9, 10, 11}, {12, 13, 14, 15, 16}};
  std::vector<bool> basal = {true, true, false};
  std::vector<bool> equatorial = {false, false, false};
  const auto tum = topoparam::tumLayerStack(cloud, rings, basal, equatorial, 2, 3.7);
  REQUIRE(tum.sequence.size() == 4);
  REQUIRE(tum.sequence[0] == 'H');
  REQUIRE(tum.sequence[2] == 'H');
  REQUIRE(tum.phiC == 0.0);
  // the five-ring is not a plane: no C letter
  REQUIRE(tum.sequence.find('C') == std::string::npos);

  // A clathrate 5-ring is not a stacking plane even if a caller flags it.
  std::vector<bool> flaggedFive = {true, true, true};
  const auto five = topoparam::tumLayerStack(cloud, rings, flaggedFive,
                                            equatorial, 2, 3.7);
  REQUIRE(five.sequence == tum.sequence);
  REQUIRE(five.phiC == 0.0);

  std::vector<bool> eq = {false, true, false};
  const auto mixed = topoparam::tumLayerStack(cloud, rings, basal, eq, 2, 3.7);
  REQUIRE(mixed.phiC > 0.0);
  REQUIRE(mixed.phiC < 1.0);
  REQUIRE(mixed.sequence.find('C') != std::string::npos);
  REQUIRE(mixed.sequence.find('H') != std::string::npos);
}

TEST_CASE("ions stay off the type-filtered oxygen neighbour graph",
          "[order_parameter]") {
  auto cloud = makeCloud({{0, 0, 0}, {2.8, 0, 0}, {1.4, 2.4, 0},
                          {0, 0, 2.8}, {10, 10, 10}},
                         20.0);
  cloud.pts[4].type = 3;
  cloud.pts[4].atomID = 99;
  cloud.idIndexMap.clear();
  for (int i = 0; i < cloud.nop; i++) {
    cloud.idIndexMap[cloud.pts[static_cast<std::size_t>(i)].atomID] = i;
  }
  const auto nList = nneigh::kNearestNeighbourList(cloud, 4, 5.5, 1, true);
  REQUIRE(nList.size() == static_cast<std::size_t>(cloud.nop));
  REQUIRE(nList[4].size() == 1);
  REQUIRE(nList[4][0] == 99);
  for (int i = 0; i < 4; i++) {
    REQUIRE(std::find(nList[static_cast<std::size_t>(i)].begin(),
                      nList[static_cast<std::size_t>(i)].end(),
                      99) == nList[static_cast<std::size_t>(i)].end());
  }
}

TEST_CASE("normHeightPercent uses recovered lz not tilt",
          "[order_parameter]") {
  auto cloud = makeCloud({{0.0, 0.0, 0.0}});
  cloud.box = {61.0, 12.0, 50.0, 60.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double h = topoparam::normHeightPercent(cloud, 10, 2.5);
  REQUIRE_THAT(h, Catch::Matchers::WithinAbs(50.0, 1e-10));
}
