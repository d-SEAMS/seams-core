#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <phase.hpp>
#include <seams_c_api.h>
#include <seams_input.hpp>

#include <cmath>
#include <vector>

static molSys::PointCloud<molSys::Point<double>, double>
cloudFrom(const std::vector<std::array<double, 3>> &xyz,
          const std::array<double, 3> &box, int type = 1) {
  molSys::PointCloud<molSys::Point<double>, double> c;
  c.box = {box[0], box[1], box[2]};
  c.boxLow = {0, 0, 0};
  c.currentFrame = 1;
  for (std::size_t i = 0; i < xyz.size(); i++) {
    molSys::Point<double> p;
    p.type = type;
    p.atomID = static_cast<int>(i) + 1;
    p.molID = p.atomID;
    p.x = xyz[i][0];
    p.y = xyz[i][1];
    p.z = xyz[i][2];
    c.pts.push_back(p);
    c.idIndexMap[p.atomID] = static_cast<int>(i);
  }
  c.nop = static_cast<int>(xyz.size());
  return c;
}

TEST_CASE("a hexagonal channel is not a closed 512 cage", "[phase]") {
  const double r = 2.75;
  std::vector<std::array<double, 3>> xyz;
  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 6; i++) {
      const double th = i * 3.14159265358979323846 / 3.0;
      xyz.push_back({5 + r * std::cos(th), 5 + r * std::sin(th), 2.0 + k * 2.75});
    }
  }
  auto cloud = cloudFrom(xyz, {20, 20, 20});
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  nList = nneigh::neighbourListByIndex(cloud, nList);
  const auto rings = primitive::ringNetwork(nList, 6);
  const auto closed = cage::findBySignature(rings, nList, cage::Signature::parse("512"));
  REQUIRE(closed.empty());
  REQUIRE(phase::openChannelCount(cloud, nList) > 0);
}

TEST_CASE("Ih and a flipped-proton copy have different proton keys", "[phase]") {
  auto cloud = cloudFrom({{0, 0, 0}, {1, 0, 0}, {-0.3, 0.9, 0}}, {10, 10, 10});
  cloud.pts[0].type = 1;
  cloud.pts[1].type = 2;
  cloud.pts[2].type = 2;
  cloud.pts[0].molID = 1;
  cloud.pts[1].molID = 1;
  cloud.pts[2].molID = 1;
  auto flipped = cloud;
  flipped.pts[1].x = -1.0;
  flipped.pts[2].x = 0.3;
  flipped.pts[2].y = -0.9;
  const auto k0 = phase::protonKey(cloud, 1, 2);
  const auto k1 = phase::protonKey(flipped, 1, 2);
  REQUIRE(k0 != k1);
}

TEST_CASE("mobile hydrogens report a finite MSD", "[phase]") {
  auto a = cloudFrom({{0, 0, 0}, {1, 0, 0}}, {10, 10, 10});
  a.pts[0].type = 1;
  a.pts[1].type = 2;
  auto b = a;
  b.pts[1].x = 2.0;
  const double msd = phase::hydrogenMSD(a, b, 2);
  REQUIRE(std::isfinite(msd));
  REQUIRE_THAT(msd, Catch::Matchers::WithinAbs(1.0, 1e-12));
}

TEST_CASE("ice XXI library hits a 152-site BCT cell and misses ice Ih",
          "[phase]") {
  std::vector<std::array<double, 3>> xyz;
  const double a = 20.197;
  const double c = 7.891;
  const int nx = 8;
  const int ny = 8;
  const int nz = 3;
  xyz.reserve(152);
  for (int i = 0; i < nx && static_cast<int>(xyz.size()) < 152; i++) {
    for (int j = 0; j < ny && static_cast<int>(xyz.size()) < 152; j++) {
      for (int k = 0; k < nz && static_cast<int>(xyz.size()) < 152; k++) {
        xyz.push_back({(i + 0.5) * a / nx, (j + 0.5) * a / ny, (k + 0.5) * c / nz});
      }
    }
  }
  REQUIRE(xyz.size() == 152);
  auto xxi = cloudFrom(xyz, {a, a, c});
  const auto hit = phase::iceXXILibrary(xxi);
  REQUIRE(hit.match);
  molSys::PointCloud<molSys::Point<double>, double> ih;
  ih = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, ih, 1);
  REQUIRE_FALSE(phase::iceXXILibrary(ih).match);
}

TEST_CASE("dense null bins as HDA/MDA; ice I does not", "[phase]") {
  molSys::PointCloud<molSys::Point<double>, double> sI;
  sI = sinp::readLammpsTrjO("traj/genice_sI.lammpstrj", 1, sI, 1);
  const auto rhoI = phase::localDensity(sI, 3.5);
  double meanI = 0.0;
  for (double r : rhoI) {
    meanI += r;
  }
  meanI /= static_cast<double>(rhoI.size());
  std::vector<std::array<double, 3>> packed;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      for (int k = 0; k < 5; k++) {
        packed.push_back({i * 2.3, j * 2.3, k * 2.3});
      }
    }
  }
  auto nullc = cloudFrom(packed, {11.5, 11.5, 11.5});
  const auto rhoN = phase::localDensity(nullc, 3.5);
  double meanN = 0.0;
  for (double r : rhoN) {
    meanN += r;
  }
  meanN /= static_cast<double>(rhoN.size());
  REQUIRE(meanN > meanI);
  const double cut = 0.5 * (meanI + meanN);
  REQUIRE(phase::glassFromDensity(meanI, cut, meanN - 1e-9) ==
          phase::GlassKind::ice);
  REQUIRE(phase::glassFromDensity(meanN, cut, meanN - 1e-9) ==
          phase::GlassKind::hda);
}

TEST_CASE("LAMMPS compute dump column is seams_chill_plus", "[phase][lammps]") {
  const double xyz[] = {0, 0, 0, 2.75, 0, 0, 1.375, 2.38, 0, 1.375, 0.79, 2.25};
  const double box[] = {0, 0, 0, 10, 10, 10};
  int labels[4] = {-1, -1, -1, -1};
  REQUIRE(seams_chill_plus(xyz, 4, box, labels) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(labels[i] >= 0);
  }
}
