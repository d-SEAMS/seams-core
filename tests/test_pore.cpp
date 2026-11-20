#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <neighbours.hpp>
#include <numbers>
#include <pore.hpp>
#include <stdexcept>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
constexpr double kPi = std::numbers::pi_v<double>;

void addAtom(Cloud &cloud, int type, double x, double y, double z) {
  molSys::Point<double> pt;
  pt.type = type;
  pt.atomID = static_cast<int>(cloud.pts.size()) + 1;
  pt.molID = pt.atomID;
  pt.x = x;
  pt.y = y;
  pt.z = z;
  cloud.idIndexMap[pt.atomID] = static_cast<int>(cloud.pts.size());
  cloud.pts.push_back(pt);
  cloud.nop = static_cast<int>(cloud.pts.size());
}

Cloud ortho(double lx, double ly, double lz) {
  Cloud cloud;
  cloud.box = {lx, ly, lz};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

Cloud triclinic(double lx, double ly, double lz, double xy, double xz, double yz) {
  Cloud cloud;
  REQUIRE(nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, 0.0, 0.0, 0.0, lx, ly, lz, xy,
                                        xz, yz));
  return cloud;
}

void placeFrac(Cloud &cloud, int type, double lx, double ly, double lz, double xy, double xz,
               double yz, double sx, double sy, double sz) {
  const double x = lx * sx + xy * sy + xz * sz;
  const double y = ly * sy + yz * sz;
  const double z = lz * sz;
  addAtom(cloud, type, x, y, z);
}

double integrate(const pore::Profile &profile) {
  double acc = 0.0;
  for (std::size_t i = 0; i < profile.volume.size(); ++i) {
    acc += profile.number_density[i] * profile.volume[i];
  }
  return acc;
}

void requireFlat(const pore::Profile &profile, double density) {
  for (std::size_t i = 0; i < profile.number_density.size(); ++i) {
    REQUIRE_THAT(profile.number_density[i], Catch::Matchers::WithinRel(density, 1e-9));
    REQUIRE_THAT(profile.mass_density[i], Catch::Matchers::WithinRel(density, 1e-9));
    REQUIRE(profile.volume[i] > 0.0);
  }
}

} // namespace

TEST_CASE("uniform gas is flat along every fractional axis", "[pore]") {
  const int n = 8;
  const double lx = 10.0;
  const double ly = 8.0;
  const double lz = 12.0;
  const double xy = 2.0;
  const double xz = 1.0;
  const double yz = -1.5;
  const double volume = lx * ly * lz;
  for (int axis = 0; axis < 3; ++axis) {
    auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
    for (int i = 0; i < n; ++i) {
      double s[3] = {0.5, 0.5, 0.5};
      s[axis] = (static_cast<double>(i) + 0.5) / static_cast<double>(n);
      placeFrac(cloud, 1, lx, ly, lz, xy, xz, yz, s[0], s[1], s[2]);
    }
    const auto profile = pore::fractionalAxis(cloud, 0, n, axis, 0.0, nullptr);
    REQUIRE(profile.coordinate.size() == static_cast<std::size_t>(n));
    requireFlat(profile, static_cast<double>(n) / volume);
    REQUIRE_THAT(integrate(profile), Catch::Matchers::WithinAbs(static_cast<double>(n), 1e-9));
    REQUIRE_THAT(profile.coordinate.front(),
                 Catch::Matchers::WithinAbs(0.5 / static_cast<double>(n), 1e-12));
  }
}

TEST_CASE("orthorhombic fractional profile matches the known slab volume", "[pore]") {
  auto cloud = ortho(4.0, 5.0, 10.0);
  for (int i = 0; i < 5; ++i) {
    addAtom(cloud, 1, 2.0, 2.5, (static_cast<double>(i) + 0.5) * 2.0);
  }
  const auto profile = pore::fractionalAxis(cloud, 0, 5, 2, 0.0, nullptr);
  requireFlat(profile, 5.0 / 200.0);
  for (double volume : profile.volume) {
    REQUIRE_THAT(volume, Catch::Matchers::WithinAbs(40.0, 1e-12));
  }
}

TEST_CASE("origin shift moves a boundary atom into the predicted bin", "[pore]") {
  auto cloud = ortho(10.0, 10.0, 10.0);
  addAtom(cloud, 1, 1.0, 1.0, 0.0);
  const auto unshifted = pore::fractionalAxis(cloud, 0, 2, 2, 0.0, nullptr);
  const auto shifted = pore::fractionalAxis(cloud, 0, 2, 2, 0.5, nullptr);
  REQUIRE(unshifted.number_density[0] > 0.0);
  REQUIRE_THAT(unshifted.number_density[1], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE_THAT(shifted.number_density[0], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE(shifted.number_density[1] > 0.0);
  const auto wrapped = pore::fractionalAxis(cloud, 0, 2, 2, 1.5, nullptr);
  REQUIRE_THAT(wrapped.number_density[1],
               Catch::Matchers::WithinRel(shifted.number_density[1], 1e-12));
}

TEST_CASE("a lattice translation leaves the fractional profile unchanged", "[pore]") {
  const double lx = 9.0;
  const double ly = 7.0;
  const double lz = 11.0;
  const double xy = 1.5;
  const double xz = -0.4;
  const double yz = 0.8;
  auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
  for (int i = 0; i < 6; ++i) {
    const double sz = (static_cast<double>(i) + 0.5) / 6.0;
    placeFrac(cloud, 1, lx, ly, lz, xy, xz, yz, 0.3, 0.7, sz);
  }
  const auto before = pore::fractionalAxis(cloud, 0, 6, 2, 0.25, nullptr);
  for (auto &pt : cloud.pts) {
    pt.x += lx + xy + xz;
    pt.y += ly + yz;
    pt.z += lz;
  }
  const auto after = pore::fractionalAxis(cloud, 0, 6, 2, 0.25, nullptr);
  for (std::size_t i = 0; i < before.number_density.size(); ++i) {
    REQUIRE_THAT(after.number_density[i],
                 Catch::Matchers::WithinRel(before.number_density[i], 1e-9));
  }
}

TEST_CASE("an atom outside the primary cell wraps into the fractional bin", "[pore]") {
  auto cloud = ortho(8.0, 8.0, 8.0);
  addAtom(cloud, 1, 8.0 + 2.0, 1.0, 1.0);
  const auto profile = pore::fractionalAxis(cloud, 0, 4, 0, 0.0, nullptr);
  REQUIRE(profile.number_density[1] > 0.0);
  REQUIRE_THAT(profile.number_density[0], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE_THAT(profile.number_density[2], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE_THAT(integrate(profile), Catch::Matchers::WithinAbs(1.0, 1e-12));
}

TEST_CASE("type filter and masses set number and mass density", "[pore]") {
  auto cloud = ortho(2.0, 2.0, 2.0);
  addAtom(cloud, 1, 1.0, 1.0, 0.5);
  addAtom(cloud, 2, 1.0, 1.0, 1.5);
  const double masses[2] = {3.0, 5.0};
  const auto kept = pore::fractionalAxis(cloud, 1, 2, 2, 0.0, masses);
  REQUIRE_THAT(integrate(kept), Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(kept.number_density[0], Catch::Matchers::WithinRel(1.0 / 4.0, 1e-12));
  REQUIRE_THAT(kept.mass_density[0], Catch::Matchers::WithinRel(3.0 / 4.0, 1e-12));
  REQUIRE_THAT(kept.number_density[1], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE_THROWS_AS(pore::fractionalAxis(cloud, 0, 0, 2, 0.0, nullptr), std::out_of_range);
  REQUIRE_THROWS_AS(pore::fractionalAxis(cloud, 0, 2, 3, 0.0, nullptr), std::out_of_range);
  const double bad[2] = {-1.0, 1.0};
  REQUIRE_THROWS_AS(pore::fractionalAxis(cloud, 0, 2, 2, 0.0, bad), std::invalid_argument);
}

TEST_CASE("uniform gas is flat along an arbitrary direction", "[pore]") {
  const double side = 6.0;
  auto cloud = ortho(side, side, side);
  addAtom(cloud, 1, 0.2 * side, 0.2 * side, 0.5 * side);
  addAtom(cloud, 1, 0.8 * side, 0.8 * side, 0.4 * side);
  const double direction[3] = {1.0, 1.0, 0.0};
  const auto profile = pore::alongDirection(cloud, 0, 2, direction, 0.0, nullptr);
  const double density = 2.0 / (side * side * side);
  requireFlat(profile, density);
  REQUIRE_THAT(profile.volume[0], Catch::Matchers::WithinRel(0.5 * side * side * side, 1e-9));
  REQUIRE_THAT(profile.volume[1], Catch::Matchers::WithinRel(0.5 * side * side * side, 1e-9));
  const double zdir[3] = {0.0, 0.0, 1.0};
  auto tilted = triclinic(10.0, 8.0, 12.0, 2.0, 1.0, -1.5);
  const int n = 8;
  for (int i = 0; i < n; ++i) {
    const double sz = (static_cast<double>(i) + 0.5) / static_cast<double>(n);
    placeFrac(tilted, 1, 10.0, 8.0, 12.0, 2.0, 1.0, -1.5, 0.4, 0.6, sz);
  }
  const auto along = pore::alongDirection(tilted, 0, n, zdir, 0.0, nullptr);
  const auto frac = pore::fractionalAxis(tilted, 0, n, 2, 0.0, nullptr);
  for (int i = 0; i < n; ++i) {
    REQUIRE_THAT(along.number_density[static_cast<std::size_t>(i)],
                 Catch::Matchers::WithinRel(frac.number_density[static_cast<std::size_t>(i)], 1e-9));
    REQUIRE_THAT(along.volume[static_cast<std::size_t>(i)],
                 Catch::Matchers::WithinRel(frac.volume[static_cast<std::size_t>(i)], 1e-9));
  }
  const double zero[3] = {0.0, 0.0, 0.0};
  REQUIRE_THROWS_AS(pore::alongDirection(cloud, 0, 4, zero, 0.0, nullptr), std::invalid_argument);
}

TEST_CASE("atoms on a cylinder shell have the analytic annular density", "[pore]") {
  const double lx = 10.0;
  const double ly = 10.0;
  const double lz = 20.0;
  auto cloud = ortho(lx, ly, lz);
  const double radius = 2.5;
  const int n = 12;
  for (int i = 0; i < n; ++i) {
    const double theta = 2.0 * kPi * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
    const double z = (static_cast<double>(i) + 0.5) / static_cast<double>(n) * lz;
    addAtom(cloud, 1, 5.0 + radius * std::cos(theta), 5.0 + radius * std::sin(theta), z);
  }
  const double centre[3] = {0.5, 0.5, 0.0};
  const auto profiles =
      pore::poreProfiles(cloud, 0, 6, 4, n, 6.0, 2, nullptr, centre, 0, 0.0, nullptr);
  REQUIRE_THAT(integrate(profiles.radial), Catch::Matchers::WithinAbs(static_cast<double>(n), 1e-8));
  for (int i = 0; i < 6; ++i) {
    if (i == 2) {
      continue;
    }
    REQUIRE_THAT(profiles.radial.number_density[static_cast<std::size_t>(i)],
                 Catch::Matchers::WithinAbs(0.0, 0.0));
  }
  const double shell = kPi * (9.0 - 4.0) * lz;
  REQUIRE_THAT(profiles.radial.volume[2], Catch::Matchers::WithinRel(shell, 1e-9));
  REQUIRE_THAT(profiles.radial.number_density[2],
               Catch::Matchers::WithinRel(static_cast<double>(n) / shell, 1e-9));
  const double sector = (lx * ly * lz) / 4.0;
  for (double volume : profiles.azimuthal.volume) {
    REQUIRE_THAT(volume, Catch::Matchers::WithinRel(sector, 1e-8));
  }
  requireFlat(profiles.azimuthal, static_cast<double>(n) / (lx * ly * lz));
  requireFlat(profiles.axial, static_cast<double>(n) / (lx * ly * lz));
}

TEST_CASE("triclinic shell volume is the Euclidean annulus times the edge", "[pore]") {
  const double lx = 12.0;
  const double ly = 10.0;
  const double lz = 16.0;
  const double xy = 3.0;
  const double xz = 1.0;
  const double yz = -2.0;
  auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
  const double clen = std::hypot(std::hypot(xz, yz), lz);
  const double radius = 1.5;
  const double centre[3] = {0.5, 0.5, 0.25};
  const int n = 8;
  // Place atoms by walking a Cartesian circle perpendicular to c, then
  // the profile must recover that radius.
  const double axis[3] = {xz / clen, yz / clen, lz / clen};
  double e1[3] = {1.0, 0.0, 0.0};
  const double proj = e1[0] * axis[0] + e1[1] * axis[1] + e1[2] * axis[2];
  e1[0] -= proj * axis[0];
  e1[1] -= proj * axis[1];
  e1[2] -= proj * axis[2];
  const double n1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  e1[0] /= n1;
  e1[1] /= n1;
  e1[2] /= n1;
  const double e2[3] = {axis[1] * e1[2] - axis[2] * e1[1], axis[2] * e1[0] - axis[0] * e1[2],
                        axis[0] * e1[1] - axis[1] * e1[0]};
  const double cx = lx * 0.5 + xy * 0.5 + xz * 0.25;
  const double cy = ly * 0.5 + yz * 0.25;
  const double cz = lz * 0.25;
  for (int i = 0; i < n; ++i) {
    const double theta = 2.0 * kPi * (static_cast<double>(i) + 0.5) / static_cast<double>(n);
    addAtom(cloud, 1, cx + radius * (std::cos(theta) * e1[0] + std::sin(theta) * e2[0]),
            cy + radius * (std::cos(theta) * e1[1] + std::sin(theta) * e2[1]),
            cz + radius * (std::cos(theta) * e1[2] + std::sin(theta) * e2[2]));
  }
  const auto profiles =
      pore::poreProfiles(cloud, 0, 4, 4, 1, 4.0, 2, nullptr, centre, 0, 0.0, nullptr);
  const double shell = kPi * (4.0 - 1.0) * clen;
  REQUIRE_THAT(profiles.radial.volume[1], Catch::Matchers::WithinRel(shell, 1e-8));
  REQUIRE_THAT(profiles.radial.number_density[1],
               Catch::Matchers::WithinRel(static_cast<double>(n) / shell, 1e-8));
  REQUIRE_THAT(profiles.radial.number_density[0], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE_THAT(integrate(profiles.radial), Catch::Matchers::WithinAbs(static_cast<double>(n), 1e-8));
  REQUIRE_THAT(profiles.axial.volume[0], Catch::Matchers::WithinRel(lx * ly * lz, 1e-9));
}

TEST_CASE("auto axis is the longest edge and auto centre is the wall centroid", "[pore]") {
  auto long_z = ortho(5.0, 6.0, 20.0);
  addAtom(long_z, 8, 1.0, 1.0, 1.0);
  const auto axis = pore::poreFrame(long_z, 8, -1, nullptr, nullptr);
  REQUIRE(axis.axis_index == 2);

  auto skewed = triclinic(3.0, 1.0, 4.0, 10.0, 0.0, 0.0);
  addAtom(skewed, 8, 1.0, 0.2, 1.0);
  const auto along_b = pore::poreFrame(skewed, 8, -1, nullptr, nullptr);
  REQUIRE(along_b.axis_index == 1);

  auto walls = ortho(10.0, 10.0, 10.0);
  addAtom(walls, 4, 7.0, 5.0, 1.0);
  addAtom(walls, 4, 3.0, 5.0, 2.0);
  addAtom(walls, 4, 5.0, 7.0, 3.0);
  addAtom(walls, 4, 5.0, 3.0, 4.0);
  const auto centre = pore::poreFrame(walls, 4, 2, nullptr, nullptr);
  REQUIRE_THAT(centre.centre_frac[0], Catch::Matchers::WithinAbs(0.5, 1e-8));
  REQUIRE_THAT(centre.centre_frac[1], Catch::Matchers::WithinAbs(0.5, 1e-8));

  auto corner = ortho(10.0, 10.0, 10.0);
  addAtom(corner, 4, 0.5, 5.0, 1.0);
  addAtom(corner, 4, 9.5, 5.0, 2.0);
  addAtom(corner, 4, 0.0, 6.5, 3.0);
  addAtom(corner, 4, 0.0, 3.5, 4.0);
  const auto wrapped = pore::poreFrame(corner, 4, 2, nullptr, nullptr);
  const double sx = wrapped.centre_frac[0] < 0.5 ? wrapped.centre_frac[0]
                                                 : wrapped.centre_frac[0] - 1.0;
  REQUIRE_THAT(sx, Catch::Matchers::WithinAbs(0.0, 1e-6));
  REQUIRE_THAT(wrapped.centre_frac[1], Catch::Matchers::WithinAbs(0.5, 1e-6));

  Cloud empty;
  empty.box = {1.0, 1.0, 1.0};
  empty.boxLow = {0.0, 0.0, 0.0};
  REQUIRE_THROWS_AS(pore::poreFrame(empty, 0, 2, nullptr, nullptr), std::invalid_argument);
  const double diagonal[3] = {1.0, 1.0, 0.0};
  REQUIRE_THROWS_AS(pore::poreProfiles(walls, 0, 2, 2, 2, 0.0, -1, diagonal, nullptr, 4, 0.0,
                                       nullptr),
                    std::invalid_argument);
}

TEST_CASE("host guest nearest image and contacts follow the triclinic cell", "[pore]") {
  const double lx = 10.0;
  const double ly = 8.0;
  const double lz = 9.0;
  const double xy = 2.0;
  const double xz = 0.5;
  const double yz = -1.0;
  auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
  placeFrac(cloud, 1, lx, ly, lz, xy, xz, yz, 0.5, 0.1, 0.5);
  placeFrac(cloud, 3, lx, ly, lz, xy, xz, yz, 0.5, 0.9, 0.5);
  const double dsy = 0.2;
  const double expected = std::hypot(xy * dsy, ly * dsy);
  const auto report = pore::hostGuest(cloud, 1, 0, expected);
  REQUIRE(report.guest.size() == 1);
  REQUIRE(report.nearest[0] == 1);
  REQUIRE(report.nearest_type[0] == 3);
  REQUIRE_THAT(report.distance[0], Catch::Matchers::WithinRel(expected, 1e-9));
  REQUIRE(report.contacts.size() == 1);
  REQUIRE(report.contacts[0] == 1);
  const auto outside = pore::hostGuest(cloud, 1, 3, expected * 0.5);
  REQUIRE(outside.nearest[0] == 1);
  REQUIRE(outside.contacts[0] == 0);

  auto counted = ortho(20.0, 20.0, 20.0);
  addAtom(counted, 1, 5.0, 5.0, 5.0);
  addAtom(counted, 3, 6.0, 5.0, 5.0);
  addAtom(counted, 3, 5.0, 8.0, 5.0);
  addAtom(counted, 4, 5.0, 5.0, 7.0);
  const auto by_type = pore::hostGuest(counted, 1, 0, 2.5);
  REQUIRE(by_type.nearest[0] == 1);
  REQUIRE(by_type.nearest_type[0] == 3);
  REQUIRE_THAT(by_type.distance[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE(by_type.framework_types.size() == 2);
  REQUIRE(by_type.framework_types[0] == 3);
  REQUIRE(by_type.framework_types[1] == 4);
  REQUIRE(by_type.contacts[0] == 1);
  REQUIRE(by_type.contacts[1] == 1);

  auto ortho_cloud = ortho(10.0, 10.0, 10.0);
  addAtom(ortho_cloud, 1, 0.2, 5.0, 5.0);
  addAtom(ortho_cloud, 7, 9.9, 5.0, 5.0);
  const auto image = pore::hostGuest(ortho_cloud, 1, 7, 1.0);
  REQUIRE_THAT(image.distance[0], Catch::Matchers::WithinAbs(0.3, 1e-12));
  REQUIRE(image.contacts[0] == 1);

  auto alone = ortho(4.0, 4.0, 4.0);
  addAtom(alone, 1, 1.0, 1.0, 1.0);
  const auto none = pore::hostGuest(alone, 1, 0, 1.0);
  REQUIRE(none.nearest[0] == -1);
  REQUIRE(none.distance[0] < 0.0);
  REQUIRE(none.contacts.empty());
  REQUIRE_THROWS_AS(pore::hostGuest(alone, 1, 0, -0.1), std::out_of_range);
  REQUIRE_THROWS_AS(pore::hostGuest(alone, 0, 1, 1.0), std::out_of_range);

  auto coincident = ortho(5.0, 5.0, 5.0);
  addAtom(coincident, 1, 1.0, 1.0, 1.0);
  addAtom(coincident, 2, 1.0, 1.0, 1.0);
  const auto touch = pore::hostGuest(coincident, 1, 2, 0.0);
  REQUIRE_THAT(touch.distance[0], Catch::Matchers::WithinAbs(0.0, 0.0));
  REQUIRE(touch.contacts[0] == 1);
}

TEST_CASE("a shifted guest keeps the same nearest framework site", "[pore]") {
  auto cloud = ortho(12.0, 9.0, 8.0);
  addAtom(cloud, 1, 1.0, 2.0, 3.0);
  addAtom(cloud, 5, 4.0, 2.0, 3.0);
  addAtom(cloud, 6, 1.0, 6.0, 3.0);
  const auto before = pore::hostGuest(cloud, 1, 0, 5.0);
  for (auto &pt : cloud.pts) {
    pt.x += 12.0;
    pt.y += 9.0;
    pt.z -= 8.0;
  }
  const auto after = pore::hostGuest(cloud, 1, 0, 5.0);
  REQUIRE(after.nearest[0] == before.nearest[0]);
  REQUIRE(after.nearest_type[0] == before.nearest_type[0]);
  REQUIRE_THAT(after.distance[0], Catch::Matchers::WithinAbs(before.distance[0], 1e-9));
  REQUIRE(after.contacts == before.contacts);
}
