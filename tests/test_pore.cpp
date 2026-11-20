#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <neighbours.hpp>
#include <numbers>
#include <pore.hpp>
#include <stdexcept>
#include <vector>

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

double volumeSum(const pore::Profile &profile) {
  double sum = 0.0;
  for (double volume : profile.volume) {
    sum += volume;
  }
  return sum;
}

struct Face2 {
  double eu[2];
  double ev[2];
  double axial;
};

Face2 faceOf(double lx, double ly, double lz, double xy, double xz, double yz, int axis) {
  const double lattice[3][3] = {
      {lx, 0.0, 0.0},
      {xy, ly, 0.0},
      {xz, yz, lz},
  };
  const double *along = lattice[axis];
  const double axial2 = along[0] * along[0] + along[1] * along[1] + along[2] * along[2];
  const double axial = std::sqrt(axial2);
  const double eaxis[3] = {along[0] / axial, along[1] / axial, along[2] / axial};
  const int iu = (axis + 1) % 3;
  const int iv = (axis + 2) % 3;
  auto dot = [](const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };
  const double pu = dot(lattice[iu], along) / axial2;
  const double pv = dot(lattice[iv], along) / axial2;
  double raw_u[3];
  double raw_v[3];
  for (int k = 0; k < 3; ++k) {
    raw_u[k] = lattice[iu][k] - pu * along[k];
    raw_v[k] = lattice[iv][k] - pv * along[k];
  }
  const double nu = std::sqrt(dot(raw_u, raw_u));
  double e1[3];
  for (int k = 0; k < 3; ++k) {
    e1[k] = raw_u[k] / nu;
  }
  double e2[3] = {eaxis[1] * e1[2] - eaxis[2] * e1[1], eaxis[2] * e1[0] - eaxis[0] * e1[2],
                  eaxis[0] * e1[1] - eaxis[1] * e1[0]};
  const double n2 = std::sqrt(dot(e2, e2));
  for (int k = 0; k < 3; ++k) {
    e2[k] /= n2;
  }
  Face2 face;
  face.axial = axial;
  face.eu[0] = dot(raw_u, e1);
  face.eu[1] = dot(raw_u, e2);
  face.ev[0] = dot(raw_v, e1);
  face.ev[1] = dot(raw_v, e2);
  return face;
}

double nearestAxisRadius(const Face2 &face, double du, double dv) {
  double best = 1e300;
  for (int i = -3; i <= 3; ++i) {
    for (int j = -3; j <= 3; ++j) {
      const double x = (du + static_cast<double>(i)) * face.eu[0] +
                       (dv + static_cast<double>(j)) * face.ev[0];
      const double y = (du + static_cast<double>(i)) * face.eu[1] +
                       (dv + static_cast<double>(j)) * face.ev[1];
      best = std::min(best, std::hypot(x, y));
    }
  }
  return best;
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

TEST_CASE("shifted direction bins sum to the cell volume", "[pore]") {
  auto cube = ortho(10.0, 10.0, 10.0);
  addAtom(cube, 1, 5.0, 5.0, 5.0);
  const double zdir[3] = {0.0, 0.0, 1.0};
  for (double origin : {0.0, 0.3, 0.6, 0.9}) {
    const auto profile = pore::alongDirection(cube, 0, 4, zdir, origin, nullptr);
    for (double volume : profile.volume) {
      REQUIRE_THAT(volume, Catch::Matchers::WithinAbs(250.0, 1e-9));
    }
    REQUIRE_THAT(volumeSum(profile), Catch::Matchers::WithinAbs(1000.0, 1e-9));
  }

  const double directions[][3] = {
      {1.0, 1.0, 0.0}, {1.0, 2.0, 3.0}, {-1.0, 0.5, 2.0}, {1.0, -1.0, 1.0}, {1.0, 1e-7, 0.0},
  };
  const double origins[] = {0.0, 0.37, 0.61, 0.9};
  const double skew_lx = 10.0;
  const double skew_ly = 12.0;
  const double skew_lz = 14.0;
  const double skew_v = skew_lx * skew_ly * skew_lz;
  auto skew = triclinic(skew_lx, skew_ly, skew_lz, 4.5, -3.0, 5.5);
  addAtom(skew, 1, 1.0, 1.0, 1.0);
  auto tilted = triclinic(10.0, 8.0, 12.0, 2.0, 1.0, -1.5);
  addAtom(tilted, 1, 1.0, 1.0, 1.0);
  const double tilted_v = 10.0 * 8.0 * 12.0;
  for (const auto &direction : directions) {
    for (double origin : origins) {
      const auto on_cube = pore::alongDirection(cube, 0, 7, direction, origin, nullptr);
      const auto on_skew = pore::alongDirection(skew, 0, 7, direction, origin, nullptr);
      const auto on_tilt = pore::alongDirection(tilted, 0, 5, direction, origin, nullptr);
      REQUIRE_THAT(volumeSum(on_cube), Catch::Matchers::WithinRel(1000.0, 1e-9));
      REQUIRE_THAT(volumeSum(on_skew), Catch::Matchers::WithinRel(skew_v, 1e-9));
      REQUIRE_THAT(volumeSum(on_tilt), Catch::Matchers::WithinRel(tilted_v, 1e-9));
      for (double volume : on_cube.volume) {
        REQUIRE(volume > 0.0);
      }
    }
  }
}

TEST_CASE("a uniform fluid is flat for every pore binning", "[pore]") {
  const double side = 10.0;
  const double volume = side * side * side;
  const int bins = 4;
  const int n = 4;
  for (double origin : {0.0, 0.3, 0.6, 0.9}) {
    auto cloud = ortho(side, side, side);
    for (int ix = 0; ix < n; ++ix) {
      for (int iy = 0; iy < n; ++iy) {
        for (int iz = 0; iz < bins; ++iz) {
          const double sz = std::fmod(origin + (static_cast<double>(iz) + 0.5) / bins, 1.0);
          placeFrac(cloud, 1, side, side, side, 0.0, 0.0, 0.0,
                    (static_cast<double>(ix) + 0.5) / n, (static_cast<double>(iy) + 0.5) / n, sz);
        }
      }
    }
    const double direction[3] = {0.0, 0.0, 1.0};
    const double density = static_cast<double>(n * n * bins) / volume;
    const auto along = pore::alongDirection(cloud, 0, bins, direction, origin, nullptr);
    requireFlat(along, density);
    const auto axis = pore::fractionalAxis(cloud, 0, bins, 2, origin, nullptr);
    requireFlat(axis, density);
  }

  auto cube = ortho(side, side, side);
  const int grid = 16;
  for (int ix = 0; ix < grid; ++ix) {
    for (int iy = 0; iy < grid; ++iy) {
      for (int iz = 0; iz < grid; ++iz) {
        placeFrac(cube, 1, side, side, side, 0.0, 0.0, 0.0,
                  (static_cast<double>(ix) + 0.5) / grid, (static_cast<double>(iy) + 0.5) / grid,
                  (static_cast<double>(iz) + 0.5) / grid);
      }
    }
  }
  const double density = static_cast<double>(grid * grid * grid) / volume;
  const double diagonal[3] = {1.0, 2.0, 3.0};
  const auto diagonal_profile = pore::alongDirection(cube, 0, bins, diagonal, 0.61, nullptr);
  REQUIRE_THAT(volumeSum(diagonal_profile), Catch::Matchers::WithinRel(volume, 1e-9));
  // A cubic lattice projected on (1, 2, 3) is not a perfect Riemann sum of
  // the slice: the four densities sit about 2% off N/V and the gap does
  // not close when the grid is refined, because the direction is rational
  // in the lattice. A cumulative bin volume, the failure at a nonzero
  // origin, moves the density by tens of percent. The volume sum above is
  // the exact pin.
  for (std::size_t i = 0; i < diagonal_profile.number_density.size(); ++i) {
    REQUIRE(diagonal_profile.volume[i] > 0.0);
    REQUIRE_THAT(diagonal_profile.number_density[i], Catch::Matchers::WithinRel(density, 0.03));
    REQUIRE_THAT(diagonal_profile.mass_density[i],
                 Catch::Matchers::WithinRel(diagonal_profile.number_density[i], 1e-12));
  }

  const double centre[3] = {0.5, 0.5, 0.5};
  const auto two = pore::poreProfiles(cube, 0, 4, 2, bins, 0.0, 2, nullptr, centre, 0, 0.0, nullptr);
  REQUIRE_THAT(two.azimuthal.volume[0], Catch::Matchers::WithinRel(0.5 * volume, 1e-8));
  REQUIRE_THAT(two.azimuthal.volume[1], Catch::Matchers::WithinRel(0.5 * volume, 1e-8));
  requireFlat(two.azimuthal, density);
  requireFlat(two.axial, density);
  REQUIRE_THAT(volumeSum(two.radial), Catch::Matchers::WithinRel(volume, 1e-8));
  // On this 16^3 lattice the outer ring holds 640 points against 595
  // expected from the shell volume, 7.6% high, and the inner disks are
  // 4.5% low. The volumes of those disks are the analytic pi r^2 L.
  // A zero shell, the half-turn wedge failure, is many orders larger.
  for (std::size_t i = 0; i < two.radial.volume.size(); ++i) {
    if (two.radial.volume[i] > 1e-8 * volume) {
      REQUIRE_THAT(two.radial.number_density[i], Catch::Matchers::WithinRel(density, 0.08));
    }
  }

  const double lx = 10.0;
  const double ly = 12.0;
  const double lz = 14.0;
  const double xy = 4.5;
  const double xz = -3.0;
  const double yz = 5.5;
  const double skew_v = lx * ly * lz;
  auto skew = triclinic(lx, ly, lz, xy, xz, yz);
  for (int ix = 0; ix < grid; ++ix) {
    for (int iy = 0; iy < grid; ++iy) {
      for (int iz = 0; iz < grid; ++iz) {
        placeFrac(skew, 1, lx, ly, lz, xy, xz, yz, (static_cast<double>(ix) + 0.5) / grid,
                  (static_cast<double>(iy) + 0.5) / grid, (static_cast<double>(iz) + 0.5) / grid);
      }
    }
  }
  const double skew_density = static_cast<double>(grid * grid * grid) / skew_v;
  const auto skew_bins =
      pore::poreProfiles(skew, 0, 4, 2, bins, 0.0, 2, nullptr, centre, 0, 0.61, nullptr);
  REQUIRE_THAT(volumeSum(skew_bins.azimuthal), Catch::Matchers::WithinRel(skew_v, 1e-8));
  REQUIRE_THAT(skew_bins.azimuthal.volume[0], Catch::Matchers::WithinRel(0.5 * skew_v, 1e-6));
  REQUIRE_THAT(skew_bins.azimuthal.volume[1], Catch::Matchers::WithinRel(0.5 * skew_v, 1e-6));
  // The fractional grid is not aligned with the skewed azimuth or with
  // the shifted axial slabs, so the counts are not exactly equal. A
  // zero sector volume produces densities near 1e15.
  for (const pore::Profile *profile : {&skew_bins.azimuthal, &skew_bins.axial}) {
    REQUIRE_THAT(volumeSum(*profile), Catch::Matchers::WithinRel(skew_v, 1e-8));
    for (std::size_t i = 0; i < profile->volume.size(); ++i) {
      REQUIRE(profile->volume[i] > 0.0);
      REQUIRE_THAT(profile->number_density[i], Catch::Matchers::WithinRel(skew_density, 0.25));
      REQUIRE_THAT(profile->mass_density[i],
                   Catch::Matchers::WithinRel(profile->number_density[i], 1e-12));
    }
  }
  REQUIRE_THAT(volumeSum(skew_bins.radial), Catch::Matchers::WithinRel(skew_v, 1e-6));
  for (std::size_t i = 0; i < skew_bins.radial.volume.size(); ++i) {
    if (skew_bins.radial.volume[i] > 1e-6 * skew_v) {
      REQUIRE_THAT(skew_bins.radial.number_density[i], Catch::Matchers::WithinRel(skew_density, 0.25));
    }
  }
}

TEST_CASE("two azimuthal bins have the half-cell volume", "[pore]") {
  auto cube = ortho(10.0, 10.0, 10.0);
  addAtom(cube, 1, 5.0, 5.0, 5.0);
  const double centre[3] = {0.5, 0.5, 0.5};
  for (int bins : {1, 2, 3, 4}) {
    const auto profiles =
        pore::poreProfiles(cube, 0, 2, bins, 1, 0.0, 2, nullptr, centre, 0, 0.0, nullptr);
    REQUIRE_THAT(volumeSum(profiles.azimuthal), Catch::Matchers::WithinRel(1000.0, 1e-8));
    if (bins == 2) {
      REQUIRE_THAT(profiles.azimuthal.volume[0], Catch::Matchers::WithinRel(500.0, 1e-8));
      REQUIRE_THAT(profiles.azimuthal.volume[1], Catch::Matchers::WithinRel(500.0, 1e-8));
    }
  }
}

TEST_CASE("a corner radius is the nearest pore-axis image", "[pore]") {
  const double lx = 20.0;
  const double ly = 10.0 * std::sqrt(3.0);
  const double lz = 8.0;
  const double xy = 10.0;
  auto cloud = triclinic(lx, ly, lz, xy, 0.0, 0.0);
  const double centre[3] = {0.5, 0.5, 0.5};
  // Fractional (0, 0) is a corner of the centred hexagonal face. The
  // centred-face distance is half the long diagonal. The nearest image
  // of the axis is half the short diagonal.
  placeFrac(cloud, 1, lx, ly, lz, xy, 0.0, 0.0, 0.0, 0.0, 0.5);
  const Face2 face = faceOf(lx, ly, lz, xy, 0.0, 0.0, 2);
  const double centred = std::hypot(-0.5 * face.eu[0] - 0.5 * face.ev[0],
                                    -0.5 * face.eu[1] - 0.5 * face.ev[1]);
  const double nearest = nearestAxisRadius(face, -0.5, -0.5);
  REQUIRE(centred > nearest + 5.0);
  REQUIRE_THAT(nearest, Catch::Matchers::WithinAbs(10.0, 1e-8));
  const auto profiles =
      pore::poreProfiles(cloud, 0, 24, 4, 1, 24.0, 2, nullptr, centre, 0, 0.0, nullptr);
  int occupied = -1;
  for (int i = 0; i < 24; ++i) {
    if (profiles.radial.number_density[static_cast<std::size_t>(i)] > 0.0) {
      REQUIRE(occupied < 0);
      occupied = i;
    }
  }
  REQUIRE(occupied >= 0);
  const double dr = 1.0;
  const double lo = dr * static_cast<double>(occupied);
  const double hi = dr * static_cast<double>(occupied + 1);
  REQUIRE(nearest >= lo - 1e-8);
  REQUIRE(nearest < hi + 1e-8);
  REQUIRE(centred > hi);
  REQUIRE_THAT(volumeSum(profiles.radial), Catch::Matchers::WithinRel(lx * ly * lz, 1e-8));
}

TEST_CASE("a disk edge on a face vertex does not jump the shell volume", "[pore]") {
  const double lx = 20.0;
  const double ly = 10.0 * std::sqrt(3.0);
  const double lz = 30.0;
  const double xy = -10.0;
  const double xz = 6.0;
  const double yz = -8.0;
  auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
  addAtom(cloud, 1, 0.0, 0.0, 0.0);
  const Face2 face = faceOf(lx, ly, lz, xy, xz, yz, 2);
  const double centre[3] = {0.5, 0.5, 0.5};
  double vertex = 1e300;
  for (double s : {-1.0, 1.0}) {
    for (double t : {-1.0, 1.0}) {
      vertex = std::min(vertex, std::hypot(0.5 * (s * face.eu[0] + t * face.ev[0]),
                                           0.5 * (s * face.eu[1] + t * face.ev[1])));
    }
  }
  std::array<pore::Profile, 3> shells;
  int index = 0;
  for (double scale : {1.0 - 1e-9, 1.0, 1.0 + 1e-9}) {
    const auto profiles = pore::poreProfiles(cloud, 0, 2, 1, 1, 2.0 * vertex * scale, 2, nullptr,
                                             centre, 0, 0.0, nullptr);
    shells[static_cast<std::size_t>(index++)] = profiles.radial;
    REQUIRE_THAT(volumeSum(profiles.radial), Catch::Matchers::WithinRel(lx * ly * lz, 1e-6));
  }
  for (int bin = 0; bin < 2; ++bin) {
    const double left = shells[0].volume[static_cast<std::size_t>(bin)];
    const double mid = shells[1].volume[static_cast<std::size_t>(bin)];
    const double right = shells[2].volume[static_cast<std::size_t>(bin)];
    // A vertex that lands on the bin edge used to move more than a
    // thousand cubic angstroms between the two shells. The smooth change
    // across 1e-9 in the radius is far smaller.
    REQUIRE_THAT(mid, Catch::Matchers::WithinAbs(left, 1e-2));
    REQUIRE_THAT(mid, Catch::Matchers::WithinAbs(right, 1e-2));
  }

  // The radial shell clips the Voronoi cell, whose vertices share one
  // radius, so the public profiles above stay continuous even if a
  // circle through a single vertex is mis-counted. The centred face
  // does not: its short and long corners differ, and the gate's
  // counterexample is that polygon. A root that falls on the vertex
  // used to be dropped and the whole edge counted.
  std::vector<std::array<double, 2>> corners;
  const double sign[4][2] = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
  for (const auto &s : sign) {
    corners.push_back({0.5 * (s[0] * face.eu[0] + s[1] * face.ev[0]),
                       0.5 * (s[0] * face.eu[1] + s[1] * face.ev[1])});
  }
  const double turn = corners[0][0] * corners[1][1] - corners[0][1] * corners[1][0];
  if (turn < 0.0) {
    std::reverse(corners.begin(), corners.end());
  }
  const double cover = 2.0 * vertex + 1.0;
  std::array<double, 3> areas{};
  int area_index = 0;
  for (double scale : {1.0 - 1e-9, 1.0, 1.0 + 1e-9}) {
    areas[static_cast<std::size_t>(area_index++)] =
        pore::clippedDiskArea(corners, vertex * scale, cover, 0.0);
  }
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(areas[0], 1e-4));
  REQUIRE_THAT(areas[1], Catch::Matchers::WithinAbs(areas[2], 1e-4));
}

TEST_CASE("host guest contacts match a brute-force minimum image", "[pore]") {
  const double ox = 1.0;
  const double oy = -2.0;
  const double oz = 0.5;
  const double lx = 10.0;
  const double ly = 12.0;
  const double lz = 14.0;
  const double xy = 4.5;
  const double xz = -3.0;
  const double yz = 5.5;
  const double cutoff = 6.0;
  auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
  // The dump origin is not zero. Shift every Cartesian point.
  auto place = [&](int type, double sx, double sy, double sz) {
    const double x = ox + lx * sx + xy * sy + xz * sz;
    const double y = oy + ly * sy + yz * sz;
    const double z = oz + lz * sz;
    addAtom(cloud, type, x, y, z);
  };
  std::vector<std::array<double, 3>> guest_frac;
  std::vector<std::array<double, 3>> frame_frac;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 2; ++k) {
        const double s[3] = {(static_cast<double>(i) + 0.17) / 4.0,
                             (static_cast<double>(j) + 0.41) / 4.0,
                             (static_cast<double>(k) + 0.63) / 2.0};
        place(1, s[0], s[1], s[2]);
        guest_frac.push_back({s[0], s[1], s[2]});
      }
    }
  }
  const int nguest = static_cast<int>(guest_frac.size());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        const double s[3] = {(static_cast<double>(i) + 0.5) / 3.0,
                             (static_cast<double>(j) + 0.2) / 3.0,
                             (static_cast<double>(k) + 0.8) / 3.0};
        place(2, s[0], s[1], s[2]);
        frame_frac.push_back({s[0], s[1], s[2]});
      }
    }
  }
  int wrap_disagrees = 0;
  std::vector<double> expect_distance(static_cast<std::size_t>(nguest), 0.0);
  std::vector<int> expect_index(static_cast<std::size_t>(nguest), -1);
  std::vector<int> expect_contacts(static_cast<std::size_t>(nguest), 0);
  for (int g = 0; g < nguest; ++g) {
    const auto gf = guest_frac[static_cast<std::size_t>(g)];
    double best = 1e300;
    int best_index = -1;
    int contacts = 0;
    for (int f = 0; f < static_cast<int>(frame_frac.size()); ++f) {
      const auto ff = frame_frac[static_cast<std::size_t>(f)];
      double ds[3] = {gf[0] - ff[0], gf[1] - ff[1], gf[2] - ff[2]};
      for (double &component : ds) {
        component -= std::round(component);
      }
      const double wx = lx * ds[0] + xy * ds[1] + xz * ds[2];
      const double wy = ly * ds[1] + yz * ds[2];
      const double wz = lz * ds[2];
      const double wrapped = std::sqrt(wx * wx + wy * wy + wz * wz);
      double image = 1e300;
      for (int ia = -4; ia <= 4; ++ia) {
        for (int ib = -4; ib <= 4; ++ib) {
          for (int ic = -4; ic <= 4; ++ic) {
            const double dx = wx - (static_cast<double>(ia) * lx + static_cast<double>(ib) * xy +
                                    static_cast<double>(ic) * xz);
            const double dy = wy - (static_cast<double>(ib) * ly + static_cast<double>(ic) * yz);
            const double dz = wz - static_cast<double>(ic) * lz;
            image = std::min(image, std::sqrt(dx * dx + dy * dy + dz * dz));
          }
        }
      }
      if (wrapped > image + 1e-6) {
        wrap_disagrees += 1;
      }
      if (image < best) {
        best = image;
        best_index = nguest + f;
      }
      if (image <= cutoff) {
        contacts += 1;
      }
    }
    expect_distance[static_cast<std::size_t>(g)] = best;
    expect_index[static_cast<std::size_t>(g)] = best_index;
    expect_contacts[static_cast<std::size_t>(g)] = contacts;
  }
  REQUIRE(wrap_disagrees > 0);
  const auto report = pore::hostGuest(cloud, 1, 0, cutoff);
  REQUIRE(report.guest.size() == static_cast<std::size_t>(nguest));
  REQUIRE(report.framework_types.size() == 1);
  REQUIRE(report.framework_types[0] == 2);
  for (int g = 0; g < nguest; ++g) {
    REQUIRE(report.nearest[static_cast<std::size_t>(g)] == expect_index[static_cast<std::size_t>(g)]);
    REQUIRE_THAT(report.distance[static_cast<std::size_t>(g)],
                 Catch::Matchers::WithinAbs(expect_distance[static_cast<std::size_t>(g)], 1e-8));
    REQUIRE(report.contacts[static_cast<std::size_t>(g)] == expect_contacts[static_cast<std::size_t>(g)]);
  }
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
