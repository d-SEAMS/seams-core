#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <neighbours.hpp>
#include <numbers>
#include <pore.hpp>
#include <random>
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

struct Reduced2 {
  double b1[2];
  double b2[2];
  double cover;
  double inscribed;
};

Reduced2 reduceLattice(const double u[2], const double v[2]) {
  Reduced2 out{{u[0], u[1]}, {v[0], v[1]}, 0.0, 0.0};
  for (int step = 0; step < 32; ++step) {
    const double n1 = out.b1[0] * out.b1[0] + out.b1[1] * out.b1[1];
    const double n2 = out.b2[0] * out.b2[0] + out.b2[1] * out.b2[1];
    if (n1 > n2) {
      std::swap(out.b1[0], out.b2[0]);
      std::swap(out.b1[1], out.b2[1]);
    }
    const double base = out.b1[0] * out.b1[0] + out.b1[1] * out.b1[1];
    const double m = std::round((out.b1[0] * out.b2[0] + out.b1[1] * out.b2[1]) / base);
    if (m == 0.0) {
      if (out.b1[0] * out.b2[0] + out.b1[1] * out.b2[1] < 0.0) {
        out.b2[0] = -out.b2[0];
        out.b2[1] = -out.b2[1];
      }
      const double cross = std::abs(out.b1[0] * out.b2[1] - out.b1[1] * out.b2[0]);
      const double e1 = std::hypot(out.b1[0], out.b1[1]);
      const double e2 = std::hypot(out.b2[0], out.b2[1]);
      const double e3 = std::hypot(out.b1[0] - out.b2[0], out.b1[1] - out.b2[1]);
      out.cover = e1 * e2 * e3 / (2.0 * cross);
      out.inscribed = 0.5 * std::min(e1, std::min(e2, e3));
      return out;
    }
    out.b2[0] -= m * out.b1[0];
    out.b2[1] -= m * out.b1[1];
  }
  return out;
}

double coveringRadius(const double u[2], const double v[2]) { return reduceLattice(u, v).cover; }

// Polar area of the Voronoi cell of a reduced lattice clipped by a disk.
// The six Voronoi-relevant vectors are ±b1, ±b2, ±(b1−b2).
double voronoiPolar(const Reduced2 &lat, double radius) {
  const double w[6][2] = {
      {lat.b1[0], lat.b1[1]},
      {-lat.b1[0], -lat.b1[1]},
      {lat.b2[0], lat.b2[1]},
      {-lat.b2[0], -lat.b2[1]},
      {lat.b1[0] - lat.b2[0], lat.b1[1] - lat.b2[1]},
      {lat.b2[0] - lat.b1[0], lat.b2[1] - lat.b1[1]},
  };
  const int steps = 1 << 14;
  double sum = 0.0;
  for (int i = 0; i < steps; ++i) {
    const double theta = (static_cast<double>(i) + 0.5) * 2.0 * kPi / static_cast<double>(steps);
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    double rho = 1e300;
    for (const auto &edge : w) {
      const double along = c * edge[0] + s * edge[1];
      if (along > 0.0) {
        rho = std::min(rho, 0.5 * (edge[0] * edge[0] + edge[1] * edge[1]) / along);
      }
    }
    const double used = std::min(radius, rho);
    sum += 0.5 * used * used;
  }
  return sum * (2.0 * kPi / static_cast<double>(steps));
}

template <class VolumeAt>
void requireVolumeNearCover(double cover, double volume, VolumeAt &&volumeAt) {
  const double scales[] = {1.0 - 1e-15, 1.0 - 1e-13, 1.0 - 1e-12, 1.0};
  for (double scale : scales) {
    REQUIRE_THAT(volumeAt(cover * scale), Catch::Matchers::WithinRel(volume, 1e-8));
  }
  // The exit crossing and the entry crossing swap roles when a vertex
  // lies on the circle. That shows up at R and within a few ULPs of R.
  // The correct caps there are far smaller than 1e-8 of the cell.
  double radius = cover;
  for (int step = 0; step < 6; ++step) {
    radius = std::nextafter(radius, 0.0);
  }
  for (int step = -6; step <= 6; ++step) {
    REQUIRE_THAT(volumeAt(radius), Catch::Matchers::WithinRel(volume, 1e-8));
    radius = std::nextafter(radius, 1.0e300);
  }
}

std::vector<std::array<double, 2>> centredFace(const Face2 &face) {
  const double sign[4][2] = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}};
  std::vector<std::array<double, 2>> corners;
  for (const auto &s : sign) {
    corners.push_back({0.5 * (s[0] * face.eu[0] + s[1] * face.ev[0]),
                       0.5 * (s[0] * face.eu[1] + s[1] * face.ev[1])});
  }
  const double turn = corners[0][0] * corners[1][1] - corners[0][1] * corners[1][0];
  if (turn < 0.0) {
    std::reverse(corners.begin(), corners.end());
  }
  return corners;
}

double ringArea(const std::vector<std::array<double, 2>> &poly) {
  double area = 0.0;
  const int n = static_cast<int>(poly.size());
  for (int i = 0; i < n; ++i) {
    const auto &a = poly[static_cast<std::size_t>(i)];
    const auto &b = poly[static_cast<std::size_t>((i + 1) % n)];
    area += a[0] * b[1] - a[1] * b[0];
  }
  return 0.5 * std::abs(area);
}

// Midpoint polar integral of the disk clipped to a convex ring about the origin.
// 2^14 steps. The full disk is several times the face on the cells below, so
// a 1e-7 relative band still rejects that answer.
double polarClip(const std::vector<std::array<double, 2>> &poly, double radius) {
  const int n = static_cast<int>(poly.size());
  const int steps = 1 << 14;
  double sum = 0.0;
  for (int i = 0; i < steps; ++i) {
    const double theta = (static_cast<double>(i) + 0.5) * 2.0 * kPi / static_cast<double>(steps);
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    double rho = 1e300;
    for (int edge = 0; edge < n; ++edge) {
      const auto &a = poly[static_cast<std::size_t>(edge)];
      const auto &b = poly[static_cast<std::size_t>((edge + 1) % n)];
      const double ex = b[0] - a[0];
      const double ey = b[1] - a[1];
      const double den = c * ey - s * ex;
      if (!(std::abs(den) > 1e-18)) {
        continue;
      }
      const double t = (a[0] * ey - a[1] * ex) / den;
      const double u = (s * a[0] - c * a[1]) / den;
      if (t > 0.0 && u >= -1e-9 && u <= 1.0 + 1e-9) {
        rho = std::min(rho, t);
      }
    }
    const double r = std::min(radius, rho);
    sum += 0.5 * r * r;
  }
  return sum * (2.0 * kPi / static_cast<double>(steps));
}

double monteCarloShell(const Face2 &face, double radius, std::mt19937_64 &gen) {
  std::uniform_real_distribution<double> unit(-0.5, 0.5);
  const int samples = 40000;
  int inside = 0;
  for (int i = 0; i < samples; ++i) {
    if (nearestAxisRadius(face, unit(gen), unit(gen)) <= radius) {
      inside += 1;
    }
  }
  const double face_area = std::abs(face.eu[0] * face.ev[1] - face.eu[1] * face.ev[0]);
  return (static_cast<double>(inside) / static_cast<double>(samples)) * face_area * face.axial;
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

TEST_CASE("a shell at the covering radius is the cell, not the disk", "[pore]") {
  const double cells[][6] = {
      {20.0, 10.0 * std::sqrt(3.0), 8.0, 10.0, 0.0, 0.0},
      {10.0, 12.0, 14.0, 4.5, -3.0, 5.5},
      {10.0, 8.0, 12.0, 2.0, 1.0, -1.5},
  };
  const double centre[3] = {0.5, 0.5, 0.5};
  for (const auto &cell : cells) {
    auto cloud = triclinic(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    addAtom(cloud, 1, 0.1, 0.2, 0.3);
    const Face2 face = faceOf(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], 2);
    const double cover = coveringRadius(face.eu, face.ev);
    const double volume = cell[0] * cell[1] * cell[2];
    REQUIRE(cover > 0.0);
    // Vertices of the perpendicular lattice sit on this radius. A shell
    // at R, just inside it, and a few ULPs either side is the cell. The
    // dropped-root disk is pi r^2 L, which is 1.209 V, 1.274 V, and
    // 1.329 V on these three cells. Using the entry hit in place of the
    // exit hit is several percent off at these radii.
    requireVolumeNearCover(cover, volume, [&](double radius) {
      const auto profiles =
          pore::poreProfiles(cloud, 0, 1, 1, 1, radius, 2, nullptr, centre, 0, 0.0, nullptr);
      return profiles.radial.volume[0];
    });
  }
}

TEST_CASE("random cells match a Monte Carlo shell at the covering radius", "[pore]") {
  std::mt19937_64 gen(20261120);
  std::uniform_real_distribution<double> length(6.0, 16.0);
  std::uniform_real_distribution<double> tilt(-0.9, 0.9);
  const double centre[3] = {0.5, 0.5, 0.5};
  for (int trial = 0; trial < 24; ++trial) {
    const double lx = length(gen);
    const double ly = length(gen);
    const double lz = length(gen);
    const double xy = tilt(gen) * lx;
    const double xz = tilt(gen) * lx;
    const double yz = tilt(gen) * ly;
    auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
    addAtom(cloud, 1, 0.1, 0.2, 0.3);
    const Face2 face = faceOf(lx, ly, lz, xy, xz, yz, 2);
    const double cover = coveringRadius(face.eu, face.ev);
    const double volume = lx * ly * lz;
    REQUIRE(cover > 0.0);
    requireVolumeNearCover(cover, volume, [&](double radius) {
      const auto profiles =
          pore::poreProfiles(cloud, 0, 1, 1, 1, radius, 2, nullptr, centre, 0, 0.0, nullptr);
      return profiles.radial.volume[0];
    });
    // 0.6 R is inside the cell on a round lattice and clips a long one.
    // If it already covers the cell, half the inradius is a disk strictly
    // smaller than the face, so the sample is not the cell volume.
    const Reduced2 reduced = reduceLattice(face.eu, face.ev);
    double partial = 0.6 * cover;
    double polar = voronoiPolar(reduced, partial) * face.axial;
    if (!(polar < 0.85 * volume)) {
      partial = 0.5 * reduced.inscribed;
      polar = voronoiPolar(reduced, partial) * face.axial;
    }
    REQUIRE(polar > 0.0);
    REQUIRE(polar < 0.85 * volume);
    const auto partial_profiles =
        pore::poreProfiles(cloud, 0, 1, 1, 1, partial, 2, nullptr, centre, 0, 0.0, nullptr);
    const double shell = partial_profiles.radial.volume[0];
    // 2^14 midpoint steps. On seed 20261120 the largest relative gap is
    // 1.2e-8, and the shell is between 0.43 and 0.60 of the cell. 1e-6
    // still rejects a shell that returns the whole cell.
    REQUIRE_THAT(shell, Catch::Matchers::WithinRel(polar, 1e-6));
    const double sampled = monteCarloShell(face, partial, gen);
    // 40000 samples of a shell that covers only part of the face. The
    // largest relative gap on this seed is 0.010. This 0.02 band is a
    // stochastic sanity check (about 3.5 standard deviations), not the
    // accuracy gate: the deterministic polar comparison above is the 1e-6
    // gate and rejects biases that this band can hide.
    REQUIRE_THAT(shell, Catch::Matchers::WithinRel(sampled, 0.02));
  }
}

TEST_CASE("a centred face clipped at its long corner keeps the polygon", "[pore]") {
  auto check = [](double lx, double ly, double xy) {
    const Face2 face = faceOf(lx, ly, 10.0, xy, 0.0, 0.0, 2);
    const auto corners = centredFace(face);
    const double face_area = ringArea(corners);
    double far = 0.0;
    for (const auto &corner : corners) {
      far = std::max(far, std::hypot(corner[0], corner[1]));
    }
    REQUIRE(face_area > 0.0);
    REQUIRE(far > 0.0);
    requireVolumeNearCover(far, face_area, [&](double radius) {
      const double got = pore::clippedDiskArea(corners, radius, far * 2.0, 0.0);
      // At the long corner and a few ULPs around it the clip is the
      // polygon. The dropped-root result is the disk, and the entry hit
      // in place of the exit hit misses the polygon by several percent.
      REQUIRE_THAT(got, Catch::Matchers::WithinRel(polarClip(corners, radius), 1e-7));
      return got;
    });
  };
  check(13.543853, 15.493012, -9.327033);
  std::mt19937_64 gen(11);
  std::uniform_real_distribution<double> length(6.0, 16.0);
  std::uniform_real_distribution<double> tilt(-0.8, 0.8);
  for (int trial = 0; trial < 16; ++trial) {
    const double lx = length(gen);
    const double ly = length(gen);
    check(lx, ly, tilt(gen) * std::min(lx, ly));
  }
}

TEST_CASE("only axis -1 is automatic", "[pore]") {
  auto cloud = ortho(8.0, 9.0, 14.0);
  addAtom(cloud, 1, 1.0, 1.0, 1.0);
  REQUIRE_THROWS_AS(pore::poreFrame(cloud, 0, -7, nullptr, nullptr), std::out_of_range);
  REQUIRE_THROWS_AS(
      pore::poreProfiles(cloud, 0, 1, 1, 1, 1.0, -7, nullptr, nullptr, 0, 0.0, nullptr),
      std::out_of_range);
  REQUIRE_THROWS_AS(
      pore::poreProfiles(cloud, 0, 1, 1, 1, -1.0, 2, nullptr, nullptr, 0, 0.0, nullptr),
      std::out_of_range);
  const auto frame = pore::poreFrame(cloud, 0, -1, nullptr, nullptr);
  REQUIRE(frame.axis_index == 2);
}

TEST_CASE("a guest type equal to the framework type is rejected", "[pore]") {
  auto cloud = ortho(10.0, 10.0, 10.0);
  addAtom(cloud, 3, 1.0, 1.0, 1.0);
  addAtom(cloud, 4, 2.0, 2.0, 2.0);
  REQUIRE_THROWS_AS(pore::hostGuest(cloud, 3, 3, 1.0), std::invalid_argument);
  const auto report = pore::hostGuest(cloud, 3, 0, 5.0);
  REQUIRE(report.guest.size() == 1);
  REQUIRE(report.nearest.size() == 1);
  REQUIRE(report.nearest[0] == 1);
}

TEST_CASE("azimuth zero points along the next cell edge", "[pore]") {
  auto cloud = ortho(10.0, 12.0, 14.0);
  addAtom(cloud, 1, 6.0, 6.0, 7.0);
  const double centre[3] = {0.5, 0.5, 0.5};
  const auto profiles =
      pore::poreProfiles(cloud, 0, 1, 4, 1, 0.0, 2, nullptr, centre, 0, 0.0, nullptr);
  REQUIRE(profiles.azimuthal.coordinate.size() == 4);
  int occupied = -1;
  for (int bin = 0; bin < 4; ++bin) {
    if (profiles.azimuthal.number_density[static_cast<std::size_t>(bin)] > 0.0) {
      REQUIRE(occupied < 0);
      occupied = bin;
    }
  }
  REQUIRE(occupied == 2);
  REQUIRE_THAT(profiles.azimuthal.coordinate[2], Catch::Matchers::WithinAbs(kPi / 4.0, 1e-12));
}

TEST_CASE("minimum image does not read past a three length box", "[pore]") {
  Cloud cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const auto wrapped = gen::triclinicMinImage(cloud, 0.2, 0.0, 0.0, 9.9, 0.0, 0.0);
  REQUIRE_THAT(wrapped[0], Catch::Matchers::WithinAbs(0.3, 1e-12));
  REQUIRE_THAT(wrapped[1], Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(wrapped[2], Catch::Matchers::WithinAbs(0.0, 1e-12));
  cloud.box = {0.0, 10.0, 10.0};
  const auto raw = gen::triclinicMinImage(cloud, 1.0, 12.0, -6.0, 0.0, 0.0, 0.0);
  REQUIRE_THAT(raw[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(raw[1], Catch::Matchers::WithinAbs(2.0, 1e-12));
  REQUIRE_THAT(raw[2], Catch::Matchers::WithinAbs(4.0, 1e-12));
}

TEST_CASE("the Selling minimum image matches a brute-force search", "[pore]") {
  const double lx = 10.0;
  const double ly = 12.0;
  const double lz = 14.0;
  const double xy = 4.5;
  const double xz = -3.0;
  const double yz = 5.5;
  auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
  const auto box = gen::makeFracBox(cloud);
  REQUIRE(box.ok);
  int disagree = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 2; ++k) {
      const double s[3] = {(static_cast<double>(i) + 0.17) / 4.0, (static_cast<double>(j) + 0.41) / 4.0,
                           (static_cast<double>(k) + 0.63) / 2.0};
      const double t[3] = {(static_cast<double>((i + 1) % 3) + 0.5) / 3.0,
                           (static_cast<double>((j + 1) % 3) + 0.2) / 3.0,
                           (static_cast<double>(k % 3) + 0.8) / 3.0};
      const double a[3] = {lx * s[0] + xy * s[1] + xz * s[2], ly * s[1] + yz * s[2], lz * s[2]};
      const double b[3] = {lx * t[0] + xy * t[1] + xz * t[2], ly * t[1] + yz * t[2], lz * t[2]};
      double ds[3] = {s[0] - t[0], s[1] - t[1], s[2] - t[2]};
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
            const double dx = wx - (ia * lx + ib * xy + ic * xz);
            const double dy = wy - (ib * ly + ic * yz);
            const double dz = wz - ic * lz;
            image = std::min(image, std::sqrt(dx * dx + dy * dy + dz * dz));
          }
        }
      }
      if (wrapped > image + 1e-6) {
        disagree += 1;
      }
      const auto wrapped_vec = gen::fracDelta(box, a[0], a[1], a[2], b[0], b[1], b[2]);
      const std::array<double, 3> raw = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
      const auto searched = gen::mic::euclideanIfShorter(box.lx, box.ly, box.lz, box.xy, box.xz,
                                                        box.yz, raw, wrapped_vec);
      const double found =
          std::sqrt(searched[0] * searched[0] + searched[1] * searched[1] + searched[2] * searched[2]);
      REQUIRE_THAT(found, Catch::Matchers::WithinAbs(image, 1e-8));
      const auto via = gen::triclinicMinImage(cloud, a[0], a[1], a[2], b[0], b[1], b[2]);
      const double library =
          std::sqrt(via[0] * via[0] + via[1] * via[1] + via[2] * via[2]);
      REQUIRE_THAT(library, Catch::Matchers::WithinAbs(image, 1e-8));
      }
    }
  }
  REQUIRE(disagree > 0);
}

double bruteImage(double lx, double ly, double lz, double xy, double xz, double yz, const double d[3],
                  int n) {
  double best = 1e300;
  for (int i = -n; i <= n; ++i) {
    for (int j = -n; j <= n; ++j) {
      for (int k = -n; k <= n; ++k) {
        const double x = d[0] + static_cast<double>(i) * lx + static_cast<double>(j) * xy +
                         static_cast<double>(k) * xz;
        const double y = d[1] + static_cast<double>(j) * ly + static_cast<double>(k) * yz;
        const double z = d[2] + static_cast<double>(k) * lz;
        best = std::min(best, std::sqrt(x * x + y * y + z * z));
      }
    }
  }
  return best;
}

struct LongVec {
  long double x;
  long double y;
  long double z;
};

long double longDot(LongVec a, LongVec b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

LongVec longSubtract(LongVec a, LongVec b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

LongVec longScale(long double scale, LongVec value) {
  return {scale * value.x, scale * value.y, scale * value.z};
}

void reduceReferenceBasis(LongVec basis[3]) {
  for (;;) {
    LongVec orthogonal[3];
    long double mu[3][3] = {};
    long double norm[3] = {};
    for (int i = 0; i < 3; ++i) {
      orthogonal[i] = basis[i];
      for (int j = 0; j < i; ++j) {
        mu[i][j] = longDot(basis[i], orthogonal[j]) / norm[j];
        orthogonal[i] =
            longSubtract(orthogonal[i], longScale(mu[i][j], orthogonal[j]));
      }
      norm[i] = longDot(orthogonal[i], orthogonal[i]);
    }
    bool changed = false;
    for (int k = 1; k < 3 && !changed; ++k) {
      for (int j = k - 1; j >= 0; --j) {
        if (std::abs(mu[k][j]) <= 0.5L + 1e-12L) {
          continue;
        }
        const long double q = std::round(mu[k][j]);
        if (q != 0.0L) {
          basis[k] = longSubtract(basis[k], longScale(q, basis[j]));
          changed = true;
          break;
        }
      }
      if (!changed &&
          norm[k] < (0.99L - mu[k][k - 1] * mu[k][k - 1]) * norm[k - 1]) {
        std::swap(basis[k], basis[k - 1]);
        changed = true;
      }
    }
    if (!changed) {
      return;
    }
  }
}

double exactImage(double lx, double ly, double lz, double xy, double xz,
                  double yz, const double d[3]) {
  LongVec basis[3] = {{lx, 0.0L, 0.0L}, {xy, ly, 0.0L}, {xz, yz, lz}};
  reduceReferenceBasis(basis);
  LongVec orthogonal[3];
  long double mu[3][3] = {};
  long double norm[3] = {};
  for (int i = 0; i < 3; ++i) {
    orthogonal[i] = basis[i];
    for (int j = 0; j < i; ++j) {
      mu[i][j] = longDot(basis[i], orthogonal[j]) / norm[j];
      orthogonal[i] =
          longSubtract(orthogonal[i], longScale(mu[i][j], orthogonal[j]));
    }
    norm[i] = longDot(orthogonal[i], orthogonal[i]);
  }
  const LongVec target = {d[0], d[1], d[2]};
  LongVec babai = target;
  for (int i = 2; i >= 0; --i) {
    const long double q = std::round(longDot(babai, orthogonal[i]) / norm[i]);
    babai = longSubtract(babai, longScale(q, basis[i]));
  }
  long double best2 = longDot(babai, babai);
  long double coordinate[3];
  for (int i = 0; i < 3; ++i) {
    coordinate[i] = longDot(target, orthogonal[i]) / norm[i];
  }
  const long double width2 = std::sqrt(best2 / norm[2]);
  for (long long i2 = static_cast<long long>(std::ceil(coordinate[2] - width2));
       i2 <= static_cast<long long>(std::floor(coordinate[2] + width2)); ++i2) {
    const long double error2 = coordinate[2] - static_cast<long double>(i2);
    const long double radius2 = error2 * error2 * norm[2];
    if (radius2 > best2) {
      continue;
    }
    const long double centre1 =
        coordinate[1] - static_cast<long double>(i2) * mu[2][1];
    const long double width1 = std::sqrt((best2 - radius2) / norm[1]);
    for (long long i1 = static_cast<long long>(std::ceil(centre1 - width1));
         i1 <= static_cast<long long>(std::floor(centre1 + width1)); ++i1) {
      const long double error1 = centre1 - static_cast<long double>(i1);
      const long double radius1 = radius2 + error1 * error1 * norm[1];
      if (radius1 > best2) {
        continue;
      }
      const long double centre0 = coordinate[0] -
                                  static_cast<long double>(i2) * mu[2][0] -
                                  static_cast<long double>(i1) * mu[1][0];
      const long double width0 = std::sqrt((best2 - radius1) / norm[0]);
      for (long long i0 = static_cast<long long>(std::ceil(centre0 - width0));
           i0 <= static_cast<long long>(std::floor(centre0 + width0)); ++i0) {
        const LongVec lattice = {static_cast<long double>(i0) * basis[0].x +
                                     static_cast<long double>(i1) * basis[1].x +
                                     static_cast<long double>(i2) * basis[2].x,
                                 static_cast<long double>(i0) * basis[0].y +
                                     static_cast<long double>(i1) * basis[1].y +
                                     static_cast<long double>(i2) * basis[2].y,
                                 static_cast<long double>(i0) * basis[0].z +
                                     static_cast<long double>(i1) * basis[1].z +
                                     static_cast<long double>(i2) * basis[2].z};
        const LongVec delta = longSubtract(target, lattice);
        best2 = std::min(best2, longDot(delta, delta));
      }
    }
  }
  return std::sqrt(static_cast<double>(best2));
}

TEST_CASE("minimum image on a long tilted cell matches brute force", "[pore]") {
  const double channels[][6] = {
      {12.0, 12.0, 100.0, 0.5, 0.0, 0.0},
      {12.0, 12.0, 400.0, 0.5, 0.0, 0.0},
      {12.0, 12.0, 1000.0, 0.5, 0.0, 0.0},
      {10.0, 10.0, 2000.0, 1e-6, 0.0, 0.0},
      {10000.0, 12.0, 12.0, 3.0, 3.0, 0.0},
      {1000000.0, 1.0, 1.0, 0.5, 0.5, 0.5},
  };
  std::mt19937_64 namedGen(4242);
  std::uniform_real_distribution<double> fractional(-0.5, 0.5);
  for (const auto &cell : channels) {
    auto cloud =
        triclinic(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
    for (int pair = 0; pair < 64; ++pair) {
      const double sf[3] = {fractional(namedGen), fractional(namedGen),
                            fractional(namedGen)};
      const double tf[3] = {fractional(namedGen), fractional(namedGen),
                            fractional(namedGen)};
      const double a[3] = {cell[0] * sf[0] + cell[3] * sf[1] + cell[4] * sf[2],
                           cell[1] * sf[1] + cell[5] * sf[2], cell[2] * sf[2]};
      const double b[3] = {cell[0] * tf[0] + cell[3] * tf[1] + cell[4] * tf[2],
                           cell[1] * tf[1] + cell[5] * tf[2], cell[2] * tf[2]};
      const auto got =
          gen::triclinicMinImage(cloud, b[0], b[1], b[2], a[0], a[1], a[2]);
      const double found =
          std::sqrt(got[0] * got[0] + got[1] * got[1] + got[2] * got[2]);
      const double d[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
      const double ref =
          exactImage(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], d);
      REQUIRE_THAT(found, Catch::Matchers::WithinAbs(ref, 1e-8));
    }
  }

  // The gate's seed and cell count. Restricting the tilts to LAMMPS limits
  // makes every generated cell physically representable while retaining the
  // large aspect ratios that exposed the fractional-start bug.
  std::mt19937_64 gen(4242);
  std::uniform_real_distribution<double> logLength(-3.0, 3.0);
  std::uniform_real_distribution<double> half(-0.5, 0.5);
  int longer = 0;
  for (int trial = 0; trial < 4000; ++trial) {
    const double lx = std::pow(10.0, logLength(gen));
    const double ly = std::pow(10.0, logLength(gen));
    const double lz = std::pow(10.0, logLength(gen));
    const double xy = half(gen) * lx;
    const double xz = half(gen) * lx;
    const double yz = half(gen) * ly;
    Cloud cloud;
    REQUIRE(nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, 0.4, -1.2,
                                          2.0, lx, ly, lz, xy, xz, yz));
    const double sf[3] = {half(gen), half(gen), half(gen)};
    const double tf[3] = {half(gen), half(gen), half(gen)};
    const double left[3] = {0.4 + lx * sf[0] + xy * sf[1] + xz * sf[2],
                            -1.2 + ly * sf[1] + yz * sf[2], 2.0 + lz * sf[2]};
    const double right[3] = {0.4 + lx * tf[0] + xy * tf[1] + xz * tf[2],
                             -1.2 + ly * tf[1] + yz * tf[2], 2.0 + lz * tf[2]};
    const auto got = gen::triclinicMinImage(cloud, left[0], left[1], left[2],
                                            right[0], right[1], right[2]);
    const double found =
        std::sqrt(got[0] * got[0] + got[1] * got[1] + got[2] * got[2]);
    const double d[3] = {left[0] - right[0], left[1] - right[1],
                         left[2] - right[2]};
    const double ref = exactImage(lx, ly, lz, xy, xz, yz, d);
    const auto wrapped =
        gen::fracDelta(gen::makeFracBox(cloud), left[0], left[1], left[2],
                       right[0], right[1], right[2]);
    const double wrap =
        std::sqrt(wrapped[0] * wrapped[0] + wrapped[1] * wrapped[1] +
                  wrapped[2] * wrapped[2]);
    if (wrap > ref + 1e-9) {
      longer += 1;
    }
    REQUIRE_THAT(found, Catch::Matchers::WithinAbs(ref, 1e-8));
  }

  // Preserve the original over-tilted seed set as well as the gate sweep.
  std::mt19937_64 legacyGen(77);
  std::uniform_real_distribution<double> length(5.0, 15.0);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> mild(-0.5, 0.5);
  std::uniform_real_distribution<double> steep(-2.5, 2.5);
  for (int trial = 0; trial < 400; ++trial) {
    const bool over = trial % 2 == 0;
    const double lx = length(legacyGen);
    const double ly = length(legacyGen);
    const double lz = length(legacyGen);
    const double xy = (over ? steep(legacyGen) : mild(legacyGen)) * lx;
    const double xz = (over ? steep(legacyGen) : mild(legacyGen)) * lx;
    const double yz = (over ? steep(legacyGen) : mild(legacyGen)) * ly;
    Cloud cloud;
    REQUIRE(nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, 0.4, -1.2,
                                          2.0, lx, ly, lz, xy, xz, yz));
    for (int pair = 0; pair < 25; ++pair) {
      auto place = [&](double s0, double s1, double s2, double out[3]) {
        out[0] = 0.4 + lx * s0 + xy * s1 + xz * s2;
        out[1] = -1.2 + ly * s1 + yz * s2;
        out[2] = 2.0 + lz * s2;
      };
      double left[3];
      double right[3];
      place(unit(legacyGen) * 3.0 - 1.0, unit(legacyGen) * 3.0 - 1.0,
            unit(legacyGen) * 3.0 - 1.0, left);
      place(unit(legacyGen), unit(legacyGen), unit(legacyGen), right);
      const auto got = gen::triclinicMinImage(cloud, left[0], left[1], left[2],
                                              right[0], right[1], right[2]);
      const double found =
          std::sqrt(got[0] * got[0] + got[1] * got[1] + got[2] * got[2]);
      const double d[3] = {left[0] - right[0], left[1] - right[1],
                           left[2] - right[2]};
      const double ref = bruteImage(lx, ly, lz, xy, xz, yz, d, 9);
      REQUIRE_THAT(found, Catch::Matchers::WithinAbs(ref, 1e-8));
    }
  }
  REQUIRE(longer > 0);
}

TEST_CASE("minimum image rejects a numerically degenerate tilted cell",
          "[pore]") {
  const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
  const double ly = std::pow(phi, -120.0);
  const std::array<double, 3> displacement = {0.4, 0.4 * ly, 0.4};
  REQUIRE_THROWS_AS(gen::mic::euclideanIfShorter(1.0, ly, 1.0, 0.5, 0.25, 0.0,
                                                 displacement, displacement),
                    std::runtime_error);

  const std::array<double, 3> finite = {0.25, 0.25, 0.25};
  const std::array<double, 3> nonfinite = {
      std::numeric_limits<double>::infinity(), 0.0, 0.0};
  REQUIRE_THROWS_AS(gen::mic::euclideanIfShorter(1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                                                 finite, nonfinite),
                    std::runtime_error);

  const std::array<double, 3> unresolved = {1.0e16, 0.0, 0.0};
  REQUIRE_THROWS_AS(gen::mic::euclideanIfShorter(1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                                                 unresolved, unresolved),
                    std::runtime_error);
}

TEST_CASE("minimum image is exact on a 1e8-aspect LAMMPS-legal cell", "[pore]") {
  // The first pair used to stop the slicer with "binary64 limits". Both pairs
  // must match the reference to within a few ulps of the wrapped length.
  const std::array<double, 3> pairs[] = {
      {-4859743.2640920915, -0.0050198905816330854, -0.031675054668108597},
      {2500000.0, 0.02, -0.03},
  };
  for (const auto &w : pairs) {
    const auto got =
        gen::mic::euclideanIfShorter(1.0e7, 0.1, 0.1, 0.025, 0.025, 0.0, w, w);
    const double d[3] = {w[0], w[1], w[2]};
    const double ref = exactImage(1.0e7, 0.1, 0.1, 0.025, 0.025, 0.0, d);
    const double wrapped = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    REQUIRE_THAT(std::sqrt(got[0] * got[0] + got[1] * got[1] + got[2] * got[2]),
                 Catch::Matchers::WithinAbs(
                     ref, 16.0 * std::numeric_limits<double>::epsilon() * wrapped));
  }
}

TEST_CASE("minimum image basis cache is not poisoned by a throw", "[pore]") {
  // Fill all four per-thread cache slots with valid tilted cells.
  const double fillers[][6] = {
      {11.0, 10.0, 10.0, 1.0, 0.0, 0.0},
      {12.0, 10.0, 10.0, 1.0, 0.0, 0.0},
      {13.0, 10.0, 10.0, 1.0, 0.0, 0.0},
      {14.0, 10.0, 10.0, 1.0, 0.0, 0.0},
  };
  const std::array<double, 3> probe = {4.0, 4.0, 4.0};
  for (const auto &cell : fillers) {
    const auto res = gen::mic::euclideanIfShorter(
        cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], probe, probe);
    REQUIRE(std::isfinite(res[0]));
  }

  // Trigger a throw inside obtuseSuperbasis (lagrangeReduce coefficient >= 2^52)
  // and verify that repeating the identical call throws again instead of
  // serving the evicted slot's previous superbasis.
  const std::array<double, 3> vec = {0.4, 0.4, 0.4};
  REQUIRE_THROWS_AS(gen::mic::euclideanIfShorter(1.0, 1.0, 1.0, 5.0e15, 0.0,
                                                 0.0, vec, vec),
                    std::runtime_error);
  REQUIRE_THROWS_AS(gen::mic::euclideanIfShorter(1.0, 1.0, 1.0, 5.0e15, 0.0,
                                                 0.0, vec, vec),
                    std::runtime_error);

  // Also test a range-overflow throw inside obtuseSuperbasis.
  REQUIRE_THROWS_AS(gen::mic::euclideanIfShorter(1.0e200, 1.0, 1.0, 0.0, 0.0,
                                                 0.0, vec, vec),
                    std::runtime_error);
  REQUIRE_THROWS_AS(gen::mic::euclideanIfShorter(1.0e200, 1.0, 1.0, 0.0, 0.0,
                                                 0.0, vec, vec),
                    std::runtime_error);

  // Subsequent valid cell still produces the exact lattice image.
  const std::array<double, 3> raw = {431.0, 5.0, 9.0};
  const auto cloud = triclinic(1000.0, 12.0, 20.4, 500.0, -350.0, 0.0);
  const auto box = gen::makeFracBox(cloud);
  const auto wrapped = gen::fracDelta(box, raw[0], raw[1], raw[2], 0.0, 0.0, 0.0);
  const auto got = gen::mic::euclideanIfShorter(1000.0, 12.0, 20.4, 500.0,
                                                -350.0, 0.0, raw, wrapped);
  const double d[3] = {wrapped[0], wrapped[1], wrapped[2]};
  const double ref = exactImage(1000.0, 12.0, 20.4, 500.0, -350.0, 0.0, d);
  REQUIRE_THAT(std::sqrt(got[0] * got[0] + got[1] * got[1] + got[2] * got[2]),
               Catch::Matchers::WithinAbs(ref, 1e-10));
}

TEST_CASE("minimum image handles ordinary cells without spurious throws",
          "[pore]") {
  const double s3 = std::sqrt(3.0);
  const double cells[][6] = {
      {20.0, 20.0 * s3 / 2.0, 20.0 * std::sqrt(2.0 / 3.0), 10.0, 10.0,
       20.0 / (2.0 * s3)},
      {2.5527, 2.5527 * s3 / 2.0, 2.5527 * std::sqrt(2.0 / 3.0), 2.5527 / 2.0,
       2.5527 / 2.0, 2.5527 / (2.0 * s3)},
      {14.142135623730951, 12.24744871391589, 11.547005383792516,
       7.0710678118654755, 7.0710678118654755, 4.08248290463863},
      {20.0, 20.0, 20.0, 1.0e-16, 0.0, 0.0},
      {20.0, 20.0, 20.0, 1.0e-12, 0.0, 0.0},
      {20.0, 20.0, 20.0, 1.0e-12, 1.0e-12, 1.0e-12},
      {1000.0, 12.0, 20.4, 500.0, -350.0, 0.0},
      {1000.0, 12.0, 20.4, 500.0, -350.0, -6.0},
      {25.0, 21.650635094610966, 30.0, 12.5, 3.3, 1.1},
      {25.0, 21.65063509461097, 30.0, 12.5, 3.3, 1.1},
  };
  std::mt19937_64 gen(20261121);
  std::uniform_real_distribution<double> u(-0.5, 0.5);
  for (const auto &c : cells) {
    auto cloud = triclinic(c[0], c[1], c[2], c[3], c[4], c[5]);
    const auto box = gen::makeFracBox(cloud);
    REQUIRE(box.ok);
    for (int pair = 0; pair < 32; ++pair) {
      const double f[3] = {pair == 0 ? 0.49 : u(gen),
                           pair == 0 ? 0.47 : u(gen),
                           pair == 0 ? 0.46 : u(gen)};
      const double p[3] = {c[0] * f[0] + c[3] * f[1] + c[4] * f[2],
                           c[1] * f[1] + c[5] * f[2], c[2] * f[2]};
      const auto wrapped = gen::fracDelta(box, p[0], p[1], p[2], 0.0, 0.0, 0.0);
      const std::array<double, 3> raw = {p[0], p[1], p[2]};
      const auto gotMic = gen::mic::euclideanIfShorter(
          c[0], c[1], c[2], c[3], c[4], c[5], raw, wrapped);
      const auto gotApi =
          gen::triclinicMinImage(cloud, p[0], p[1], p[2], 0.0, 0.0, 0.0);
      const double d[3] = {wrapped[0], wrapped[1], wrapped[2]};
      const double ref = exactImage(c[0], c[1], c[2], c[3], c[4], c[5], d);
      const double normMic = std::sqrt(gotMic[0] * gotMic[0] +
                                       gotMic[1] * gotMic[1] +
                                       gotMic[2] * gotMic[2]);
      const double normApi = std::sqrt(gotApi[0] * gotApi[0] +
                                       gotApi[1] * gotApi[1] +
                                       gotApi[2] * gotApi[2]);
      REQUIRE_THAT(normMic, Catch::Matchers::WithinAbs(ref, 1e-10));
      REQUIRE_THAT(normApi, Catch::Matchers::WithinAbs(ref, 1e-10));
      // Check lattice membership of gotMic - wrapped in H.
      const double dv[3] = {gotMic[0] - wrapped[0], gotMic[1] - wrapped[1],
                            gotMic[2] - wrapped[2]};
      const double cz = dv[2] / c[2];
      const double cy = (dv[1] - c[5] * cz) / c[1];
      const double cx = (dv[0] - c[3] * cy - c[4] * cz) / c[0];
      REQUIRE(std::abs(cx - std::round(cx)) < 1e-9);
      REQUIRE(std::abs(cy - std::round(cy)) < 1e-9);
      REQUIRE(std::abs(cz - std::round(cz)) < 1e-9);
    }
  }
}

// Exact minimum image over the upper-triangular cell, enumerating k, then j,
// then the two nearest i, inside the radius `bound`. Returns -1 when no
// lattice point lies inside that radius.
long double triangularImage(const double c[6], const double d[3], long double bound) {
  const long double lx = c[0];
  const long double ly = c[1];
  const long double lz = c[2];
  const long double xy = c[3];
  const long double xz = c[4];
  const long double yz = c[5];
  long double best = bound * bound;
  bool found = false;
  const long double radius = std::sqrt(best);
  const auto k0 = static_cast<long long>(std::floor((d[2] - radius) / lz)) - 1;
  const auto k1 = static_cast<long long>(std::ceil((d[2] + radius) / lz)) + 1;
  for (long long k = k0; k <= k1; ++k) {
    const long double z = d[2] - static_cast<long double>(k) * lz;
    const long double restZ = best - z * z;
    if (restZ < 0.0L) {
      continue;
    }
    const long double yc = d[1] - static_cast<long double>(k) * yz;
    const long double ry = std::sqrt(restZ);
    const auto j0 = static_cast<long long>(std::floor((yc - ry) / ly)) - 1;
    const auto j1 = static_cast<long long>(std::ceil((yc + ry) / ly)) + 1;
    for (long long j = j0; j <= j1; ++j) {
      const long double y = yc - static_cast<long double>(j) * ly;
      if (restZ - y * y < 0.0L) {
        continue;
      }
      const long double xc = d[0] - static_cast<long double>(j) * xy -
                             static_cast<long double>(k) * xz;
      const auto ic = static_cast<long long>(std::floor(xc / lx));
      for (long long i = ic; i <= ic + 1; ++i) {
        const long double x = xc - static_cast<long double>(i) * lx;
        const long double n = x * x + y * y + z * z;
        if (n <= best) {
          best = n;
          found = true;
        }
      }
    }
  }
  return found ? std::sqrt(best) : -1.0L;
}

// Library image of the fractional wrap `w` on cell `c`, checked against the
// triangular reference: the call returns, the result is a lattice translate
// of `w`, and its length matches to `rel`.
void requireExactImage(const double c[6], const std::array<double, 3> &w, double rel) {
  const auto got = gen::mic::euclideanIfShorter(c[0], c[1], c[2], c[3], c[4], c[5], w, w);
  const double found = std::sqrt(got[0] * got[0] + got[1] * got[1] + got[2] * got[2]);
  const double d[3] = {w[0], w[1], w[2]};
  const long double ref = triangularImage(c, d, static_cast<long double>(found) * (1.0L + 1e-12L));
  REQUIRE(ref >= 0.0L);
  REQUIRE_THAT(found, Catch::Matchers::WithinRel(static_cast<double>(ref), rel));
  const double dv[3] = {got[0] - w[0], got[1] - w[1], got[2] - w[2]};
  const double cz = dv[2] / c[2];
  const double cy = (dv[1] - c[5] * cz) / c[1];
  const double cx = (dv[0] - c[3] * cy - c[4] * cz) / c[0];
  REQUIRE(std::abs(cx - std::round(cx)) < 1e-9);
  REQUIRE(std::abs(cy - std::round(cy)) < 1e-9);
  REQUIRE(std::abs(cz - std::round(cz)) < 1e-9);
}

TEST_CASE("minimum image does not stall on cells where the slicer used to throw",
          "[pore]") {
  // Cells and wraps printed by the round-6 gate (nine significant digits).
  // Each threw "minimum-image slicer reached binary64 limits".
  const double legal[][9] = {
      {14.745736, 71.0592255, 3.85588635, 4.0013901, 3.48405873, 32.153036, 3.74272441,
       51.6061307, 1.92794318},
      {17.2569651, 4.67960859, 1.8394437, -4.97263683, -5.72362751, 0.627733896, 3.2803504,
       2.65367124, 0.91972185},
      {663.227964, 189.35351, 25.6360714, -285.603005, -145.440404, 74.3425649, 116.092277,
       131.848038, 12.8180357},
      {13.5686657, 3.46816642, 0.58590604, 4.27697286, 3.17231146, -1.25982259, 6.78433287,
       0.0, 0.0},
      {3.70116856, 0.868685911, 11.5064988, -0.23454385, -0.965981861, 0.395936338,
       -0.117271925, 0.434342955, 0.0},
      {186.992411, 313.299545, 10.4083736, 44.9603314, 42.4777125, 108.680464, 114.735062,
       54.3402322, 5.2041868},
      {0.576046042, 0.959331047, 0.711803081, -0.145219671, 0.106501407, -0.0340919453,
       -0.307382152, 0.462619551, 0.35590154},
      {11.9522972, 10.3268806, 0.687091899, -1.75540131, 2.32532548, 3.01871018, 5.97614862,
       0.0, 0.0},
      {2.07161518, 0.0703767252, 0.226021188, -0.649221026, 0.529175895, 0.0335878335,
       1.03580759, 0.0, 0.0},
      {58.196444, 53.3197415, 54.9325006, 29.098222, -29.098222, -26.6598708, 43.647333,
       26.6598708, 0.0},
  };
  for (const auto &row : legal) {
    requireExactImage(row, {row[6], row[7], row[8]}, 1e-14);
  }
  // The gate's brute-force 30.3714560195 is for the unrounded corner cell.
  const auto corner = gen::mic::euclideanIfShorter(
      58.196444, 53.3197415, 54.9325006, 29.098222, -29.098222, -26.6598708,
      {43.647333, 26.6598708, 0.0}, {43.647333, 26.6598708, 0.0});
  REQUIRE_THAT(std::sqrt(corner[0] * corner[0] + corner[1] * corner[1] + corner[2] * corner[2]),
               Catch::Matchers::WithinAbs(30.3714560195, 1e-6));

  // A 1e-14 tilt on a 10000 x 12 x 12 cell, where most of these wraps are
  // already the minimum image.
  const double tiny[6] = {10000.0, 12.0, 12.0, 1e-14, 1e-14, 0.0};
  const std::array<double, 3> wraps[] = {
      {5000.0, 0.0, 6.0},  {5000.0, 5.88, 6.0},  {5000.0, 6.0, 0.0},   {5000.0, 6.0, 6.0},
      {-4900.0, 0.0, 6.0}, {-4900.0, 6.0, 0.0},  {-4900.0, 6.0, 6.0},  {-2500.0, 0.0, 6.0},
      {-2500.0, 6.0, 0.0}, {-2500.0, 6.0, 6.0},  {2500.0, 5.88, 6.0},  {2500.0, 6.0, 5.88},
  };
  for (const auto &w : wraps) {
    requireExactImage(tiny, w, 1e-14);
  }
  const auto first = gen::mic::euclideanIfShorter(tiny[0], tiny[1], tiny[2], tiny[3], tiny[4],
                                                  tiny[5], wraps[0], wraps[0]);
  REQUIRE_THAT(std::sqrt(first[0] * first[0] + first[1] * first[1] + first[2] * first[2]),
               Catch::Matchers::WithinRel(5000.0036, 1e-9));
}

TEST_CASE("minimum image matches a triangular brute force on random legal cells",
          "[pore]") {
  // 1200 LAMMPS-legal cells with edge aspect at most 40; a quarter put one
  // tilt on the +-edge/2 limit and a quarter put all three there. Each cell
  // probes the seven face, edge, and body half-lattice points plus 17 random
  // fractions.
  std::mt19937_64 gen(20261120);
  std::uniform_real_distribution<double> logScale(-1.0, 3.0);
  std::uniform_real_distribution<double> logAspect(0.0, std::log10(40.0));
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> half(-0.5, 0.5);
  const double halves[][3] = {{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5},
                              {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
                              {0.5, 0.5, 0.5}};
  int threw = 0;
  long pairs = 0;
  for (int cell = 0; cell < 1200; ++cell) {
    const double base = std::pow(10.0, logScale(gen));
    double edge[3] = {base, base * std::pow(10.0, logAspect(gen) * unit(gen)),
                      base * std::pow(10.0, logAspect(gen))};
    std::shuffle(edge, edge + 3, gen);
    double c[6] = {edge[0], edge[1], edge[2], half(gen) * edge[0], half(gen) * edge[0],
                   half(gen) * edge[1]};
    if (cell % 4 == 1) {
      c[3] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[0];
    }
    if (cell % 4 == 2) {
      c[3] = 0.5 * edge[0];
      c[4] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[0];
      c[5] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[1];
    }
    for (int p = 0; p < 24; ++p) {
      double f[3];
      for (int k = 0; k < 3; ++k) {
        f[k] = p < 7 ? halves[p][k] : half(gen);
      }
      const std::array<double, 3> w = {c[0] * f[0] + c[3] * f[1] + c[4] * f[2],
                                       c[1] * f[1] + c[5] * f[2], c[2] * f[2]};
      ++pairs;
      try {
        requireExactImage(c, w, 1e-14);
      } catch (const std::runtime_error &) {
        ++threw;
      }
    }
  }
  REQUIRE(pairs == 1200 * 24);
  REQUIRE(threw == 0);
}

// sum m[i] * x[i] with error-free products and a compensated sum (Dot2).
long double dot2(const double m[4], const double x[4]) {
  double s = 0.0;
  double e = 0.0;
  for (int i = 0; i < 4; ++i) {
    const double p = m[i] * x[i];
    const double pe = std::fma(m[i], x[i], -p);
    const double t = s + p;
    const double b = t - s;
    e += pe + ((s - (t - b)) + (p - b));
    s = t;
  }
  return static_cast<long double>(s) + e;
}

// w - (n0 a + n1 b + n2 c) on the upper-triangular cell, each component to
// about one rounding.
std::array<long double, 3> latticeImage(const double c[6], const double w[3],
                                        const long long n[3]) {
  const double m[4] = {1.0, -static_cast<double>(n[0]), -static_cast<double>(n[1]),
                       -static_cast<double>(n[2])};
  const double xs[4] = {w[0], c[0], c[3], c[4]};
  const double ys[4] = {w[1], 0.0, c[1], c[5]};
  const double zs[4] = {w[2], 0.0, 0.0, c[2]};
  return {dot2(m, xs), dot2(m, ys), dot2(m, zs)};
}

template <class T> long double lengthOf(const std::array<T, 3> &v) {
  return std::sqrt(static_cast<long double>(v[0]) * v[0] +
                   static_cast<long double>(v[1]) * v[1] +
                   static_cast<long double>(v[2]) * v[2]);
}

// Shortest |w - n H| no longer than `radius`, or -1 when there is none or the
// search is too wide. LLL with an integer transform picks the enumeration
// basis; every candidate is scored from its integer coefficients in H.
long double coefficientImage(const double c[6], const double w[3], long double radius) {
  long long u[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  long double b[3][3];
  const auto rebuild = [&](int k) {
    const double zero[3] = {0.0, 0.0, 0.0};
    const long long n[3] = {-u[k][0], -u[k][1], -u[k][2]};
    const auto v = latticeImage(c, zero, n);
    for (int m = 0; m < 3; ++m) {
      b[k][m] = v[m];
    }
  };
  const auto dot = [](const long double *p, const long double *q) {
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
  };
  long double gs[3][3];
  long double mu[3][3] = {};
  long double norm[3];
  const auto gram = [&]() {
    for (int i = 0; i < 3; ++i) {
      for (int m = 0; m < 3; ++m) {
        gs[i][m] = b[i][m];
      }
      for (int j = 0; j < i; ++j) {
        mu[i][j] = dot(b[i], gs[j]) / norm[j];
        for (int m = 0; m < 3; ++m) {
          gs[i][m] -= mu[i][j] * gs[j][m];
        }
      }
      norm[i] = dot(gs[i], gs[i]);
    }
  };
  for (int k = 0; k < 3; ++k) {
    rebuild(k);
  }
  for (int iter = 0; iter < 100000; ++iter) {
    gram();
    bool changed = false;
    for (int k = 1; k < 3 && !changed; ++k) {
      for (int j = k - 1; j >= 0; --j) {
        if (std::abs(mu[k][j]) <= 0.5L + 1e-9L) {
          continue;
        }
        const long long q = std::llround(mu[k][j]);
        for (int m = 0; m < 3; ++m) {
          u[k][m] -= q * u[j][m];
        }
        rebuild(k);
        changed = true;
        break;
      }
      if (!changed && norm[k] < (0.99L - mu[k][k - 1] * mu[k][k - 1]) * norm[k - 1]) {
        std::swap(u[k], u[k - 1]);
        rebuild(k);
        rebuild(k - 1);
        changed = true;
      }
    }
    if (!changed) {
      break;
    }
  }
  gram();
  const long double t[3] = {w[0], w[1], w[2]};
  long double coordinate[3];
  for (int i = 0; i < 3; ++i) {
    coordinate[i] = dot(t, gs[i]) / norm[i];
  }
  long double best = -1.0L;
  const long double best2 = radius * radius;
  long nodes = 0;
  const long double width2 = std::sqrt(best2 / norm[2]);
  if (!(width2 < 1e15L)) {
    return -1.0L;
  }
  for (long long x2 = static_cast<long long>(std::ceil(coordinate[2] - width2));
       x2 <= static_cast<long long>(std::floor(coordinate[2] + width2)); ++x2) {
    const long double e2 = coordinate[2] - static_cast<long double>(x2);
    const long double r2 = e2 * e2 * norm[2];
    if (r2 > best2) {
      continue;
    }
    const long double c1 = coordinate[1] - static_cast<long double>(x2) * mu[2][1];
    const long double width1 = std::sqrt((best2 - r2) / norm[1]);
    if (!(width1 < 1e15L)) {
      return -1.0L;
    }
    for (long long x1 = static_cast<long long>(std::ceil(c1 - width1));
         x1 <= static_cast<long long>(std::floor(c1 + width1)); ++x1) {
      const long double e1 = c1 - static_cast<long double>(x1);
      const long double r1 = r2 + e1 * e1 * norm[1];
      if (r1 > best2) {
        continue;
      }
      const long double c0 = coordinate[0] - static_cast<long double>(x2) * mu[2][0] -
                             static_cast<long double>(x1) * mu[1][0];
      const long double width0 = std::sqrt((best2 - r1) / norm[0]);
      if (!(width0 < 1e15L)) {
        return -1.0L;
      }
      // One extra point either side absorbs rounding in the long double
      // Gram-Schmidt coordinates.
      for (long long x0 = static_cast<long long>(std::ceil(c0 - width0)) - 1;
           x0 <= static_cast<long long>(std::floor(c0 + width0)) + 1; ++x0) {
        if (++nodes > 50000000) {
          return -1.0L;
        }
        const long long n[3] = {x0 * u[0][0] + x1 * u[1][0] + x2 * u[2][0],
                                x0 * u[0][1] + x1 * u[1][1] + x2 * u[2][1],
                                x0 * u[0][2] + x1 * u[1][2] + x2 * u[2][2]};
        const long double l = lengthOf(latticeImage(c, w, n));
        if (best < 0.0L || l < best) {
          best = l;
        }
      }
    }
  }
  return best;
}

struct ImageError {
  double ulps = 0.0;
  double relative = 0.0;
};

// Error of the library image of the wrap `w` against the coefficient
// reference, in ulps of |w| and relative to the exact minimum.
ImageError imageError(const double c[6], const std::array<double, 3> &w) {
  constexpr long double eps = std::numeric_limits<double>::epsilon();
  const auto got = gen::mic::euclideanIfShorter(c[0], c[1], c[2], c[3], c[4], c[5], w, w);
  const long double found = lengthOf(got);
  const long double wrap = lengthOf(w);
  const double d[3] = {w[0], w[1], w[2]};
  const long double ref =
      coefficientImage(c, d, (found + 64.0L * eps * wrap + 1e-300L) * (1.0L + 1e-12L));
  REQUIRE(ref >= 0.0L);
  const long double err = std::abs(found - ref);
  return {static_cast<double>(err / (eps * wrap)),
          static_cast<double>(ref > 0.0L ? err / ref : err)};
}

TEST_CASE("minimum image error is within 4 ulps of the wrap and 1e-14 relative",
          "[pore]") {
  // The bounds stated in changelog.d/+pore-profiles.added.md, over the sets
  // they were measured on.
  long pairs = 0;
  const auto judge = [&](const double c[6], const std::array<double, 3> &w) {
    const ImageError e = imageError(c, w);
    ++pairs;
    REQUIRE(e.ulps <= 4.0);
    REQUIRE(e.relative <= 1e-14);
  };
  // An aspect-1e8 cell whose wrap is about 4e7 times its minimum image.
  const double gate[6] = {874115888.13459587, 26.140636200847869, 8.7411588813459584,
                          -55749478.818535045, 437057944.06729794, 4.6514428281782481};
  judge(gate, {437057944.06729794, 0.0, 0.0});
  const double elongated[6] = {1.0e7, 0.1, 0.1, 0.025, 0.025, 0.0};
  judge(elongated, {-4859743.2640920915, -0.0050198905816330854, -0.031675054668108597});
  judge(elongated, {2500000.0, 0.02, -0.03});

  const double halves[][3] = {{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5},
                              {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5},
                              {0.5, 0.5, 0.5}};
  const auto wrapOf = [](const double c[6], const double f[3]) {
    return std::array<double, 3>{c[0] * f[0] + c[3] * f[1] + c[4] * f[2],
                                 c[1] * f[1] + c[5] * f[2], c[2] * f[2]};
  };
  {
    // The aspect-40 generator of the triangular brute-force test above.
    std::mt19937_64 gen(20261120);
    std::uniform_real_distribution<double> logScale(-1.0, 3.0);
    std::uniform_real_distribution<double> logAspect(0.0, std::log10(40.0));
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> half(-0.5, 0.5);
    for (int cell = 0; cell < 1200; ++cell) {
      const double base = std::pow(10.0, logScale(gen));
      double edge[3] = {base, base * std::pow(10.0, logAspect(gen) * unit(gen)),
                        base * std::pow(10.0, logAspect(gen))};
      std::shuffle(edge, edge + 3, gen);
      double c[6] = {edge[0], edge[1], edge[2], half(gen) * edge[0], half(gen) * edge[0],
                     half(gen) * edge[1]};
      if (cell % 4 == 1) {
        c[3] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[0];
      }
      if (cell % 4 == 2) {
        c[3] = 0.5 * edge[0];
        c[4] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[0];
        c[5] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[1];
      }
      for (int p = 0; p < 24; ++p) {
        double f[3];
        for (int k = 0; k < 3; ++k) {
          f[k] = p < 7 ? halves[p][k] : half(gen);
        }
        judge(c, wrapOf(c, f));
      }
    }
  }
  // 3000 LAMMPS-legal cells per edge aspect: one edge is the base, one is
  // aspect times it, the third lies between. By index mod 4 the tilts are
  // random, one is at +-edge/2, all three are, or xz is half the long lx.
  unsigned long seed = 20261122;
  for (const double aspect : {1e3, 1e5, 1e8}) {
    std::mt19937_64 gen(seed++);
    std::uniform_real_distribution<double> logScale(-1.0, 2.0);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> half(-0.5, 0.5);
    const double la = std::log10(aspect);
    for (int cell = 0; cell < 3000; ++cell) {
      const double base = std::pow(10.0, logScale(gen));
      double edge[3] = {base, base * aspect, base * std::pow(10.0, la * unit(gen))};
      if (cell % 4 == 3) {
        std::swap(edge[0], edge[1]);
        if (unit(gen) < 0.5) {
          std::swap(edge[1], edge[2]);
        }
      } else {
        std::shuffle(edge, edge + 3, gen);
      }
      double c[6] = {edge[0], edge[1], edge[2], half(gen) * edge[0], half(gen) * edge[0],
                     half(gen) * edge[1]};
      if (cell % 4 == 1) {
        c[3 + (cell / 4) % 3] =
            (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[(cell / 4) % 3 == 2 ? 1 : 0];
      }
      if (cell % 4 == 2) {
        c[3] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[0];
        c[4] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[0];
        c[5] = (unit(gen) < 0.5 ? 0.5 : -0.5) * edge[1];
      }
      if (cell % 4 == 3) {
        c[4] = 0.5 * edge[0];
      }
      for (int p = 0; p < 16; ++p) {
        double f[3];
        for (int k = 0; k < 3; ++k) {
          f[k] = p < 7 ? halves[p][k] : half(gen);
        }
        judge(c, wrapOf(c, f));
      }
    }
  }
  REQUIRE(pairs == 3 + 1200 * 24 + 3 * 3000 * 16);
}

TEST_CASE("elongated minimum images below 2^52 are returned", "[pore]") {
  // The Babai coordinate on the short reduced vector and the large
  // coefficient on a different vector used to be multiplied together and
  // rejected. Each product that is actually formed stays near 1e8.
  const auto coefficient = [](const double c[6], const std::array<double, 3> &w,
                              const std::array<double, 3> &got) {
    const long double cz = (static_cast<long double>(w[2]) - got[2]) / c[2];
    const long double cy =
        ((static_cast<long double>(w[1]) - got[1]) - static_cast<long double>(c[5]) * cz) / c[1];
    const long double cx = ((static_cast<long double>(w[0]) - got[0]) -
                            static_cast<long double>(c[3]) * cy -
                            static_cast<long double>(c[4]) * cz) /
                           c[0];
    return std::array<long double, 3>{cx, cy, cz};
  };
  const auto expect = [&](const double c[6], const double f[3], long long nx, long long ny,
                          long long nz) {
    const std::array<double, 3> w = {c[0] * f[0] + c[3] * f[1] + c[4] * f[2],
                                     c[1] * f[1] + c[5] * f[2], c[2] * f[2]};
    const auto got = gen::mic::euclideanIfShorter(c[0], c[1], c[2], c[3], c[4], c[5], w, w);
    const auto n = coefficient(c, w, got);
    REQUIRE(std::llround(n[0]) == nx);
    REQUIRE(std::llround(n[1]) == ny);
    REQUIRE(std::llround(n[2]) == nz);
    REQUIRE(std::abs(n[0] - static_cast<long double>(nx)) < 1e-4L);
    REQUIRE(std::abs(n[1] - static_cast<long double>(ny)) < 1e-4L);
    REQUIRE(std::abs(n[2] - static_cast<long double>(nz)) < 1e-4L);
  };
  const double skewed[6] = {120000000, 120000000, 1, 60000000, 1, 1};
  const double skewedFrac[3] = {0.5, 0.5, -0.5};
  expect(skewed, skewedFrac, 1, 1, -50000000);
  const double needle[6] = {500000000, 1, 1, 0.5, 0, 0};
  const double needleFrac[3] = {0.5, 0, 0};
  expect(needle, needleFrac, 1, -100000000, 0);
}

TEST_CASE("minimum image is scale-invariant in LAMMPS cgs and si units", "[pore]") {
  const double base[6] = {30.0, 28.0, 26.0, 14.0, 11.0, 12.0};
  const std::array<double, 3> baseW = {25.0, 18.0, 12.0};
  const auto refMic = gen::mic::euclideanIfShorter(
      base[0], base[1], base[2], base[3], base[4], base[5], baseW, baseW);
  const double refNorm =
      std::sqrt(refMic[0] * refMic[0] + refMic[1] * refMic[1] + refMic[2] * refMic[2]);

  // Real/metal (A), nano (nm), micro (um), cgs (cm, 1e-8 A), si (m, 1e-10 A),
  // plus extreme binary64 scales.
  for (const double scale : {1.0, 1e-1, 1e-4, 1e-8, 1e-10, 1e-100, 1e+100}) {
    const double lx = base[0] * scale;
    const double ly = base[1] * scale;
    const double lz = base[2] * scale;
    const double xy = base[3] * scale;
    const double xz = base[4] * scale;
    const double yz = base[5] * scale;
    const std::array<double, 3> w = {baseW[0] * scale, baseW[1] * scale, baseW[2] * scale};
    const auto gotMic = gen::mic::euclideanIfShorter(lx, ly, lz, xy, xz, yz, w, w);
    for (int k = 0; k < 3; ++k) {
      REQUIRE_THAT(gotMic[k] / scale, Catch::Matchers::WithinAbs(refMic[k], 1e-12 * refNorm));
    }

    const auto cloud = triclinic(lx, ly, lz, xy, xz, yz);
    const auto gotPub = gen::triclinicMinImage(cloud, w[0], w[1], w[2], 0.0, 0.0, 0.0);
    for (int k = 0; k < 3; ++k) {
      REQUIRE_THAT(gotPub[k] / scale, Catch::Matchers::WithinAbs(refMic[k], 1e-12 * refNorm));
    }
  }
}

TEST_CASE("minimum image basis build is bounded on high-aspect legal needles", "[pore]") {
  const double ly = 1.6877791608068398;
  const double lz = 1.7939748045334107;
  const double yz = 0.12097869737185198;
  for (const double L : {1.12982562035e11, 1.0e15}) {
    for (const double xy : {0.0, 1.0e-12}) {
      const std::array<double, 3> w = {0.3 * L, 0.2 * ly, 0.4 * lz};
      const auto got = gen::mic::euclideanIfShorter(L, ly, lz, xy, 0.0, yz, w, w);
      REQUIRE_THAT(got[0], Catch::Matchers::WithinRel(0.3 * L, 1e-14));
      REQUIRE_THAT(got[1], Catch::Matchers::WithinAbs(0.2 * ly, 1e-12));
      REQUIRE_THAT(got[2], Catch::Matchers::WithinAbs(0.4 * lz, 1e-12));
    }
  }
}
