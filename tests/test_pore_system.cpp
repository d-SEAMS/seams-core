// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <seams.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

double volumeSum(const seams::FractionalProfile &profile) {
  double sum = 0.0;
  for (double volume : profile.volume) {
    sum += volume;
  }
  return sum;
}

std::array<double, 3> cartesian(const std::array<double, 9> &box, double sx, double sy, double sz) {
  const double lx = box[3];
  const double ly = box[4];
  const double lz = box[5];
  const double xy = box[6];
  const double xz = box[7];
  const double yz = box[8];
  return {box[0] + lx * sx + xy * sy + xz * sz, box[1] + ly * sy + yz * sz, box[2] + lz * sz};
}

double coveringRadius(double u0, double u1, double v0, double v1) {
  double b1[2] = {u0, u1};
  double b2[2] = {v0, v1};
  for (int step = 0; step < 32; ++step) {
    if (b1[0] * b1[0] + b1[1] * b1[1] > b2[0] * b2[0] + b2[1] * b2[1]) {
      std::swap(b1[0], b2[0]);
      std::swap(b1[1], b2[1]);
    }
    const double base = b1[0] * b1[0] + b1[1] * b1[1];
    const double m = std::round((b1[0] * b2[0] + b1[1] * b2[1]) / base);
    if (m == 0.0) {
      if (b1[0] * b2[0] + b1[1] * b2[1] < 0.0) {
        b2[0] = -b2[0];
        b2[1] = -b2[1];
      }
      const double cross = std::abs(b1[0] * b2[1] - b1[1] * b2[0]);
      return std::hypot(b1[0], b1[1]) * std::hypot(b2[0], b2[1]) *
             std::hypot(b1[0] - b2[0], b1[1] - b2[1]) / (2.0 * cross);
    }
    b2[0] -= m * b1[0];
    b2[1] -= m * b1[1];
  }
  return 0.0;
}

double poreCover(double lx, double ly, double lz, double xy, double xz, double yz) {
  const double c2 = xz * xz + yz * yz + lz * lz;
  const double along[3] = {xz, yz, lz};
  const double axial = std::sqrt(c2);
  const double eaxis[3] = {along[0] / axial, along[1] / axial, along[2] / axial};
  const double pu = (lx * along[0]) / c2;
  const double pv = (xy * along[0] + ly * along[1]) / c2;
  const double raw_u[3] = {lx - pu * along[0], -pu * along[1], -pu * along[2]};
  const double raw_v[3] = {xy - pv * along[0], ly - pv * along[1], -pv * along[2]};
  const double nu = std::sqrt(raw_u[0] * raw_u[0] + raw_u[1] * raw_u[1] + raw_u[2] * raw_u[2]);
  const double e1[3] = {raw_u[0] / nu, raw_u[1] / nu, raw_u[2] / nu};
  double e2[3] = {eaxis[1] * e1[2] - eaxis[2] * e1[1], eaxis[2] * e1[0] - eaxis[0] * e1[2],
                  eaxis[0] * e1[1] - eaxis[1] * e1[0]};
  const double n2 = std::sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
  const double ev0 = raw_v[0] * e1[0] + raw_v[1] * e1[1] + raw_v[2] * e1[2];
  const double ev1 = (raw_v[0] * e2[0] + raw_v[1] * e2[1] + raw_v[2] * e2[2]) / n2;
  return coveringRadius(nu, 0.0, ev0, ev1);
}

}  // namespace

TEST_CASE("System fractional density keeps shifted bins on the cell volume", "[pore][system]") {
  seams::System sys;
  const std::vector<double> xyz = {5.0, 5.0, 5.0};
  sys.set_ortho(xyz, 10.0, 10.0, 10.0);
  seams::FractionalDensityOptions options;
  options.bins = 4;
  options.use_direction = true;
  options.direction = {0.0, 0.0, 1.0};
  options.origin = 0.6;
  const auto profile = sys.fractional_density(options);
  REQUIRE(profile.volume.size() == 4);
  for (double volume : profile.volume) {
    REQUIRE_THAT(volume, Catch::Matchers::WithinAbs(250.0, 1e-9));
  }
  REQUIRE_THAT(volumeSum(profile), Catch::Matchers::WithinAbs(1000.0, 1e-9));
}

TEST_CASE("System pore profiles give half-cell azimuthal bins", "[pore][system]") {
  seams::System sys;
  std::vector<double> xyz;
  std::vector<int> types;
  const int grid = 6;
  for (int ix = 0; ix < grid; ++ix) {
    for (int iy = 0; iy < grid; ++iy) {
      for (int iz = 0; iz < grid; ++iz) {
        xyz.push_back((static_cast<double>(ix) + 0.5) / grid * 10.0);
        xyz.push_back((static_cast<double>(iy) + 0.5) / grid * 10.0);
        xyz.push_back((static_cast<double>(iz) + 0.5) / grid * 10.0);
        types.push_back(1);
      }
    }
  }
  sys.set_ortho(xyz, 10.0, 10.0, 10.0);
  seams::PoreProfileOptions options;
  options.radial_bins = 4;
  options.azimuthal_bins = 2;
  options.axial_bins = 4;
  options.axis = 2;
  options.use_centre = true;
  options.centre_fractional = {0.5, 0.5, 0.5};
  const auto profiles = sys.pore_profiles(options);
  const double density = static_cast<double>(grid * grid * grid) / 1000.0;
  REQUIRE_THAT(profiles.azimuthal.volume[0], Catch::Matchers::WithinRel(500.0, 1e-8));
  REQUIRE_THAT(profiles.azimuthal.volume[1], Catch::Matchers::WithinRel(500.0, 1e-8));
  for (double rho : profiles.azimuthal.number_density) {
    REQUIRE_THAT(rho, Catch::Matchers::WithinRel(density, 1e-8));
  }
  REQUIRE_THAT(volumeSum(profiles.axial), Catch::Matchers::WithinRel(1000.0, 1e-12));
  REQUIRE_THAT(volumeSum(profiles.radial), Catch::Matchers::WithinRel(1000.0, 1e-8));

  seams::PoreFrameOptions frame_options;
  frame_options.axis = 2;
  frame_options.use_centre = true;
  frame_options.centre_fractional = {0.5, 0.5, 0.5};
  const auto frame = sys.pore_frame(frame_options);
  REQUIRE(frame.axis_index == 2);
  REQUIRE_THAT(frame.axis[2], Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(frame.centre_fractional[0], Catch::Matchers::WithinAbs(0.5, 1e-12));
  REQUIRE_THAT(frame.centre_fractional[1], Catch::Matchers::WithinAbs(0.5, 1e-12));
}

TEST_CASE("System host guest uses the triclinic minimum image", "[pore][system]") {
  const std::array<double, 9> box = {1.0, -2.0, 0.5, 10.0, 12.0, 14.0, 4.5, -3.0, 5.5};
  const double cutoff = 6.0;
  std::vector<double> xyz;
  std::vector<int> types;
  std::vector<std::array<double, 3>> guests;
  std::vector<std::array<double, 3>> frames;
  auto push = [&](int type, double sx, double sy, double sz) {
    const auto r = cartesian(box, sx, sy, sz);
    xyz.insert(xyz.end(), r.begin(), r.end());
    types.push_back(type);
  };
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 2; ++k) {
        const double s[3] = {(static_cast<double>(i) + 0.17) / 3.0,
                             (static_cast<double>(j) + 0.41) / 3.0,
                             (static_cast<double>(k) + 0.63) / 2.0};
        push(1, s[0], s[1], s[2]);
        guests.push_back({s[0], s[1], s[2]});
      }
    }
  }
  const int nguest = static_cast<int>(guests.size());
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double s[3] = {(static_cast<double>(i) + 0.5) / 3.0, (static_cast<double>(j) + 0.2) / 3.0,
                           0.8};
      push(2, s[0], s[1], s[2]);
      frames.push_back({s[0], s[1], s[2]});
    }
  }
  const double lx = box[3];
  const double ly = box[4];
  const double lz = box[5];
  const double xy = box[6];
  const double xz = box[7];
  const double yz = box[8];
  int wrap_disagrees = 0;
  std::vector<double> expect_distance(static_cast<std::size_t>(nguest));
  std::vector<int> expect_contacts(static_cast<std::size_t>(nguest));
  for (int g = 0; g < nguest; ++g) {
    double best = 1e300;
    int contacts = 0;
    for (const auto &site : frames) {
      double ds[3] = {guests[static_cast<std::size_t>(g)][0] - site[0],
                      guests[static_cast<std::size_t>(g)][1] - site[1],
                      guests[static_cast<std::size_t>(g)][2] - site[2]};
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
        wrap_disagrees += 1;
      }
      best = std::min(best, image);
      if (image <= cutoff) {
        contacts += 1;
      }
    }
    expect_distance[static_cast<std::size_t>(g)] = best;
    expect_contacts[static_cast<std::size_t>(g)] = contacts;
  }
  REQUIRE(wrap_disagrees > 0);
  seams::System sys;
  sys.set_frame(xyz, box, types);
  seams::HostGuestOptions options;
  options.guest = 1;
  options.framework = 2;
  options.cutoff = cutoff;
  const auto report = sys.host_guest(options);
  REQUIRE(report.guest_index.size() == static_cast<std::size_t>(nguest));
  REQUIRE(report.framework_types.size() == 1);
  REQUIRE(report.framework_types[0] == 2);
  for (int g = 0; g < nguest; ++g) {
    REQUIRE_THAT(report.nearest_distance[static_cast<std::size_t>(g)],
                 Catch::Matchers::WithinAbs(expect_distance[static_cast<std::size_t>(g)], 1e-8));
    REQUIRE(report.contacts[static_cast<std::size_t>(g)] == expect_contacts[static_cast<std::size_t>(g)]);
  }
}

TEST_CASE("seams_pore_profiles at the covering radius returns the cell volume", "[pore][system]") {
  const double named[][6] = {
      {20.0, 10.0 * std::sqrt(3.0), 8.0, 10.0, 0.0, 0.0},
      {10.0, 12.0, 14.0, 4.5, -3.0, 5.5},
      {10.0, 8.0, 12.0, 2.0, 1.0, -1.5},
  };
  auto expect = [](double lx, double ly, double lz, double xy, double xz, double yz) {
    const double cover = poreCover(lx, ly, lz, xy, xz, yz);
    const double volume = lx * ly * lz;
    REQUIRE(cover > 0.0);
    const std::array<double, 9> box = {0.0, 0.0, 0.0, lx, ly, lz, xy, xz, yz};
    const std::vector<double> xyz = {0.1, 0.2, 0.3};
    const std::vector<int> types = {1};
    seams::System sys;
    sys.set_frame(xyz, box, types);
    seams::PoreProfileOptions options;
    options.axis = 2;
    options.use_centre = true;
    options.centre_fractional = {0.5, 0.5, 0.5};
    const double scales[] = {1.0 - 1e-15, 1.0 - 1e-13, 1.0 - 1e-12, 1.0};
    for (double scale : scales) {
      options.r_max = cover * scale;
      const auto profiles = sys.pore_profiles(options);
      // Caps beyond the circle are about 1e-11 of the cell at 1e-12, and
      // far smaller at one ULP. The disk pi r^2 L is 1.209 V on the hex
      // cell. The entry hit in place of the exit hit is several percent off
      // at R and within a few ULPs of R.
      REQUIRE_THAT(profiles.radial.volume[0], Catch::Matchers::WithinRel(volume, 1e-8));
    }
    double radius = cover;
    for (int step = 0; step < 6; ++step) {
      radius = std::nextafter(radius, 0.0);
    }
    for (int step = -6; step <= 6; ++step) {
      options.r_max = radius;
      const auto profiles = sys.pore_profiles(options);
      REQUIRE_THAT(profiles.radial.volume[0], Catch::Matchers::WithinRel(volume, 1e-8));
      radius = std::nextafter(radius, 1.0e300);
    }
  };
  for (const auto &cell : named) {
    expect(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
  }
  std::mt19937_64 gen(20261120);
  std::uniform_real_distribution<double> length(6.0, 16.0);
  std::uniform_real_distribution<double> tilt(-0.9, 0.9);
  for (int trial = 0; trial < 8; ++trial) {
    const double lx = length(gen);
    const double ly = length(gen);
    const double lz = length(gen);
    expect(lx, ly, lz, tilt(gen) * lx, tilt(gen) * lx, tilt(gen) * ly);
  }
}

TEST_CASE("System rejects a negative axis other than -1 and a repeated type", "[pore][system]") {
  seams::System sys;
  const std::vector<double> xyz = {1.0, 1.0, 1.0, 4.0, 1.0, 1.0};
  const std::vector<int> types = {2, 5};
  const std::array<double, 9> box = {0.0, 0.0, 0.0, 10.0, 9.0, 8.0, 0.0, 0.0, 0.0};
  sys.set_frame(xyz, box, types);
  seams::PoreProfileOptions profiles;
  profiles.axis = -7;
  profiles.r_max = 1.0;
  try {
    sys.pore_profiles(profiles);
    FAIL("axis -7 was accepted");
  } catch (const seams::Error &error) {
    REQUIRE(error.code() == SEAMS_ERR_RANGE);
  }
  profiles.axis = -1;
  profiles.r_max = 0.0;
  const auto made = sys.pore_profiles(profiles);
  REQUIRE(made.radial.volume.size() == 1);
  REQUIRE_THAT(made.radial.volume[0], Catch::Matchers::WithinRel(10.0 * 9.0 * 8.0, 1e-8));
  seams::HostGuestOptions same;
  same.guest = 2;
  same.framework = 2;
  same.cutoff = 3.0;
  try {
    sys.host_guest(same);
    FAIL("a repeated guest and framework type was accepted");
  } catch (const seams::Error &error) {
    REQUIRE(error.code() == SEAMS_ERR_ARG);
  }
  same.framework = 0;
  const auto report = sys.host_guest(same);
  REQUIRE(report.guest_index.size() == 1);
  REQUIRE(report.framework_types.size() == 1);
  REQUIRE(report.framework_types[0] == 5);
}

TEST_CASE("seams_host_guest on a long tilted cell matches brute force", "[pore][system]") {
  const double cells[][4] = {
      {12.0, 12.0, 100.0, 0.5}, {12.0, 12.0, 400.0, 0.5}, {12.0, 12.0, 1000.0, 0.5},
      {10.0, 10.0, 2000.0, 1.0e-6},
  };
  for (const auto &cell : cells) {
    const double lx = cell[0];
    const double ly = cell[1];
    const double lz = cell[2];
    const double xy = cell[3];
    const std::vector<double> xyz = {1.0, 1.0, 1.0, 5.0, 6.0, 0.45 * lz};
    const std::vector<int> types = {1, 2};
    const std::array<double, 9> box = {0.0, 0.0, 0.0, lx, ly, lz, xy, 0.0, 0.0};
    seams::System sys;
    sys.set_frame(xyz, box, types);
    seams::HostGuestOptions options;
    options.guest = 1;
    options.framework = 2;
    options.cutoff = 3.0;
    const auto report = sys.host_guest(options);
    double best = 1e300;
    const double d[3] = {xyz[3] - xyz[0], xyz[4] - xyz[1], xyz[5] - xyz[2]};
    for (int i = -3; i <= 3; ++i) {
      for (int j = -3; j <= 3; ++j) {
        for (int k = -3; k <= 3; ++k) {
          const double x = d[0] + i * lx + j * xy;
          const double y = d[1] + j * ly;
          const double z = d[2] + k * lz;
          best = std::min(best, std::sqrt(x * x + y * y + z * z));
        }
      }
    }
    REQUIRE(report.nearest_distance.size() == 1);
    REQUIRE_THAT(report.nearest_distance[0], Catch::Matchers::WithinAbs(best, 1e-8));
  }
}

TEST_CASE("seams_host_guest refines minimage on an x-elongated cell",
          "[pore][system]") {
  constexpr double lx = 10000.0;
  constexpr double ly = 12.0;
  constexpr double lz = 12.0;
  constexpr double xy = 3.0;
  constexpr double xz = 3.0;
  const std::vector<double> xyz = {5000.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const std::vector<int> types = {1, 2};
  const std::array<double, 9> box = {0.0, 0.0, 0.0, lx, ly, lz, xy, xz, 0.0};
  seams::System sys;
  sys.set_frame(xyz, box, types);
  seams::HostGuestOptions options;
  options.guest = 1;
  options.framework = 2;
  options.cutoff = 3.0;
  const auto report = sys.host_guest(options);

  double best = std::numeric_limits<double>::infinity();
  for (int j = -500; j <= 500; ++j) {
    for (int k = -500; k <= 500; ++k) {
      const double shiftedX = xyz[0] + j * xy + k * xz;
      const int i = static_cast<int>(std::round(-shiftedX / lx));
      const double x = shiftedX + i * lx;
      const double y = j * ly;
      const double z = k * lz;
      best = std::min(best, std::sqrt(x * x + y * y + z * z));
    }
  }
  REQUIRE(report.nearest_distance.size() == 1);
  REQUIRE_THAT(report.nearest_distance[0],
               Catch::Matchers::WithinAbs(best, 1e-8));
}

namespace {

// Provably exact minimum image for upper-triangular H (lx, ly, lz, xy, xz, yz):
// z bounds k, y bounds j, and for fixed (j, k) rounding selects optimal i.
double bruteTriclinicDist(double lx, double ly, double lz, double xy, double xz,
                          double yz, double dx, double dy, double dz) {
  const double kw = std::round(-dz / lz);
  const double zw = dz + kw * lz;
  const double y0w = dy + kw * yz;
  const double jw = std::round(-y0w / ly);
  const double yw = y0w + jw * ly;
  const double x0w = dx + jw * xy + kw * xz;
  const double xw = x0w + std::round(-x0w / lx) * lx;
  double best = xw * xw + yw * yw + zw * zw;
  const double R = std::sqrt(best);
  const long k0 = static_cast<long>(std::floor((-dz - R) / lz)) - 1;
  const long k1 = static_cast<long>(std::ceil((-dz + R) / lz)) + 1;
  for (long k = k0; k <= k1; ++k) {
    const double z = dz + static_cast<double>(k) * lz;
    if (z * z > best) {
      continue;
    }
    const double rem = std::sqrt(best - z * z);
    const double y0 = dy + static_cast<double>(k) * yz;
    const long j0 = static_cast<long>(std::floor((-y0 - rem) / ly)) - 1;
    const long j1 = static_cast<long>(std::ceil((-y0 + rem) / ly)) + 1;
    for (long j = j0; j <= j1; ++j) {
      const double y = y0 + static_cast<double>(j) * ly;
      const double x0 = dx + static_cast<double>(j) * xy + static_cast<double>(k) * xz;
      const double ic = std::round(-x0 / lx);
      for (double i = ic - 1.0; i <= ic + 1.0; i += 1.0) {
        const double x = x0 + i * lx;
        best = std::min(best, x * x + y * y + z * z);
      }
    }
  }
  return std::sqrt(best);
}

}  // namespace

TEST_CASE("System host_guest does not poison the basis cache after an error",
          "[pore][system]") {
  const double fillers[][6] = {
      {11.0, 10.0, 10.0, 1.0, 0.0, 0.0},
      {12.0, 10.0, 10.0, 1.0, 0.0, 0.0},
      {13.0, 10.0, 10.0, 1.0, 0.0, 0.0},
      {14.0, 10.0, 10.0, 1.0, 0.0, 0.0},
  };
  for (const auto &c : fillers) {
    seams::System sys;
    const std::array<double, 9> box = {0.0, 0.0, 0.0, c[0], c[1], c[2], c[3], c[4], c[5]};
    const std::vector<double> xyz = {5.5, 5.0, 5.0, 0.0, 0.0, 0.0};
    const std::vector<int> types = {1, 2};
    sys.set_frame(xyz, box, types);
    seams::HostGuestOptions opt;
    opt.guest = 1;
    opt.framework = 2;
    opt.cutoff = 2.0;
    const auto rep = sys.host_guest(opt);
    REQUIRE(rep.nearest_distance.size() == 1);
  }

  // xy = 5e15 with lx = 1 preserves lx = 1 exactly in binary64 bound spans
  // while forcing Lagrange reduction to throw on coefficient >= 2^52.
  seams::System badSys;
  const std::array<double, 9> badBox = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5.0e15, 0.0, 0.0};
  const std::vector<double> badXyz = {0.5, 0.5, 0.5, 0.0, 0.0, 0.0};
  const std::vector<int> badTypes = {1, 2};
  badSys.set_frame(badXyz, badBox, badTypes);
  seams::HostGuestOptions opt;
  opt.guest = 1;
  opt.framework = 2;
  opt.cutoff = 1.0;
  REQUIRE_THROWS_AS(badSys.host_guest(opt), seams::Error);
  REQUIRE_THROWS_AS(badSys.host_guest(opt), seams::Error);
}

TEST_CASE("System host_guest and seams_rdf succeed on ordinary triclinic cells",
          "[pore][system]") {
  const double s3 = std::sqrt(3.0);
  struct NamedCell {
    double lx, ly, lz, xy, xz, yz;
    double rmax;
  };
  const NamedCell cells[] = {
      {20.0, 20.0 * s3 / 2.0, 20.0 * std::sqrt(2.0 / 3.0), 10.0, 10.0,
       20.0 / (2.0 * s3), 12.0},
      {14.142135623730951, 12.24744871391589, 11.547005383792516,
       7.0710678118654755, 7.0710678118654755, 4.08248290463863, 8.0},
      {20.0, 20.0, 20.0, 1.0e-16, 0.0, 0.0, 12.0},
      {20.0, 20.0, 20.0, 1.0e-12, 0.0, 0.0, 12.0},
      {20.0, 20.0, 20.0, 1.0e-12, 1.0e-12, 1.0e-12, 12.0},
      {1000.0, 12.0, 20.4, 500.0, -350.0, 0.0, 12.0},
      {25.0, 21.650635094610966, 30.0, 12.5, 3.3, 1.1, 12.0},
      {58.196444, 53.3197415, 54.9325006, 29.098222, -29.098222, -26.6598708, 40.0},
      {14.745736, 71.0592255, 3.85588635, 4.0013901, 3.48405873, 32.153036, 20.0},
  };
  std::mt19937_64 rng(20261121);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (const auto &c : cells) {
    const std::array<double, 9> box = {0.0, 0.0, 0.0, c.lx, c.ly, c.lz, c.xy, c.xz, c.yz};
    // 1) host_guest against exact brute force.
    {
      const int nguest = 40;
      std::vector<double> xyz;
      std::vector<int> types;
      for (int i = 0; i < nguest; ++i) {
        const double f[3] = {u(rng), u(rng), u(rng)};
        xyz.insert(xyz.end(), {c.lx * f[0] + c.xy * f[1] + c.xz * f[2],
                               c.ly * f[1] + c.yz * f[2], c.lz * f[2]});
        types.push_back(1);
      }
      xyz.insert(xyz.end(), {0.0, 0.0, 0.0});
      types.push_back(2);
      seams::System sys;
      sys.set_frame(xyz, box, types);
      seams::HostGuestOptions opt;
      opt.guest = 1;
      opt.framework = 2;
      opt.cutoff = 1.0;
      const auto rep = sys.host_guest(opt);
      REQUIRE(rep.nearest_distance.size() == static_cast<std::size_t>(nguest));
      for (int i = 0; i < nguest; ++i) {
        const double ref = bruteTriclinicDist(
            c.lx, c.ly, c.lz, c.xy, c.xz, c.yz, xyz[3 * static_cast<std::size_t>(i)],
            xyz[3 * static_cast<std::size_t>(i) + 1],
            xyz[3 * static_cast<std::size_t>(i) + 2]);
        REQUIRE_THAT(rep.nearest_distance[static_cast<std::size_t>(i)],
                     Catch::Matchers::WithinAbs(ref, 1e-9));
      }
    }
    // 2) seams_rdf succeeds past the Smith ball (rmax = c.rmax, exercising
    //    periodicDistSq -> euclideanIfShorter) and matches exact brute force
    //    within the single-image ball.
    {
      const int natoms = 60;
      std::vector<double> xyz;
      std::vector<int> types;
      for (int i = 0; i < natoms; ++i) {
        const double f[3] = {u(rng), u(rng), u(rng)};
        xyz.insert(xyz.end(), {c.lx * f[0] + c.xy * f[1] + c.xz * f[2],
                               c.ly * f[1] + c.yz * f[2], c.lz * f[2]});
        types.push_back(1);
      }
      seams::System sys;
      sys.set_frame(xyz, box, types);
      seams::DistributionOptions opt;
      opt.rmax = c.rmax;
      opt.bins = 10;
      const auto dist = sys.rdf(opt);
      long gotPairs = 0;
      for (int count : dist.count) {
        gotPairs += count;
      }
      REQUIRE(gotPairs > 0);

      // Single-image ball check against brute force.
      const double bcz = c.xy * c.yz - c.ly * c.xz;
      const double wx = c.lx * c.ly * c.lz /
                        std::sqrt(c.ly * c.lz * c.ly * c.lz +
                                  c.xy * c.lz * c.xy * c.lz + bcz * bcz);
      const double wy = c.ly * c.lz / std::sqrt(c.lz * c.lz + c.yz * c.yz);
      const double wz = c.lz;
      const double innerRmax = 0.45 * std::min({wx, wy, wz});
      seams::DistributionOptions innerOpt;
      innerOpt.rmax = innerRmax;
      innerOpt.bins = 10;
      const auto innerDist = sys.rdf(innerOpt);
      long innerGot = 0;
      for (int count : innerDist.count) {
        innerGot += count;
      }
      long innerRef = 0;
      for (int i = 0; i < natoms; ++i) {
        for (int j = i + 1; j < natoms; ++j) {
          const double r = bruteTriclinicDist(
              c.lx, c.ly, c.lz, c.xy, c.xz, c.yz,
              xyz[3 * static_cast<std::size_t>(j)] - xyz[3 * static_cast<std::size_t>(i)],
              xyz[3 * static_cast<std::size_t>(j) + 1] -
                  xyz[3 * static_cast<std::size_t>(i) + 1],
              xyz[3 * static_cast<std::size_t>(j) + 2] -
                  xyz[3 * static_cast<std::size_t>(i) + 2]);
          if (r <= innerRmax) {
            innerRef += 1;
          }
        }
      }
      REQUIRE(innerGot == innerRef);
    }
  }
}

namespace {

uint64_t g_state = 0;

double rnd() {
  g_state ^= g_state << 13;
  g_state ^= g_state >> 7;
  g_state ^= g_state << 17;
  return static_cast<double>(g_state >> 11) / 9007199254740992.0;
}

void seed_rng(unsigned long long seed) { g_state = seed * 2654435761ull + 88172645463325252ull; }

void fill_cloud(std::vector<double> &xyz, std::vector<int> &types, int n, double ly, double lz,
                double xy, double xz, double yz, unsigned long long seed, double x_scale) {
  seed_rng(seed);
  xyz.assign(static_cast<std::size_t>(3 * n), 0.0);
  types.assign(static_cast<std::size_t>(n), 0);
  for (int i = 0; i < n; ++i) {
    const double fb = rnd();
    const double fc = rnd();
    const double xa = x_scale * rnd();
    xyz[static_cast<std::size_t>(3 * i)] = xa + fb * xy + fc * xz;
    xyz[static_cast<std::size_t>(3 * i + 1)] = fb * ly + fc * yz;
    xyz[static_cast<std::size_t>(3 * i + 2)] = fc * lz;
    types[static_cast<std::size_t>(i)] = 1 + (i % 2);
  }
}

} // namespace

TEST_CASE("a long ordinary cell with atoms far from the origin is not a range error",
          "[pore][system]") {
  std::vector<double> xyz;
  std::vector<int> types;
  fill_cloud(xyz, types, 100, 30.0, 28.0, 3.0, 1.5, 2.0, 7, 3.0);
  for (int i = 0; i < 100; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] += 9000.0;
  }
  const std::array<double, 9> box = {0.0, 0.0, 0.0, 1.0e4, 30.0, 28.0, 3.0, 1.5, 2.0};
  REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
  seams::System sys;
  sys.set_frame(xyz, box, types);
  seams::HostGuestOptions options;
  options.guest = 1;
  options.framework = 2;
  options.cutoff = 5.0;
  const auto report = sys.host_guest(options);
  REQUIRE(report.nearest_distance.size() > 0);
  for (double distance : report.nearest_distance) {
    REQUIRE(std::isfinite(distance));
    if (distance >= 0.0) {
      REQUIRE(distance < 40.0);
    }
  }
  const auto contacts = sys.contact_pairs("1=cationHead,2=anion");
  REQUIRE(contacts.pairs > 0);
}

TEST_CASE("contacts at the far end of a needle are exact or a range error", "[pore][system]") {
  const double edge = 1.0e17;
  const std::vector<double> xyz = {0.6665140243063847, 0.5750009372335341, 0.14932477255801682,
                                   1.7702453113768646, 0.8684923662106784, 0.0025740216324493304,
                                   edge,               0.9225423641868384, 0.03232585887320041,
                                   edge,               0.6289187723290619, 0.016358306193924622};
  const std::vector<int> types = {1, 2, 1, 2};
  const std::array<double, 9> box = {0.0, 0.0, 0.0, edge, 1.0, 0.4, 0.3, 0.3, 0.3};
  seams::System sys;
  sys.set_frame(xyz, box, types);
  for (const int threads : {1, 8}) {
    REQUIRE(seams_set_num_threads(threads) == SEAMS_OK);
    try {
      const auto contacts = sys.contact_pairs("1=cationHead,2=anion");
      REQUIRE(contacts.pairs == 1);
    } catch (const seams::Error &error) {
      REQUIRE(error.code() == SEAMS_ERR_RANGE);
    }
  }
}

TEST_CASE("a cell past binary64 distance precision is a typed range error", "[pore][system]") {
  std::vector<double> ends_xyz;
  std::vector<int> ends_types;
  fill_cloud(ends_xyz, ends_types, 4, 1.0, 0.4, 0.3, 0.3, 0.3, 9, 3.0);
  ends_xyz[6] = 1.0e200;
  ends_xyz[9] = 1.0e200;
  const std::array<double, 9> ends_box = {0.0, 0.0, 0.0, 1.0e200, 1.0, 0.4, 0.3, 0.3, 0.3};
  const std::vector<double> cube_xyz = {0.0, 0.0, 0.0, 24.0e154, 0.0, 0.0};
  const std::vector<int> cube_types = {1, 2};
  const std::array<double, 9> cube_box = {0.0, 0.0, 0.0, 24.0e154, 24.0e154, 24.0e154, 0.0, 0.0, 0.0};
  const std::vector<double> *xyz[2] = {&ends_xyz, &cube_xyz};
  const std::vector<int> *types[2] = {&ends_types, &cube_types};
  const std::array<double, 9> *boxes[2] = {&ends_box, &cube_box};
  for (int k = 0; k < 2; ++k) {
    seams::System sys;
    sys.set_frame(*xyz[k], *boxes[k], *types[k]);
    bool ranged = false;
    try {
      sys.contact_pairs("1=cationHead,2=anion");
    } catch (const seams::Error &error) {
      REQUIRE(error.code() == SEAMS_ERR_RANGE);
      ranged = true;
    }
    REQUIRE(ranged);
    ranged = false;
    seams::HostGuestOptions options;
    options.guest = 1;
    options.framework = 2;
    options.cutoff = 5.0;
    try {
      sys.host_guest(options);
    } catch (const seams::Error &error) {
      REQUIRE(error.code() == SEAMS_ERR_RANGE);
      ranged = true;
    }
    REQUIRE(ranged);
  }
}

TEST_CASE("the contact-pair contract names the range error and the coefficient limit",
          "[pore][system]") {
  const std::string here = __FILE__;
  const auto tests = here.rfind("/tests/");
  REQUIRE(tests != std::string::npos);
  const std::string root = here.substr(0, tests);
  std::ifstream header(root + "/include/seams.h");
  REQUIRE(header.good());
  const std::string text((std::istreambuf_iterator<char>(header)), {});
  const auto fn = text.find("seams_status seams_contact_pairs(");
  REQUIRE(fn != std::string::npos);
  const auto doc = text.rfind("/**", fn);
  REQUIRE(doc != std::string::npos);
  const std::string block = text.substr(doc, fn - doc);
  REQUIRE(block.find("SEAMS_ERR_RANGE") != std::string::npos);
  REQUIRE(block.find("2^52") != std::string::npos);
  std::ifstream page(root + "/docs/api/c.md");
  REQUIRE(page.good());
  const std::string markdown((std::istreambuf_iterator<char>(page)), {});
  const auto section = markdown.rfind("### seams_contact_pairs");
  REQUIRE(section != std::string::npos);
  const auto next = markdown.find("\n### ", section + 4);
  const std::string contact = markdown.substr(section, next == std::string::npos ? std::string::npos : next - section);
  REQUIRE(contact.find("SEAMS_ERR_RANGE") != std::string::npos);
  REQUIRE(contact.find("2^52") != std::string::npos);
}

namespace {

struct Needle {
  std::vector<double> xyz;
  std::vector<int> types;
};

// The gate cloud: mt19937_64 seed 7, x = i*0.7 + 0.3*U, and, on the seam,
// the first tenth of the atoms placed just inside x = L.
Needle make_needle(double edge, int n, bool seam) {
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  Needle cloud;
  cloud.xyz.assign(static_cast<std::size_t>(3 * n), 0.0);
  cloud.types.assign(static_cast<std::size_t>(n), 0);
  const double ly = 1.0;
  const double lz = 0.4;
  for (int i = 0; i < n; ++i) {
    double x = static_cast<double>(i) * 0.7 + 0.3 * unit(rng);
    if (seam && i < n / 10) {
      x = edge - (static_cast<double>(i) * 0.7 + 0.5 * unit(rng)) - 1.0;
    }
    cloud.xyz[static_cast<std::size_t>(3 * i)] = x;
    cloud.xyz[static_cast<std::size_t>(3 * i) + 1] = ly * unit(rng);
    cloud.xyz[static_cast<std::size_t>(3 * i) + 2] = lz * unit(rng);
    cloud.types[static_cast<std::size_t>(i)] = 1 + (i % 3);
  }
  return cloud;
}

std::array<double, 9> needle_box(double edge) {
  return {0.0, 0.0, 0.0, edge, 1.0, 0.4, 0.3, 0.3, 0.3};
}

#if defined(__SIZEOF_FLOAT128__)
using Wide = __float128;
#else
using Wide = long double;
#endif

Wide wide_abs(Wide value) { return value < 0 ? -value : value; }

// Brute-force minimum image in a wide type. A pair whose squared length
// sits on the cutoff is left out of the comparison: the cutoff test is
// strict, and a tie there is not a wrong neighbour.
bool wide_neighbour(const std::vector<double> &xyz, int i, int j, double edge, double cutoff,
                    bool &on_cutoff) {
  const Wide ly = 1;
  const Wide lz = 0.4;
  const Wide xy = 0.3;
  const Wide xz = 0.3;
  const Wide yz = 0.3;
  const Wide dx = static_cast<Wide>(xyz[static_cast<std::size_t>(3 * i)]) -
                  static_cast<Wide>(xyz[static_cast<std::size_t>(3 * j)]);
  const Wide dy = static_cast<Wide>(xyz[static_cast<std::size_t>(3 * i) + 1]) -
                  static_cast<Wide>(xyz[static_cast<std::size_t>(3 * j) + 1]);
  const Wide dz = static_cast<Wide>(xyz[static_cast<std::size_t>(3 * i) + 2]) -
                  static_cast<Wide>(xyz[static_cast<std::size_t>(3 * j) + 2]);
  Wide best = static_cast<Wide>(1.0e300) * static_cast<Wide>(1.0e300);
  for (int nx = -1; nx <= 1; ++nx) {
    const Wide bx = dx + static_cast<Wide>(nx) * static_cast<Wide>(edge);
    if (wide_abs(bx) > 40) {
      continue;
    }
    for (int nz = -12; nz <= 12; ++nz) {
      for (int ny = -12; ny <= 12; ++ny) {
        const Wide x = bx + static_cast<Wide>(ny) * xy + static_cast<Wide>(nz) * xz;
        const Wide y = dy + static_cast<Wide>(ny) * ly + static_cast<Wide>(nz) * yz;
        const Wide z = dz + static_cast<Wide>(nz) * lz;
        const Wide d2 = x * x + y * y + z * z;
        if (d2 < best) {
          best = d2;
        }
      }
    }
  }
  const Wide limit = static_cast<Wide>(cutoff) * static_cast<Wide>(cutoff);
  const Wide gap = best - limit;
  on_cutoff = wide_abs(gap) < static_cast<Wide>(1.0e-8);
  return best < limit;
}

} // namespace

TEST_CASE("seam neighbours on a long needle match the wide minimum image", "[pore][system]") {
  constexpr int n = 260;
  constexpr double edge = 5.0e15;
  constexpr double cutoff = 3.5;
  const Needle cloud = make_needle(edge, n, true);
  const auto box = needle_box(edge);
  REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
  seams::System sys;
  sys.set_frame(cloud.xyz, box, cloud.types);
  std::vector<std::pair<int, int>> listed;
  try {
    listed = sys.neighbours(seams::NeighbourOptions{cutoff, 0, seams::Graph::Cutoff}).pairs();
  } catch (const seams::Error &error) {
    // The reduced basis of this needle needs a coefficient past 2^52.
    // A typed range error is the contract; a successful wrong list is not.
    REQUIRE(error.code() == SEAMS_ERR_RANGE);
    return;
  }
  std::set<std::pair<int, int>> got(listed.begin(), listed.end());
  long miss = 0;
  long extra = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        continue;
      }
      bool on_cutoff = false;
      const bool want = wide_neighbour(cloud.xyz, i, j, edge, cutoff, on_cutoff);
      if (on_cutoff) {
        continue;
      }
      const bool have = got.count({i, j}) != 0;
      if (want && !have) {
        ++miss;
      }
      if (!want && have) {
        ++extra;
      }
    }
  }
  REQUIRE(miss == 0);
  REQUIRE(extra == 0);
}

TEST_CASE("q6 on a long needle matches the compact cell", "[pore][system]") {
  constexpr int n = 260;
  constexpr double cutoff = 3.5;
  const Needle compact = make_needle(1.0e3, n, false);
  REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
  seams::System reference;
  reference.set_frame(compact.xyz, needle_box(1.0e3), compact.types);
  const auto q_ref = reference.q4q6(cutoff, 0);
  REQUIRE(q_ref.q6.size() == static_cast<std::size_t>(n));
  for (const double edge : {1.0e6, 1.0e12, 5.0e15}) {
    const Needle cloud = make_needle(edge, n, false);
    seams::System sys;
    sys.set_frame(cloud.xyz, needle_box(edge), cloud.types);
    seams::Q4Q6 got;
    try {
      got = sys.q4q6(cutoff, 0);
    } catch (const seams::Error &error) {
      REQUIRE(error.code() == SEAMS_ERR_RANGE);
      REQUIRE(edge >= 5.0e15);
      continue;
    }
    REQUIRE(got.q6.size() == q_ref.q6.size());
    double worst = 0.0;
    for (std::size_t i = 0; i < got.q6.size(); ++i) {
      worst = std::max(worst, std::fabs(got.q6[i] - q_ref.q6[i]));
      worst = std::max(worst, std::fabs(got.q4[i] - q_ref.q4[i]));
    }
    REQUIRE(worst < 1.0e-8);
  }
}

TEST_CASE("host-guest distances on a long needle match the compact cell", "[pore][system]") {
  constexpr int n = 260;
  constexpr double cutoff = 3.5;
  const Needle compact = make_needle(1.0e3, n, false);
  REQUIRE(seams_set_num_threads(1) == SEAMS_OK);
  seams::HostGuestOptions options;
  options.guest = 1;
  options.framework = 2;
  options.cutoff = cutoff;
  seams::System reference;
  reference.set_frame(compact.xyz, needle_box(1.0e3), compact.types);
  const auto base = reference.host_guest(options);
  const Needle cloud = make_needle(1.0e12, n, false);
  seams::System sys;
  sys.set_frame(cloud.xyz, needle_box(1.0e12), cloud.types);
  const auto got = sys.host_guest(options);
  REQUIRE(got.nearest_index == base.nearest_index);
  REQUIRE(got.nearest_distance.size() == base.nearest_distance.size());
  for (std::size_t i = 0; i < got.nearest_distance.size(); ++i) {
    REQUIRE(std::fabs(got.nearest_distance[i] - base.nearest_distance[i]) < 1.0e-6);
  }
}
