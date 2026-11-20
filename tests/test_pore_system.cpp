// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <seams.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
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
