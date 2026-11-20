// Copyright (c) d-SEAMS developers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <seams.hpp>

#include <array>
#include <cmath>
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
