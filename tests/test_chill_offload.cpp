#include <bop.hpp>
#include <chill_offload.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <cmath>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("specialized CHILL+ matches host getCorrelPlus on mixed TIP4P",
          "[chill][offload]") {
  molSys::PointCloud<molSys::Point<double>, double> host;
  host = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, host, 2);
  REQUIRE(host.nop > 0);
  auto nList = nneigh::kNearestNeighbourList(host, 4, 5.5, 2, true);
  auto specCloud = host;
  chill::getCorrelPlus(host, nList, false);
  chill::getIceTypePlusNoPrint(host, nList, false);
  const auto spec = chill::specializedChillPlus(specCloud, nList);
  UNSCOPED_INFO("usedDevice=" << spec.usedDevice);
  REQUIRE(spec.iceType.size() == static_cast<std::size_t>(host.nop));
  int mismatch = 0;
  for (int i = 0; i < host.nop; i++) {
    if (spec.iceType[static_cast<std::size_t>(i)] !=
        static_cast<int>(host.pts[static_cast<std::size_t>(i)].iceType)) {
      ++mismatch;
    }
  }
  REQUIRE(mismatch == 0);
}

TEST_CASE("specialized CHILL+ matches host on primitive hexagonal triclinic Ice Ih",
          "[chill][offload][triclinic]") {
  // Primitive hexagonal Ice Ih supercell (3x3x2, 72 oxygen sites) with tilt xy = -0.5 * lx.
  const double bond = 2.75;
  const double a = bond * std::sqrt(8.0 / 3.0);
  const double cLat = bond * (8.0 / 3.0);
  const int nx = 3;
  const int ny = 3;
  const int nz = 2;
  const double lx = nx * a;
  const double ly = ny * a * std::sqrt(3.0) * 0.5;
  const double lz = nz * cLat;
  const double xy = -0.5 * ny * a;
  const double xloBound = std::min(0.0, xy);
  const double xhiBound = lx + std::max(0.0, xy);

  const double basis[4][3] = {
      {1.0 / 3.0, 2.0 / 3.0, 1.0 / 16.0},
      {2.0 / 3.0, 1.0 / 3.0, 9.0 / 16.0},
      {1.0 / 3.0, 2.0 / 3.0, 7.0 / 16.0},
      {2.0 / 3.0, 1.0 / 3.0, 15.0 / 16.0},
  };

  molSys::PointCloud<molSys::Point<double>, double> host;
  host.box = {xhiBound - xloBound, ly, lz, xy, 0.0, 0.0};
  host.boxLow = {xloBound, 0.0, 0.0};
  host.currentFrame = 1;
  int id = 1;
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        for (int b = 0; b < 4; ++b) {
          const double fa = (ix + basis[b][0]) / nx;
          const double fb = (iy + basis[b][1]) / ny;
          const double fc = (iz + basis[b][2]) / nz;
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = lx * fa + xy * fb;
          p.y = ly * fb;
          p.z = lz * fc;
          host.idIndexMap[id] = static_cast<int>(host.pts.size());
          host.pts.push_back(p);
          ++id;
        }
      }
    }
  }
  host.nop = static_cast<int>(host.pts.size());
  REQUIRE(host.nop == 72);

  const auto nList = nneigh::neighListO(3.5, host, 1);
  auto specCloud = host;
  const auto hRes = chill::hostChillPlus(host, nList);
  const auto sRes = chill::specializedChillPlus(specCloud, nList);
  REQUIRE(hRes.iceType.size() == 72);
  REQUIRE(sRes.iceType.size() == 72);

  int mismatch = 0;
  for (int i = 0; i < host.nop; ++i) {
    const auto u = static_cast<std::size_t>(i);
    REQUIRE(hRes.iceType[u] == static_cast<int>(molSys::atom_state_type::hexagonal));
    if (sRes.iceType[u] != hRes.iceType[u]) {
      ++mismatch;
    }
  }
  REQUIRE(mismatch == 0);
}
