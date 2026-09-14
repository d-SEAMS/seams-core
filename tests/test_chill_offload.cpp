#include <bop.hpp>
#include <chill_offload.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

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
