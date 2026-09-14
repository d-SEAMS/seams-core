#include <bop.hpp>
#include <chill_offload.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <iostream>

int main() {
  molSys::PointCloud<molSys::Point<double>, double> host;
  host = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, host, 2);
  if (host.nop <= 0) {
    std::cerr << "empty cloud\n";
    return 2;
  }
  auto nList = nneigh::kNearestNeighbourList(host, 4, 5.5, 2, true);
  auto specCloud = host;
  chill::getCorrelPlus(host, nList, false);
  chill::getIceTypePlusNoPrint(host, nList, false);
  const auto spec = chill::specializedChillPlus(specCloud, nList);
  int mismatch = 0;
  for (int i = 0; i < host.nop; i++) {
    if (spec.iceType[static_cast<std::size_t>(i)] !=
        static_cast<int>(host.pts[static_cast<std::size_t>(i)].iceType)) {
      ++mismatch;
    }
  }
  std::cout << "nop=" << host.nop << " mismatch=" << mismatch
            << " usedDevice=" << spec.usedDevice << "\n";
  if (mismatch != 0) {
    return 1;
  }
  return spec.usedDevice ? 0 : 3;
}
