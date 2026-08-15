/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Emits the per-atom cage classification for a trajectory frame as CSV, for
** figure generation.  Columns: x,y,z,icetype where icetype is 0 for an atom in
** no cage, 1 hexagonal cage, 2 double-diamond cage, 3 shared between the two.
**
** Build target: gen_figure_data.  Run from the input/ directory.
*/

#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  const std::string traj = argc > 1 ? argv[1] : "traj/mW_cubic.lammpstrj";
  const int frame = argc > 2 ? std::atoi(argv[2]) : 1;
  const std::string out = argc > 3 ? argv[3] : "cages.csv";

  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  yCloud = sinp::readLammpsTrjO(traj, frame, yCloud, 1);
  if (yCloud.nop == 0) {
    std::cerr << "could not read " << traj << "\n";
    return 1;
  }

  // mW is monatomic: the ring network comes from the oxygen neighbour list
  auto nList = nneigh::neighListO(3.5, yCloud, 1);
  auto nListIdx = nneigh::neighbourListByIndex(yCloud, nList);
  auto rings = primitive::ringNetwork(nListIdx, 7);

  std::vector<std::vector<int>> sixRings;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      sixRings.push_back(r);
    }
  }

  std::vector<ring::strucType> ringType(sixRings.size(),
                                        ring::strucType::unclassified);
  std::vector<cage::Cage> cageList;
  const auto index =
      ring::buildRingSearchIndex(sixRings, static_cast<int>(nListIdx.size()));

  auto listHC = ring::findHC(sixRings, ringType, nListIdx, cageList, index);
  auto listDDC = ring::findDDC(sixRings, ringType, listHC, cageList, index);
  auto listMixed = ring::findMixedRings(sixRings, ringType, listDDC, listHC);

  std::vector<cage::iceType> atomTypes(yCloud.nop, cage::iceType::dummy);
  ring::getAtomTypesTopoBulk(sixRings, ringType, atomTypes);

  int numHC = 0, numDDC = 0, mixedRings = 0, prismaticRings = 0, basalRings = 0;
  ring::getStrucNumbers(ringType, cageList, numHC, numDDC, mixedRings,
                        basalRings, prismaticRings);

  std::ofstream file(out);
  if (!file) {
    std::cerr << "could not open " << out << " for writing\n";
    return 1;
  }
  file << "x,y,z,icetype\n";
  int counts[4] = {0, 0, 0, 0};
  for (int i = 0; i < yCloud.nop; i++) {
    int code = 0;
    switch (atomTypes[i]) {
    case cage::iceType::hc:
      code = 1;
      break;
    case cage::iceType::ddc:
      code = 2;
      break;
    case cage::iceType::mixed:
    case cage::iceType::mixed2:
      code = 3;
      break;
    default:
      code = 0;
      break;
    }
    counts[code]++;
    file << yCloud.pts[i].x << "," << yCloud.pts[i].y << ","
         << yCloud.pts[i].z << "," << code << "\n";
  }

  // Also emit the ring membership of one double-diamond cage, so a figure can
  // show the unit the search actually identifies rather than only a bulk
  // colouring. A DDC is stored as its equatorial ring followed by six
  // peripheral rings.
  {
    std::ofstream cageFile(out + ".cage");
    cageFile << "cage_ring,vertex_order,atom_index,x,y,z\n";
    int emitted = 0;
    for (const auto &cage : cageList) {
      if (cage.type != cage::cageType::DoubleDiaC || cage.rings.size() != 7) {
        continue;
      }
      for (size_t r = 0; r < cage.rings.size(); r++) {
        const int ringIdx = cage.rings[r];
        if (ringIdx < 0 || ringIdx >= static_cast<int>(sixRings.size())) {
          continue;
        }
        for (size_t v = 0; v < sixRings[ringIdx].size(); v++) {
          const int a = sixRings[ringIdx][v];
          cageFile << r << "," << v << "," << a << "," << yCloud.pts[a].x << ","
                   << yCloud.pts[a].y << "," << yCloud.pts[a].z << "\n";
        }
      }
      emitted = 1;
      break; // one representative cage is enough for a figure
    }
    std::cout << "cage file  " << (emitted ? out + ".cage" : std::string("none"))
              << "\n";
  }

  std::cout << "atoms      " << yCloud.nop << "\n"
            << "box        " << yCloud.box[0] << " x " << yCloud.box[1] << " x "
            << yCloud.box[2] << "\n"
            << "six-rings  " << sixRings.size() << "\n"
            << "HC cages   " << numHC << "\n"
            << "DDC cages  " << numDDC << "\n"
            << "mixed ring " << mixedRings << "\n"
            << "-- atoms by classification --\n"
            << "  none     " << counts[0] << "\n"
            << "  HC       " << counts[1] << "\n"
            << "  DDC      " << counts[2] << "\n"
            << "  shared   " << counts[3] << "\n"
            << "wrote      " << out << "\n";

  return 0;
}
