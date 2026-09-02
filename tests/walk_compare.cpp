/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** CHILL+ labels and cage membership, with their largest ice clusters, on
** every frame of a trajectory in one pass.
**   walk_compare TRAJ [lastFrame] [atomType] [stride]
**
** Columns: frame nop
**   chill_cubic chill_hex chill_interfacial chill_clathrate
**   chill_interclathrate chill_water chill_ice chill_max chill_clus
**   cut_ice cut_max cut_clus
**   seed_ih seed_ic seed_both seed_ice seed_max seed_clus
** CHILL+ reads the cutoff graph. chill_ice is cubic plus hexagonal, the
** bulk-ice count of Nguyen and Molinero; interfacial molecules are not
** ice. cut_* is cage membership on the same cutoff graph; seed_* is the
** seeded assignment on the mutual and union four-nearest graphs.
*/

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace {

std::vector<std::vector<int>>
sixOf(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &r : rings) {
    if (r.size() == 6) {
      six.push_back(r);
    }
  }
  return six;
}

int findRoot(std::vector<int> &parent, int a) {
  while (parent[static_cast<std::size_t>(a)] != a) {
    parent[static_cast<std::size_t>(a)] =
        parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
    a = parent[static_cast<std::size_t>(a)];
  }
  return a;
}

void unite(std::vector<int> &parent, std::vector<int> &sz, int a, int b) {
  a = findRoot(parent, a);
  b = findRoot(parent, b);
  if (a == b) {
    return;
  }
  if (sz[static_cast<std::size_t>(a)] < sz[static_cast<std::size_t>(b)]) {
    std::swap(a, b);
  }
  parent[static_cast<std::size_t>(b)] = a;
  sz[static_cast<std::size_t>(a)] += sz[static_cast<std::size_t>(b)];
}

// Connected components of the flagged atoms over the index graph idx
void clusterFlags(const std::vector<char> &ice,
                  const std::vector<std::vector<int>> &idx, int &nIce,
                  int &nMax, int &nClus) {
  const int n = static_cast<int>(ice.size());
  nIce = 0;
  for (int i = 0; i < n; i++) {
    nIce += ice[static_cast<std::size_t>(i)] ? 1 : 0;
  }
  nMax = 0;
  nClus = 0;
  if (nIce == 0) {
    return;
  }
  std::vector<int> parent(static_cast<std::size_t>(n));
  std::vector<int> sz(static_cast<std::size_t>(n), 1);
  for (int i = 0; i < n; i++) {
    parent[static_cast<std::size_t>(i)] = i;
  }
  for (int i = 0; i < n; i++) {
    if (!ice[static_cast<std::size_t>(i)] ||
        static_cast<int>(idx.size()) <= i) {
      continue;
    }
    for (std::size_t k = 1; k < idx[static_cast<std::size_t>(i)].size();
         k++) {
      const int j = idx[static_cast<std::size_t>(i)][k];
      if (j > i && j < n && ice[static_cast<std::size_t>(j)]) {
        unite(parent, sz, i, j);
      }
    }
  }
  for (int i = 0; i < n; i++) {
    if (!ice[static_cast<std::size_t>(i)]) {
      continue;
    }
    if (findRoot(parent, i) == i) {
      ++nClus;
      nMax = std::max(nMax, sz[static_cast<std::size_t>(i)]);
    }
  }
}

std::vector<char> atomsFromRings(const std::vector<std::vector<int>> &six,
                                 const ring::CageAffiliation &aff, int nop) {
  std::vector<char> ice(static_cast<std::size_t>(nop), 0);
  for (std::size_t r = 0; r < six.size(); r++) {
    if (!aff.hc[r] && !aff.ddc[r]) {
      continue;
    }
    for (const int a : six[r]) {
      if (a >= 0 && a < nop) {
        ice[static_cast<std::size_t>(a)] = 1;
      }
    }
  }
  return ice;
}

struct ChillCounts {
  int cubic = 0, hex = 0, interfacial = 0, clathrate = 0, interClathrate = 0,
      water = 0;
};

ChillCounts
tallyChill(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
           std::vector<char> &bulkIce) {
  ChillCounts c;
  bulkIce.assign(static_cast<std::size_t>(cloud.nop), 0);
  for (int i = 0; i < cloud.nop; i++) {
    switch (cloud.pts[static_cast<std::size_t>(i)].iceType) {
    case molSys::atom_state_type::cubic:
    case molSys::atom_state_type::reCubic:
      ++c.cubic;
      bulkIce[static_cast<std::size_t>(i)] = 1;
      break;
    case molSys::atom_state_type::hexagonal:
    case molSys::atom_state_type::reHex:
      ++c.hex;
      bulkIce[static_cast<std::size_t>(i)] = 1;
      break;
    case molSys::atom_state_type::interfacial:
      ++c.interfacial;
      break;
    case molSys::atom_state_type::clathrate:
      ++c.clathrate;
      break;
    case molSys::atom_state_type::interClathrate:
      ++c.interClathrate;
      break;
    default:
      ++c.water;
      break;
    }
  }
  return c;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: walk_compare TRAJ [lastFrame] [atomType] [stride]\n");
    return 2;
  }
  const std::string traj = argv[1];
  const int want = argc > 2 ? std::atoi(argv[2]) : 0;
  const int typeI = argc > 3 ? std::atoi(argv[3]) : 1;
  const int stride = argc > 4 ? std::atoi(argv[4]) : 1;
  const int nframes = sinp::nLammpsFrames(traj);
  if (nframes <= 0) {
    std::fprintf(stderr, "no frames in %s\n", traj.c_str());
    return 1;
  }
  const int last = (want > 0 && want < nframes) ? want : nframes;
  const int step = stride > 0 ? stride : 1;

  std::printf("# %s nframes %d last %d stride %d\n", traj.c_str(), nframes,
              last, step);
  std::printf("# frame nop "
              "chill_cubic chill_hex chill_interfacial chill_clathrate "
              "chill_interclathrate chill_water chill_ice chill_max chill_clus "
              "cut_ice cut_max cut_clus "
              "seed_ih seed_ic seed_both seed_ice seed_max seed_clus\n");

  const double cutoff = 3.5;
  const double cand = 5.5;
  const int k = 4;
  for (int frame = 1; frame <= last; frame += step) {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI);
    if (cloud.nop == 0) {
      std::printf("%d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n", frame);
      continue;
    }
    const int nop = cloud.nop;

    // CHILL+ on the cutoff graph
    auto cutRows = nneigh::neighListO(cutoff, cloud, typeI);
    auto idxC = nneigh::neighbourListByIndex(cloud, cutRows);
    chill::getCorrelPlus(cloud, cutRows, false);
    chill::getIceTypePlusNoPrint(cloud, cutRows, false);
    std::vector<char> chillIce;
    const ChillCounts cc = tallyChill(cloud, chillIce);
    int chillN = 0, chillMax = 0, chillClus = 0;
    clusterFlags(chillIce, idxC, chillN, chillMax, chillClus);

    // Cage membership on the same cutoff graph
    auto sixC = sixOf(primitive::ringNetwork(idxC, 6));
    const auto affC = ring::cageAffiliation(sixC, idxC);
    auto cutIce = atomsFromRings(sixC, affC, nop);
    int cutN = 0, cutMax = 0, cutClus = 0;
    clusterFlags(cutIce, idxC, cutN, cutMax, cutClus);

    // Seeded assignment on the four-nearest graphs
    auto graphs = nneigh::kNearestNeighbourPair(cloud, k, cand, typeI);
    auto idxS = nneigh::neighbourListByIndex(cloud, graphs.first);
    auto idxU = nneigh::neighbourListByIndex(cloud, graphs.second);
    auto sixS = sixOf(primitive::ringNetwork(idxS, 6));
    auto sixU = sixOf(primitive::ringNetwork(idxU, 6));
    const auto seeded = ring::seededCageAffiliation(sixS, idxS, sixU, idxU);
    int ih = 0, ic = 0, both = 0;
    std::vector<char> seedIce(static_cast<std::size_t>(nop), 0);
    for (int i = 0; i < nop; i++) {
      const bool h = seeded.hc[static_cast<std::size_t>(i)];
      const bool d = seeded.ddc[static_cast<std::size_t>(i)];
      if (h && d) {
        ++both;
      } else if (h) {
        ++ih;
      } else if (d) {
        ++ic;
      }
      seedIce[static_cast<std::size_t>(i)] = (h || d) ? 1 : 0;
    }
    int seedN = 0, seedMax = 0, seedClus = 0;
    clusterFlags(seedIce, idxU, seedN, seedMax, seedClus);

    std::printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
                frame, nop, cc.cubic, cc.hex, cc.interfacial, cc.clathrate,
                cc.interClathrate, cc.water, chillN, chillMax, chillClus, cutN,
                cutMax, cutClus, ih, ic, both, seedN, seedMax, seedClus);
    std::fflush(stdout);
  }
  return 0;
}
