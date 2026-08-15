/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Seeded TUM ice score and largest ice cluster on every frame.
**   walk_seeded TRAJ [lastFrame] [atomType] [stride]
*/

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

void clusterIce(const ring::SeededAtomLabels &lab,
                const std::vector<std::vector<int>> &idxU, int &nMax,
                int &nClus) {
  const int n = static_cast<int>(lab.hc.size());
  std::vector<char> ice(static_cast<std::size_t>(n), 0);
  int nIce = 0;
  for (int i = 0; i < n; i++) {
    if (lab.hc[static_cast<std::size_t>(i)] ||
        lab.ddc[static_cast<std::size_t>(i)]) {
      ice[static_cast<std::size_t>(i)] = 1;
      ++nIce;
    }
  }
  if (nIce == 0) {
    nMax = 0;
    nClus = 0;
    return;
  }
  std::vector<int> parent(static_cast<std::size_t>(n));
  std::vector<int> sz(static_cast<std::size_t>(n), 1);
  for (int i = 0; i < n; i++) {
    parent[static_cast<std::size_t>(i)] = i;
  }
  for (int i = 0; i < n; i++) {
    if (!ice[static_cast<std::size_t>(i)] ||
        static_cast<int>(idxU.size()) <= i) {
      continue;
    }
    for (std::size_t k = 1; k < idxU[static_cast<std::size_t>(i)].size();
         k++) {
      const int j = idxU[static_cast<std::size_t>(i)][k];
      if (j > i && j < n && ice[static_cast<std::size_t>(j)]) {
        unite(parent, sz, i, j);
      }
    }
  }
  nMax = 0;
  nClus = 0;
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

void countIce(const ring::SeededAtomLabels &lab, int &ih, int &ic, int &both,
              int &water) {
  ih = ic = both = water = 0;
  const int n = static_cast<int>(lab.hc.size());
  for (int i = 0; i < n; i++) {
    const bool h = lab.hc[static_cast<std::size_t>(i)];
    const bool d = lab.ddc[static_cast<std::size_t>(i)];
    if (h && d) {
      ++both;
    } else if (h) {
      ++ih;
    } else if (d) {
      ++ic;
    } else {
      ++water;
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: walk_seeded TRAJ [lastFrame] [atomType] [stride]\n");
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
  std::printf("# frame nop ih ic both water n_ice n_max n_clus cubicity\n");

  const double cand = 5.5;
  const int k = 4;
  for (int frame = 1; frame <= last; frame += step) {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI);
    if (cloud.nop == 0) {
      std::fprintf(stderr, "empty frame %d\n", frame);
      return 1;
    }
    auto mutual = nneigh::kNearestNeighbourList(cloud, k, cand, typeI, true);
    auto uni = nneigh::kNearestNeighbourList(cloud, k, cand, typeI, false);
    auto idxS = nneigh::neighbourListByIndex(cloud, mutual);
    auto idxU = nneigh::neighbourListByIndex(cloud, uni);
    auto sixS = sixOf(primitive::ringNetwork(idxS, 6));
    auto sixU = sixOf(primitive::ringNetwork(idxU, 6));
    const auto seeded = ring::seededCageAffiliation(sixS, idxS, sixU, idxU);
    int ih = 0;
    int ic = 0;
    int both = 0;
    int water = 0;
    countIce(seeded, ih, ic, both, water);
    int nMax = 0;
    int nClus = 0;
    clusterIce(seeded, idxU, nMax, nClus);
    const int nIce = ih + ic + both;
    const double cub =
        nIce > 0 ? static_cast<double>(ic + both) / nIce : 0.0;
    std::printf("%d %d %d %d %d %d %d %d %d %.4f\n", frame, cloud.nop, ih, ic,
                both, water, nIce, nMax, nClus, cub);
  }
  return 0;
}
