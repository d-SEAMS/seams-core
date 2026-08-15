/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Ice cluster size for every published bond graph on each frame.
**   walk_graphs TRAJ [lastFrame] [atomType] [stride]
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

void clusterFlags(const std::vector<bool> &hc, const std::vector<bool> &ddc,
                  const std::vector<std::vector<int>> &idx, int &nIce,
                  int &nMax, int &nClus) {
  const int n = static_cast<int>(hc.size());
  std::vector<char> ice(static_cast<std::size_t>(n), 0);
  nIce = 0;
  for (int i = 0; i < n; i++) {
    if (hc[static_cast<std::size_t>(i)] || ddc[static_cast<std::size_t>(i)]) {
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
        static_cast<int>(idx.size()) <= i) {
      continue;
    }
    for (std::size_t k = 1; k < idx[static_cast<std::size_t>(i)].size(); k++) {
      const int j = idx[static_cast<std::size_t>(i)][k];
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

void atomsFromRings(const std::vector<std::vector<int>> &six,
                    const ring::CageAffiliation &aff, int nop,
                    std::vector<bool> &hc, std::vector<bool> &ddc) {
  hc.assign(static_cast<std::size_t>(nop), false);
  ddc.assign(static_cast<std::size_t>(nop), false);
  for (std::size_t r = 0; r < six.size(); r++) {
    for (const int a : six[r]) {
      if (a >= 0 && a < nop) {
        hc[static_cast<std::size_t>(a)] =
            hc[static_cast<std::size_t>(a)] || aff.hc[r];
        ddc[static_cast<std::size_t>(a)] =
            ddc[static_cast<std::size_t>(a)] || aff.ddc[r];
      }
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: walk_graphs TRAJ [lastFrame] [atomType] [stride]\n");
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
  std::printf("# frame "
              "cut_ice cut_max cut_clus "
              "knn_ice knn_max knn_clus "
              "uni_ice uni_max uni_clus "
              "seed_ice seed_max seed_clus\n");

  const double cutoff = 3.5;
  const double cand = 5.5;
  const int k = 4;
  for (int frame = 1; frame <= last; frame += step) {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI);
    if (cloud.nop == 0) {
      std::printf("%d 0 0 0 0 0 0 0 0 0 0 0 0\n", frame);
      continue;
    }
    const int nop = cloud.nop;
    auto cutRows = nneigh::neighListO(cutoff, cloud, typeI);
    auto knnRows = nneigh::kNearestNeighbourList(cloud, k, cand, typeI, true);
    auto uniRows = nneigh::kNearestNeighbourList(cloud, k, cand, typeI, false);
    auto idxC = nneigh::neighbourListByIndex(cloud, cutRows);
    auto idxK = nneigh::neighbourListByIndex(cloud, knnRows);
    auto idxU = nneigh::neighbourListByIndex(cloud, uniRows);
    auto sixC = sixOf(primitive::ringNetwork(idxC, 6));
    auto sixK = sixOf(primitive::ringNetwork(idxK, 6));
    auto sixU = sixOf(primitive::ringNetwork(idxU, 6));
    const auto affC = ring::cageAffiliation(sixC, idxC);
    const auto affK = ring::cageAffiliation(sixK, idxK);
    const auto affU = ring::cageAffiliation(sixU, idxU);
    const auto seeded = ring::seededCageAffiliation(sixK, idxK, sixU, idxU);

    std::vector<bool> hc, ddc;
    int ice = 0;
    int mx = 0;
    int cl = 0;
    std::printf("%d", frame);
    atomsFromRings(sixC, affC, nop, hc, ddc);
    clusterFlags(hc, ddc, idxC, ice, mx, cl);
    std::printf(" %d %d %d", ice, mx, cl);
    atomsFromRings(sixK, affK, nop, hc, ddc);
    clusterFlags(hc, ddc, idxK, ice, mx, cl);
    std::printf(" %d %d %d", ice, mx, cl);
    atomsFromRings(sixU, affU, nop, hc, ddc);
    clusterFlags(hc, ddc, idxU, ice, mx, cl);
    std::printf(" %d %d %d", ice, mx, cl);
    clusterFlags(seeded.hc, seeded.ddc, idxU, ice, mx, cl);
    std::printf(" %d %d %d\n", ice, mx, cl);
  }
  return 0;
}
