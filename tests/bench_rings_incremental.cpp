/*
** Incremental RingUpdater versus a full ringNetwork, on a jittered
** tetrahedral cell at the mW density.  Three atoms are displaced each
** synthetic frame.  Run: bcpp/tests/bench_rings_incremental [nAtoms]
** [frames] [sites].
*/

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

namespace {

constexpr double kDensity = 0.0332;

molSys::PointCloud<molSys::Point<double>, double>
makeCloud(int nAtoms, unsigned long long seed) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double boxLength = std::cbrt(nAtoms / kDensity);
  const int perSide = static_cast<int>(std::ceil(std::cbrt(nAtoms)));
  const double spacing = boxLength / perSide;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  unsigned long long state = seed;
  auto jitter = [&state, spacing]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    const double unit = static_cast<double>(state >> 11) / 9007199254740992.0;
    return (unit - 0.5) * 0.30 * spacing;
  };
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    p.x = ((i % perSide) + 0.5) * spacing + jitter();
    p.y = (((i / perSide) % perSide) + 0.5) * spacing + jitter();
    p.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing + jitter();
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  return cloud;
}

std::vector<std::vector<int>> indexList(
    molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  return nneigh::neighbourListByIndex(cloud, nList);
}

std::set<std::vector<int>> canon(const std::vector<std::vector<int>> &rings) {
  std::set<std::vector<int>> out;
  for (auto r : rings) {
    auto it = std::min_element(r.begin(), r.end());
    std::rotate(r.begin(), it, r.end());
    out.insert(r);
  }
  return out;
}

} // namespace

int main(int argc, char **argv) {
  const int nAtoms = argc > 1 ? std::atoi(argv[1]) : 8000;
  const int frames = argc > 2 ? std::atoi(argv[2]) : 5;
  const int sites = argc > 3 ? std::atoi(argv[3]) : 3;
  const int maxDepth = 6;

  auto cloud = makeCloud(nAtoms, 88172645463325252ULL);
  primitive::RingUpdater updater(maxDepth);

  auto idx = indexList(cloud);
  updater.update(idx);
  (void)primitive::ringNetwork(idx, maxDepth);

  unsigned long long state = 424242ULL;
  auto pick = [&state, nAtoms]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return static_cast<int>((state >> 11) %
                            static_cast<unsigned long long>(nAtoms));
  };

  double bestFull = 1e300;
  double bestIncr = 1e300;
  int lastRecomp = 0;
  int lastBalls = 0;
  bool match = true;

  for (int frame = 0; frame < frames; frame++) {
    for (int k = 0; k < sites; k++) {
      auto &p = cloud.pts[pick()];
      p.x += 0.8 * ((frame + k) % 2 ? 1.0 : -1.0);
      p.y += 0.4;
    }
    idx = indexList(cloud);

    const auto t0 = std::chrono::steady_clock::now();
    auto full = primitive::ringNetwork(idx, maxDepth);
    const auto t1 = std::chrono::steady_clock::now();
    auto incr = updater.update(idx);
    const auto t2 = std::chrono::steady_clock::now();

    const double fullMs =
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double incrMs =
        std::chrono::duration<double, std::milli>(t2 - t1).count();
    bestFull = std::min(bestFull, fullMs);
    bestIncr = std::min(bestIncr, incrMs);
    lastRecomp = updater.lastRecomputedSources();
    lastBalls = updater.lastBallsRefreshed();
    if (canon(full) != canon(incr)) {
      match = false;
    }
  }

  std::cout << "nAtoms     " << nAtoms << "\n"
            << "frames     " << frames << "\n"
            << "sites      " << sites << "\n"
            << "full/ms    " << std::fixed << std::setprecision(3) << bestFull
            << "\n"
            << "incr/ms    " << bestIncr << "\n"
            << "ratio      " << std::setprecision(2)
            << (bestIncr > 0.0 ? bestFull / bestIncr : 0.0) << "\n"
            << "recomp     " << lastRecomp << " / " << nAtoms << "\n"
            << "balls      " << lastBalls << " / " << nAtoms << "\n"
            << "identical  " << (match ? "yes" : "NO") << "\n";
  return match ? 0 : 1;
}
