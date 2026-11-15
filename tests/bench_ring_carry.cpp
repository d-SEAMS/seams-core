/*
** Carried rings against a full search when the old full-pass rules would
** have fired. A tetrahedral cell is searched once, then a fixed share of
** the molecules is moved. The updater is given the atom fraction 0.2 and
** the bond crossover 0.15. Those rules used to replace the frame with a
** search from every vertex. The timed update searches only the vertices
** in T and must match ringNetwork.
**
**   bench_ring_carry [nAtoms] [samples]
*/

#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
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
    return (unit - 0.5) * 0.15 * spacing;
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

std::vector<std::vector<int>>
indexList(molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  auto nList = nneigh::neighListO(3.5, cloud, 1);
  return nneigh::neighbourListByIndex(cloud, nList);
}

double msSince(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - t0)
      .count();
}

} // namespace

int main(int argc, char **argv) {
  const int nAtoms = argc > 1 ? std::atoi(argv[1]) : 4096;
  const int samples = argc > 2 ? std::atoi(argv[2]) : 8;
  const bool local = argc > 3 && std::string(argv[3]) == "local";
  auto cloud = makeCloud(nAtoms, 88172645463325252ULL);
  auto first = indexList(cloud);
  primitive::RingUpdater warm(6);
  warm.setFullPassFraction(0.2);
  warm.setBondCrossover(0.15);
  warm.update(first);

  // local: a few molecules move, which is the carried-ring case.
  // The default moves four of every five, the rate at which the old
  // full-pass rule fired on the nucleation trajectory.
  const int moved = local ? std::min(nAtoms, 1200) : nAtoms;
  for (int i = 0; i < moved; i++) {
    if (!local && i % 5 == 0) {
      continue;
    }
    cloud.pts[static_cast<std::size_t>(i)].x += 1.6;
    cloud.pts[static_cast<std::size_t>(i)].y += 0.4;
  }
  auto idx = indexList(cloud);
  primitive::RingUpdater updater(6);
  updater.setFullPassFraction(0.2);
  updater.setBondCrossover(0.15);
  updater.update(first);
  const auto &once = updater.update(idx);
  const auto fullOnce = primitive::ringNetwork(idx, 6);
  const bool same = once == fullOnce;
  const int recomputed = updater.lastRecomputedSources();
  const int changedBonds = updater.lastChangedBonds();
  const double perMol = updater.lastChangedBondFraction();
  const bool oldRule = recomputed > static_cast<int>(0.2 * nAtoms) || perMol > 0.15;

  std::vector<double> fullMs;
  std::vector<double> carryMs;
  fullMs.reserve(static_cast<std::size_t>(samples));
  carryMs.reserve(static_cast<std::size_t>(samples));
  // One updater, already holding the moved frame. Each sample walks back
  // to the first frame and forward again, and times the forward step
  // against a full search of that same frame.
  updater.update(first);
  for (int i = 0; i < samples; i++) {
    const auto tFull = std::chrono::steady_clock::now();
    auto full = primitive::ringNetwork(idx, 6);
    fullMs.push_back(msSince(tFull));
    const auto tCarry = std::chrono::steady_clock::now();
    auto carried = updater.update(idx);
    carryMs.push_back(msSince(tCarry));
    if (full != carried) {
      std::cout << "identical no\n";
      return 1;
    }
    updater.update(first);
  }
  const auto fullMed = [&]() {
    auto v = fullMs;
    std::sort(v.begin(), v.end());
    return 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
  }();
  const auto carryMed = [&]() {
    auto v = carryMs;
    std::sort(v.begin(), v.end());
    return 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
  }();
  const double fullMin = *std::min_element(fullMs.begin(), fullMs.end());
  const double carryMin = *std::min_element(carryMs.begin(), carryMs.end());
  const char *verdict =
      (carryMed < fullMed && carryMin < fullMin) ? "WIN" : "LOSE";
  std::cout.setf(std::ios::fixed);
  std::cout.precision(4);
  std::cout << "molecules " << nAtoms << " samples " << samples << "\n"
            << "changed_bonds " << changedBonds << " per_molecule " << perMol
            << "\n"
            << "recomputed " << recomputed << " of " << nAtoms << "\n"
            << "old_full_pass_rule " << (oldRule ? "yes" : "no") << "\n"
            << "identical " << (same ? "yes" : "no") << "\n"
            << "full_median " << fullMed << " full_min " << fullMin << "\n"
            << "carry_median " << carryMed << " carry_min " << carryMin << " "
            << verdict << "\n";
  for (double v : carryMs) {
    std::cout << "carry_sample " << v << "\n";
  }
  for (double v : fullMs) {
    std::cout << "full_sample " << v << "\n";
  }
  return same ? 0 : 1;
}
