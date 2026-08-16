/*
** Pipeline wall times on the same synthetic mW-density cells as
** bench_scaling: cutoff list, k-nearest pair, primitive rings, cage
** affiliation. Every size in the sweep is timed. Run from input/.
*/

#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
constexpr double kNumberDensity = 0.0332;

Cloud makeCloud(int nAtoms, int typeI) {
  Cloud cloud;
  const double boxLength = std::cbrt(nAtoms / kNumberDensity);
  const int perSide = static_cast<int>(std::ceil(std::cbrt(nAtoms)));
  const double spacing = boxLength / perSide;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  unsigned long long state = 88172645463325252ULL;
  auto jitter = [&state, spacing]() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    const double unit = static_cast<double>(state >> 11) / 9007199254740992.0;
    return (unit - 0.5) * 0.3 * spacing;
  };
  for (int i = 0; i < nAtoms; i++) {
    molSys::Point<double> point;
    point.type = typeI;
    point.atomID = i + 1;
    point.molID = i + 1;
    point.x = ((i % perSide) + 0.5) * spacing + jitter();
    point.y = (((i / perSide) % perSide) + 0.5) * spacing + jitter();
    point.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing + jitter();
    cloud.pts.push_back(point);
    cloud.idIndexMap[point.atomID] = i;
  }
  return cloud;
}

template <typename Fn> double bestMillis(Fn &&fn, int reps) {
  double best = std::numeric_limits<double>::max();
  for (int rep = 0; rep < reps; rep++) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    best = std::min(
        best, std::chrono::duration<double, std::milli>(end - start).count());
  }
  return best;
}

} // namespace

int main(int argc, char **argv) {
  std::vector<int> sizes;
  for (int arg = 1; arg < argc; arg++) {
    sizes.push_back(std::atoi(argv[arg]));
  }
  if (sizes.empty()) {
    sizes = {1000, 2000, 4000, 8000, 16000, 32000};
  }

  constexpr int typeI = 1;
  constexpr double cutoff = 3.5;

  std::cout << std::left << std::setw(10) << "nAtoms" << std::setw(14)
            << "neigh/ms" << std::setw(14) << "knn/ms" << std::setw(14)
            << "knnPair/ms" << std::setw(14) << "rings/ms" << std::setw(14)
            << "affil/ms"
            << "\n";
  std::cout << std::string(80, '-') << "\n";

  for (int nAtoms : sizes) {
    Cloud cloud = makeCloud(nAtoms, typeI);
    const int reps = nAtoms <= 4000 ? 5 : 3;

    std::vector<std::vector<int>> nList;
    const double tNeigh = bestMillis(
        [&]() { nList = nneigh::neighListO(cutoff, cloud, typeI); }, reps);
    std::vector<std::vector<int>> idx =
        nneigh::neighbourListByIndex(cloud, nList);

    const double tKnn = bestMillis(
        [&]() {
          auto knn = nneigh::kNearestNeighbourList(cloud, 4, 5.5, typeI, true);
          (void)knn;
        },
        reps);
    const double tPair = bestMillis(
        [&]() {
          auto both = nneigh::kNearestNeighbourPair(cloud, 4, 5.5, typeI);
          (void)both;
        },
        reps);

    std::vector<std::vector<int>> rings;
    const double tRings =
        bestMillis([&]() { rings = primitive::ringNetwork(idx, 6); }, reps);
    std::vector<std::vector<int>> six;
    for (const auto &r : rings) {
      if (r.size() == 6) {
        six.push_back(r);
      }
    }
    ring::CageAffiliation aff;
    const double tAff =
        bestMillis([&]() { aff = ring::cageAffiliation(six, idx); }, reps);

    std::cout << std::left << std::setw(10) << nAtoms << std::fixed
              << std::setprecision(3) << std::setw(14) << tNeigh << std::setw(14)
              << tKnn << std::setw(14) << tPair << std::setw(14) << tRings
              << std::setw(14) << tAff << "\n";
  }
  return 0;
}
