//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

/*
** Wall time of classifyTemplates, soapSpectrumAll and voronoiFeatures on
** an FCC supercell.  Not in the regular test suite.
**
**   bench_structure_desc [n]
**
** n is the number of FCC conventional-cell repeats (default 4).  OpenMP
** is pinned to one thread so the printed times are a serial baseline.
*/

#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <structure_desc.hpp>

#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud fcc(int reps, double a) {
  Cloud cloud;
  const std::array<std::array<double, 3>, 4> basis = {
      {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}}};
  int id = 1;
  for (int i = 0; i < reps; i++) {
    for (int j = 0; j < reps; j++) {
      for (int k = 0; k < reps; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = (i + b[0]) * a;
          p.y = (j + b[1]) * a;
          p.z = (k + b[2]) * a;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  const double L = reps * a;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

template <typename Fn> double millis(Fn &&fn) {
  const auto t0 = std::chrono::steady_clock::now();
  fn();
  const auto t1 = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

} // namespace

int main(int argc, char **argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 4;
#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  if (std::getenv("OMP_NUM_THREADS") == nullptr) {
    omp_set_num_threads(1);
  }
#endif

  auto cloud = fcc(n, 4.0);
  auto nList = nneigh::neighListO(3.2, cloud, 1);

  const double tClassify =
      millis([&]() { (void)chill::classifyTemplates(cloud, nList, 12); });
  const double tSoap =
      millis([&]() { (void)chill::soapSpectrumAll(cloud, nList, 3, 6, 3.2); });
  const double tVoronoi =
      millis([&]() { (void)chill::voronoiFeatures(cloud, 4.8); });

  std::cout << "n                       " << n << "\n"
            << "nop                     " << cloud.nop << "\n"
            << "classifyTemplates/ms    " << std::fixed << std::setprecision(3)
            << tClassify << "\n"
            << "soapSpectrumAll/ms      " << tSoap << "\n"
            << "voronoiFeatures/ms      " << tVoronoi << "\n";
  return 0;
}
