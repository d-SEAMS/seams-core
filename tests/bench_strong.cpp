/*
** This file is part of d-SEAMS.
**
** SPDX-License-Identifier: MIT
**
** Strong scaling of the Steinhardt kernel and of primitive rings at
** fixed N. Neighbour-list construction is timed separately and is
** still replicated on every rank; do not read a drop in t_neigh as an
** MPI or GPU win.
**
**   bench_strong [nAtoms] [reps]
**
** MPI: launch with mpirun -n R. OpenMP threads come from OMP_NUM_THREADS.
** Offload: build with -Dwith_openmp_offload=true and SEAMS_OFFLOAD=1.
*/

#include <bop.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef SEAMS_HAS_MPI
#include <mpi.h>
#endif

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
constexpr double kNumberDensity = 0.0332;

Cloud makeCloud(int nAtoms) {
  Cloud cloud;
  const double boxLength = std::cbrt(nAtoms / kNumberDensity);
  const int perSide = static_cast<int>(std::ceil(std::cbrt(nAtoms)));
  const double spacing = boxLength / perSide;
  cloud.nop = nAtoms;
  cloud.currentFrame = 1;
  cloud.box = {boxLength, boxLength, boxLength};
  cloud.boxLow = {0.0, 0.0, 0.0};
  unsigned long long state = 88172645463325252ULL;
  for (int i = 0; i < nAtoms; i++) {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    const double unit = static_cast<double>(state >> 11) / 9007199254740992.0;
    const double jit = (unit - 0.5) * 0.3 * spacing;
    molSys::Point<double> p;
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    p.x = ((i % perSide) + 0.5) * spacing + jit;
    p.y = (((i / perSide) % perSide) + 0.5) * spacing + jit;
    p.z = (((i / (perSide * perSide)) % perSide) + 0.5) * spacing + jit;
    cloud.pts.push_back(p);
    cloud.idIndexMap[p.atomID] = i;
  }
  return cloud;
}

template <typename Fn> double bestMillis(Fn &&fn, int reps) {
  double best = 1.0e300;
  for (int r = 0; r < reps; r++) {
    const auto t0 = std::chrono::steady_clock::now();
    fn();
    const auto t1 = std::chrono::steady_clock::now();
    best = std::min(
        best, std::chrono::duration<double, std::milli>(t1 - t0).count());
  }
  return best;
}

} // namespace

int main(int argc, char **argv) {
#ifdef SEAMS_HAS_MPI
  MPI_Init(&argc, &argv);
  int rank = 0;
  int nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
#else
  const int rank = 0;
  const int nranks = 1;
#endif

  const int nAtoms = argc > 1 ? std::atoi(argv[1]) : 65536;
  const int reps = argc > 2 ? std::atoi(argv[2]) : 3;
  int nThreads = 1;
#ifdef SEAMS_HAS_OPENMP
  nThreads = omp_get_max_threads();
#endif
  int nDevices = 0;
#ifdef SEAMS_HAS_OFFLOAD
  nDevices = omp_get_num_devices();
#endif

  Cloud cloud = makeCloud(nAtoms);
  std::vector<std::vector<int>> nList;
  const double tNeigh = bestMillis(
      [&]() { nList = nneigh::neighListO(3.5, cloud, 1); }, reps);

#ifdef SEAMS_HAS_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  const double tQl = bestMillis(
      [&]() {
        volatile auto q = chill::steinhardtQl(cloud, nList, 6);
        (void)q;
      },
      reps);

  std::vector<std::vector<int>> byIndex;
  const double tIndex = bestMillis(
      [&]() { byIndex = nneigh::getNewNeighbourListByIndex(cloud, 3.5); },
      reps);
  const double tRings = bestMillis(
      [&]() {
        volatile auto rings = primitive::ringNetwork(byIndex, 6);
        (void)rings;
      },
      reps);

  if (rank == 0) {
    std::cout << std::left << std::setw(10) << "nAtoms" << std::setw(8)
              << "ranks" << std::setw(8) << "thr" << std::setw(8) << "devs"
              << std::setw(16) << "neigh/ms" << std::setw(16) << "ql/ms"
              << std::setw(16) << "index/ms" << std::setw(16) << "rings/ms"
              << "\n";
    std::cout << std::left << std::setw(10) << nAtoms << std::setw(8) << nranks
              << std::setw(8) << nThreads << std::setw(8) << nDevices
              << std::setw(16) << std::fixed << std::setprecision(3) << tNeigh
              << std::setw(16) << tQl << std::setw(16) << tIndex
              << std::setw(16) << tRings << "\n";
    std::cout << "# neigh and index are replicated; ql and rings are threaded\n";
  }

#ifdef SEAMS_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}
