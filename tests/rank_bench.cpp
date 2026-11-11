// Time the Hilbert domain ring reduction. Rank 0 prints one millisecond
// per timed call. The first call is the warm-up and is not printed.
// mpiexec -n N sets the rank count. OMP_NUM_THREADS is the width inside
// each rank.

#include <domain.hpp>
#include <mol_sys.hpp>

#include <mpi.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace {

molSys::PointCloud<molSys::Point<double>, double> ice(int cells) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double side = cells * lattice;
  const double basis[8][3] = {{0, 0, 0},         {0.5, 0.5, 0},     {0.5, 0, 0.5},
                              {0, 0.5, 0.5},     {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  int id = 1;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (const auto &b : basis) {
          molSys::Point<double> p;
          p.type = 1;
          p.atomID = id;
          p.molID = id;
          p.x = (i + b[0]) * lattice;
          p.y = (j + b[1]) * lattice;
          p.z = (k + b[2]) * lattice;
          cloud.pts.push_back(p);
          cloud.idIndexMap[id] = id - 1;
          id++;
        }
      }
    }
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  cloud.currentFrame = 1;
  cloud.box = {side, side, side};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const int cells = argc > 1 ? std::atoi(argv[1]) : 24;
  const int reps = argc > 2 ? std::atoi(argv[2]) : 6;
  const auto cloud = ice(cells);
  seams::domain::reduceRings(cloud, 3.5, 6, MPI_COMM_WORLD);
  for (int n = 0; n < reps; n++) {
    const auto t0 = std::chrono::steady_clock::now();
    seams::domain::reduceRings(cloud, 3.5, 6, MPI_COMM_WORLD);
    const double ms = std::chrono::duration<double, std::milli>(
                          std::chrono::steady_clock::now() - t0)
                          .count();
    if (rank == 0) {
      std::printf("%.6f\n", ms);
    }
  }
  MPI_Finalize();
  return 0;
}
