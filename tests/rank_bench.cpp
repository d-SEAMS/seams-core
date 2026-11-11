// Time one distributed operation. Rank 0 prints, per timed call,
// wall milliseconds, ghost-exchange milliseconds, compute milliseconds,
// and the largest rank's halo fraction. The first call is the warm-up.

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <domain.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <mpi.h>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace {

molSys::PointCloud<molSys::Point<double>, double> ice(int mx, int my, int mz) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const double bond = 2.75;
  const double lattice = 4.0 * bond / std::sqrt(3.0);
  const double sideX = mx * lattice;
  const double sideY = my * lattice;
  const double sideZ = mz * lattice;
  const double basis[8][3] = {{0, 0, 0},         {0.5, 0.5, 0},     {0.5, 0, 0.5},
                              {0, 0.5, 0.5},     {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  int id = 1;
  for (int i = 0; i < mx; i++) {
    for (int j = 0; j < my; j++) {
      for (int k = 0; k < mz; k++) {
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
  cloud.box = {sideX, sideY, sideZ};
  cloud.boxLow = {0.0, 0.0, 0.0};
  return cloud;
}

double haloOf(const std::string &op, double cutoff) {
  if (op == "rings") {
    return seams::domain::ringHops(6) * cutoff;
  }
  if (op == "chill" || op == "q6") {
    // The value at a neighbour needs that neighbour's own shell.
    return 2.0 * cutoff;
  }
  return cutoff;
}

void finish(const std::string &op, seams::domain::Work &work,
            const std::vector<std::vector<int>> &rows) {
  if (op == "chill") {
    chill::applyChillPlus(work.cloud, rows, chill::ChillPlusRule::Published);
  } else if (op == "q6") {
    auto q = chill::steinhardtQl(work.cloud, rows, 6);
    if (q.ql.empty() && work.cloud.nop > 0) {
      std::fprintf(stderr, "empty q6\n");
    }
  }
}

void compute(const std::string &op, seams::domain::Work &work, double cutoff) {
  if (op == "rings") {
    auto rings = seams::domain::ownedRings(work.cloud, cutoff, 6, work.share.owned);
    if (rings.empty() && work.share.owned.size() > 1000) {
      std::fprintf(stderr, "no rings\n");
    }
    return;
  }
  // Several ranks only fill rows for the owned shell. The cell grid is
  // the one reusedShellList keeps. One rank stays on the fused list.
  bool single = true;
  for (const char owned : work.share.owned) {
    if (!owned) {
      single = false;
      break;
    }
  }
  if (single || work.cloud.nop < 40000) {
    finish(op, work, nneigh::neighListO(cutoff, work.cloud, 1));
    return;
  }
  std::vector<char> source;
  seams::domain::markSources(work.cloud, work.share, op == "neigh" ? -1.0 : cutoff,
                             source);
  finish(op, work, nneigh::reusedShellList(cutoff, work.cloud, &source));
}

} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#ifdef SEAMS_HAS_OPENMP
  if (const char *threads = std::getenv("OMP_NUM_THREADS")) {
    const int n = std::atoi(threads);
    if (n > 0) {
      omp_set_num_threads(n);
    }
  }
#endif
  int mx = 24, my = 24, mz = 24, reps = 6;
  std::string op = "rings";
  if (argc >= 6) {
    mx = std::atoi(argv[1]);
    my = std::atoi(argv[2]);
    mz = std::atoi(argv[3]);
    reps = std::atoi(argv[4]);
    op = argv[5];
  } else {
    const int cells = argc > 1 ? std::atoi(argv[1]) : 24;
    mx = my = mz = cells;
    reps = argc > 2 ? std::atoi(argv[2]) : 6;
    if (argc > 3) {
      op = argv[3];
    }
  }
  const double cutoff = 3.5;
  const double halo = haloOf(op, cutoff);
  const auto cloud = ice(mx, my, mz);
  for (int n = 0; n < reps + 1; n++) {
    MPI_Barrier(MPI_COMM_WORLD);
    const auto wall0 = std::chrono::steady_clock::now();
    auto work = seams::domain::gatherWork(cloud, halo, MPI_COMM_WORLD);
    const auto compute0 = std::chrono::steady_clock::now();
    compute(op, work, cutoff);
    const auto compute1 = std::chrono::steady_clock::now();
    MPI_Barrier(MPI_COMM_WORLD);
    const double wall = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - wall0)
                            .count();
    const double computeMs =
        std::chrono::duration<double, std::milli>(compute1 - compute0).count();
    double comm = work.commMs;
    double haloFrac = work.haloFraction;
    double computeMax = 0.0;
    MPI_Allreduce(&comm, &comm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&haloFrac, &haloFrac, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&computeMs, &computeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0 && n > 0) {
      std::printf("%.6f %.6f %.6f %.6f\n", wall, comm, computeMax, haloFrac);
    }
  }
  MPI_Finalize();
  return 0;
}
