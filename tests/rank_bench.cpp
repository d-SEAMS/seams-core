// Time one distributed operation. Rank 0 prints, per timed call,
// wall milliseconds, ghost-exchange milliseconds, compute milliseconds,
// and the largest rank's halo fraction. The first call is the warm-up.

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <domain.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include "slab_frames.hpp"

#include <mpi.h>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <thread>

namespace {

void compute(const std::string &op,
             molSys::PointCloud<molSys::Point<double>, double> &full, double cutoff) {
  if (op == "rings") {
    auto rings = seams::domain::distributedRings(full, cutoff, 6, MPI_COMM_WORLD);
    if (rings.empty() && full.nop > 1000) {
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0) {
        std::fprintf(stderr, "no rings\n");
      }
    }
    return;
  }
  if (op == "chill") {
    chill::distributedChillPlus(full, cutoff, chill::chillPlusRule(), true, MPI_COMM_WORLD);
    return;
  }
  if (op == "q6") {
    auto q = chill::distributedSteinhardt(full, cutoff, 6, MPI_COMM_WORLD);
    if (q.ql.empty() && full.nop > 0) {
      std::fprintf(stderr, "empty q6\n");
    }
    return;
  }
  auto &work = seams::domain::gatherWork(full, cutoff, MPI_COMM_WORLD, seams::domain::kWorkNeigh);
  bool single = true;
  for (const char owned : work.share.owned) {
    if (!owned) {
      single = false;
      break;
    }
  }
  const std::vector<std::vector<int>> *rows = nullptr;
  static std::vector<std::vector<int>> ownedRows;
  const bool window =
      nneigh::blockWindowRows(work.cloud, cutoff, single ? nullptr : &work.share.owned, ownedRows);
  if (window) {
    rows = &ownedRows;
  } else if (single || work.cloud.nop < 40000) {
    if (!nneigh::refillNeighListO(cutoff, work.cloud, 1, ownedRows)) {
      ownedRows = nneigh::neighListO(cutoff, work.cloud, 1);
    }
    rows = &ownedRows;
  } else {
    rows = &nneigh::reusedShellList(cutoff, work.cloud, &work.share.owned);
  }
  if (rows->size() != work.share.local.size() && work.cloud.nop > 0) {
    std::fprintf(stderr, "neighbour rows %zu\n", rows->size());
  }
  std::vector<double> effort(rows->size(), 0.0);
  for (std::size_t k = 0; k < rows->size() && k < work.share.owned.size(); k++) {
    if (work.share.owned[k]) {
      effort[k] = static_cast<double>((*rows)[k].size());
    }
  }
  seams::domain::captureLocalWork(seams::domain::kWorkNeigh, effort.data(),
                                  static_cast<int>(effort.size()), MPI_COMM_WORLD);
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
  bool slab = false;
  std::string kind = "ice";
  int which = 110;
  bool fromFile = false;
  bool frames = false;
  std::string framePath;
  if (argc >= 2 && std::string(argv[1]) == "frames") {
    frames = true;
    mx = argc > 2 ? std::atoi(argv[2]) : 48;
    my = argc > 3 ? std::atoi(argv[3]) : 48;
    mz = argc > 4 ? std::atoi(argv[4]) : 48;
    reps = argc > 5 ? std::atoi(argv[5]) : 8;
    if (argc > 6) {
      op = argv[6];
    }
  } else if (argc >= 2 && std::string(argv[1]) == "file") {
    fromFile = true;
    framePath = argc > 2 ? argv[2] : "";
    reps = argc > 3 ? std::atoi(argv[3]) : 6;
    if (argc > 4) {
      op = argv[4];
    }
  } else if (argc >= 2 && std::string(argv[1]) == "slab") {
    slab = true;
    kind = argc > 2 ? argv[2] : "ice";
    which = argc > 3 ? std::atoi(argv[3]) : 110;
    reps = argc > 4 ? std::atoi(argv[4]) : 6;
    if (argc > 5) {
      op = argv[5];
    }
  } else if (argc >= 6) {
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
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  if (fromFile) {
    cloud = sinp::readLammpsTrj(framePath, 1, cloud);
  } else if (!slab) {
    cloud = seams::frames::iceLattice(mx, my, mz);
  } else if (kind == "hydrate") {
    cloud = which >= 800 ? seams::frames::hydrateLiquidSlab(40, 24, 20)
                         : seams::frames::hydrateLiquidSlab(24, 24, 4);
  } else {
    cloud = seams::frames::iceWaterSlab(which >= 800 ? 48 : 24);
  }
  if (rank == 0) {
    std::fprintf(stderr, "frame %s nop %d cores %u\n",
                 frames ? "thermal" : (fromFile ? framePath.c_str() : (slab ? kind.c_str() : "ice")),
                 cloud.nop, std::thread::hardware_concurrency());
  }
  if (frames) {
    const auto base = cloud;
    const double amp = 0.12;
    auto place = [&](std::uint32_t seed) {
      cloud = base;
      std::mt19937 rng(seed);
      std::normal_distribution<double> shift(0.0, amp);
      const double box[3] = {cloud.box[0], cloud.box[1], cloud.box[2]};
      for (auto &p : cloud.pts) {
        const double moved[3] = {p.x + shift(rng), p.y + shift(rng), p.z + shift(rng)};
        const double *m = moved;
        double *out[3] = {&p.x, &p.y, &p.z};
        for (int d = 0; d < 3; d++) {
          double s = std::fmod(m[d], box[d]);
          if (s < 0.0) {
            s += box[d];
          }
          *out[d] = s;
        }
      }
    };
    auto once = [&](const char *tag) {
      seams::domain::setCutMode(2);
      MPI_Barrier(MPI_COMM_WORLD);
      const auto wall0 = std::chrono::steady_clock::now();
      const auto compute0 = std::chrono::steady_clock::now();
      const auto &rings = seams::domain::distributedRings(cloud, cutoff, 6, MPI_COMM_WORLD);
      const auto compute1 = std::chrono::steady_clock::now();
      MPI_Barrier(MPI_COMM_WORLD);
      const double wall = std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - wall0)
                              .count();
      const double computeMs =
          std::chrono::duration<double, std::milli>(compute1 - compute0).count();
      const auto &work = seams::domain::gatheredWork();
      double haloFrac = work.haloFraction;
      double computeMax = 0.0;
      double computeSum = 0.0;
      MPI_Allreduce(&haloFrac, &haloFrac, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&computeMs, &computeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&computeMs, &computeSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      long long localCount[7] = {0, 0, 0, 0, 0, 0, 0};
      for (const auto &ring : rings) {
        const int sz = static_cast<int>(ring.size());
        if (sz >= 0 && sz <= 6) {
          localCount[sz]++;
        }
      }
      long long count[7] = {0, 0, 0, 0, 0, 0, 0};
      const auto red0 = std::chrono::steady_clock::now();
      MPI_Allreduce(localCount, count, 7, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      const double reduceMs =
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - red0)
              .count();
      int nRanks = 1;
      MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
      const double mean = nRanks > 0 ? computeSum / nRanks : 0.0;
      const double imbalance = mean > 0.0 ? computeMax / mean : 1.0;
      const auto phase = seams::domain::lastRingPhases();
      double part[4] = {phase.rebuildMs, phase.searchMs, phase.primitiveMs,
                        phase.reduceMs + reduceMs};
      MPI_Allreduce(MPI_IN_PLACE, part, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (rank == 0) {
        std::printf("%s %.6f %.6f %.6f %lld %lld %lld %lld %.3f %.3f %.3f %.3f\n", tag, wall,
                    haloFrac, imbalance, count[3], count[4], count[5], count[6], part[0], part[1],
                    part[2], part[3]);
      }
    };
    auto timeOp = [&](const std::string &which, const char *tag) {
      seams::domain::setCutMode(2);
      MPI_Barrier(MPI_COMM_WORLD);
      const auto wall0 = std::chrono::steady_clock::now();
      compute(which, cloud, cutoff);
      MPI_Barrier(MPI_COMM_WORLD);
      const double wall = std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - wall0)
                              .count();
      const auto &work = seams::domain::gatheredWork();
      double haloFrac = work.haloFraction;
      MPI_Allreduce(&haloFrac, &haloFrac, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      chill::AnalysisPhases phase;
      if (which == "chill") {
        phase = chill::lastChillPhases();
      } else if (which == "q6") {
        phase = chill::lastSteinhardtPhases();
      }
      double part[4] = {phase.rebuildMs, phase.computeMs, phase.exchangeMs, phase.reduceMs};
      MPI_Allreduce(MPI_IN_PLACE, part, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (rank == 0) {
        std::printf("%s %s %.6f %.6f %.3f %.3f %.3f %.3f\n", which.c_str(), tag, wall, haloFrac,
                    part[0], part[1], part[2], part[3]);
      }
    };
    primitive::setPhaseTiming(true);
    const bool rings = op == "rings" || op == "all";
    const bool others = op == "all";
    if (rings) {
    place(1);
    // The first call records per-source time and recuts. The second search
    // is that new cut, and it is not one of the timed repeats.
    once("warm");
    once("warm");
    for (int n = 0; n < reps; n++) {
      once("same");
    }
    for (int n = 0; n < reps; n++) {
      place(static_cast<std::uint32_t>(100 + n));
      once("distinct");
    }
    }
    if (others) {
      for (const char *which : {"neigh", "chill", "q6"}) {
        place(1);
        timeOp(which, "warm");
        for (int n = 0; n < reps; n++) {
          place(static_cast<std::uint32_t>(100 + n));
          timeOp(which, "distinct");
        }
      }
    }
    if (!rings && !others) {
      place(1);
      timeOp(op, "warm");
      for (int n = 0; n < reps; n++) {
        place(static_cast<std::uint32_t>(100 + n));
        timeOp(op, "distinct");
      }
    }
    MPI_Finalize();
    return 0;
  }
  auto timeCut = [&](int mode, const char *tag, bool print) {
    seams::domain::setCutMode(mode);
    MPI_Barrier(MPI_COMM_WORLD);
    const auto wall0 = std::chrono::steady_clock::now();
    const auto compute0 = std::chrono::steady_clock::now();
    compute(op, cloud, cutoff);
    const auto compute1 = std::chrono::steady_clock::now();
    MPI_Barrier(MPI_COMM_WORLD);
    const double wall = std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - wall0)
                            .count();
    const double computeMs =
        std::chrono::duration<double, std::milli>(compute1 - compute0).count();
    const auto &work = seams::domain::gatheredWork();
    double comm = work.commMs + seams::domain::takeValueCommMs();
    double haloFrac = work.haloFraction;
    double computeMax = 0.0;
    double computeSum = 0.0;
    MPI_Allreduce(&comm, &comm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&haloFrac, &haloFrac, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&computeMs, &computeMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&computeMs, &computeSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    int nRanks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    const double mean = nRanks > 0 ? computeSum / nRanks : 0.0;
    const double imbalance = mean > 0.0 ? computeMax / mean : 1.0;
    if (rank == 0 && print) {
      std::printf("%s %.6f %.6f %.6f %.6f %.6f\n", tag, wall, comm, computeMax, haloFrac,
                  imbalance);
    }
  };
  // The scaling session times the linkcell bisection only.
  if (std::getenv("SEAMS_SCALE") != nullptr) {
    timeCut(2, "cell", false);
    for (int n = 0; n < reps; n++) {
      timeCut(2, "cell", true);
    }
    MPI_Finalize();
    return 0;
  }
  // One warm-up of each cut, then block and cell in the same loop.
  timeCut(1, "block", false);
  timeCut(2, "cell", false);
  for (int n = 0; n < reps; n++) {
    timeCut(1, "block", true);
    timeCut(2, "cell", true);
  }
  MPI_Finalize();
  return 0;
}
