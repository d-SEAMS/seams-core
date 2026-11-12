/*
** Timed driver for the four bulk analyses on one cubic-ice frame.
**
**   seams_bench [--cells N] [--traj PATH] [--frame N] [--type N]
**               [--cutoff R] [--case all|neigh|q3|q6|rings|ice]
**               [--warmup K] [--reps K]
**
** Default input is a periodic diamond (ice Ic) lattice, 8x8x8
** conventional cells, 4096 molecules, lattice constant 6.37 Angstrom,
** the same size as input/traj/mW_cubic.lammpstrj. Setup is outside
** the timed region. Each line is one full run of one analysis.
*/

#include <bop.hpp>
#include <cage.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <ring.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Clock = std::chrono::steady_clock;

void mix(std::uint64_t &h, std::uint64_t v) {
  h ^= v;
  h *= 0x100000001b3ULL;
}

void mixInt(std::uint64_t &h, int v) {
  mix(h, static_cast<std::uint64_t>(static_cast<std::uint32_t>(v)));
}

void mixDouble(std::uint64_t &h, double x) {
  std::uint64_t u = 0;
  std::memcpy(&u, &x, sizeof(u));
  mix(h, u);
}

std::uint64_t hashRows(const std::vector<std::vector<int>> &rows) {
  std::uint64_t h = 0xcbf29ce484222325ULL;
  mixInt(h, static_cast<int>(rows.size()));
  for (const auto &row : rows) {
    mixInt(h, static_cast<int>(row.size()));
    for (const int v : row) {
      mixInt(h, v);
    }
  }
  return h;
}

std::uint64_t hashDoubles(const std::vector<double> &xs) {
  std::uint64_t h = 0xcbf29ce484222325ULL;
  mixInt(h, static_cast<int>(xs.size()));
  for (const double x : xs) {
    mixDouble(h, x);
  }
  return h;
}

Cloud diamondIce(int cells) {
  // Conventional cubic diamond cell. Nearest neighbour is a*sqrt(3)/4.
  constexpr double a = 6.37;
  constexpr double basis[8][3] = {
      {0.00, 0.00, 0.00}, {0.50, 0.50, 0.00}, {0.50, 0.00, 0.50},
      {0.00, 0.50, 0.50}, {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
      {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75},
  };
  Cloud cloud;
  const double L = a * static_cast<double>(cells);
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  const int n = cells * cells * cells * 8;
  cloud.pts.reserve(static_cast<std::size_t>(n));
  int id = 0;
  for (int ix = 0; ix < cells; ++ix) {
    for (int iy = 0; iy < cells; ++iy) {
      for (int iz = 0; iz < cells; ++iz) {
        for (const auto &b : basis) {
          molSys::Point<double> pt;
          pt.type = 1;
          pt.molID = id;
          pt.atomID = id;
          pt.x = (static_cast<double>(ix) + b[0]) * a;
          pt.y = (static_cast<double>(iy) + b[1]) * a;
          pt.z = (static_cast<double>(iz) + b[2]) * a;
          cloud.pts.push_back(pt);
          cloud.idIndexMap[id] = id;
          ++id;
        }
      }
    }
  }
  cloud.nop = id;
  return cloud;
}

struct Prepared {
  Cloud cloud;
  std::vector<std::vector<int>> byId;
  std::vector<std::vector<int>> byIndex;
  std::vector<std::vector<int>> rings;
  std::vector<std::vector<int>> six;
};

void buildNeighbours(Prepared &prep, double cutoff) {
  prep.byId = nneigh::neighListO(cutoff, prep.cloud, 1);
  prep.byIndex = nneigh::neighbourListByIndex(prep.cloud, prep.byId);
}

void buildRings(Prepared &prep) {
  prep.rings = primitive::ringNetwork(prep.byIndex, 6);
  prep.six = ring::getSingleRingSize(prep.rings, 6);
}

std::uint64_t runNeigh(const Cloud &cloud, double cutoff) {
  const auto rows = nneigh::neighListO(cutoff, cloud, 1);
  return hashRows(rows);
}

std::uint64_t runQl(const Prepared &prep, int orderL) {
  const auto ql = chill::steinhardtQl(prep.cloud, prep.byId, orderL);
  std::uint64_t h = 0xcbf29ce484222325ULL;
  mix(h, hashDoubles(ql.ql));
  mix(h, hashDoubles(ql.qlBar));
  return h;
}

std::uint64_t runRings(const Prepared &prep) {
  const auto rings = primitive::ringNetwork(prep.byIndex, 6);
  return hashRows(rings);
}

std::uint64_t runIce(const Prepared &prep) {
  std::vector<ring::strucType> ringType(prep.six.size(),
                                        ring::strucType::unclassified);
  std::vector<cage::Cage> cages;
  const auto index = ring::buildRingSearchIndex(prep.six, prep.cloud);
  const auto hc = ring::findHC(prep.six, ringType, prep.byIndex, cages, index);
  const auto ddc =
      ring::findDDC(prep.six, ringType, hc, cages, index);
  std::uint64_t h = 0xcbf29ce484222325ULL;
  mixInt(h, static_cast<int>(hc.size()));
  for (const int id : hc) {
    mixInt(h, id);
  }
  mixInt(h, static_cast<int>(ddc.size()));
  for (const int id : ddc) {
    mixInt(h, id);
  }
  mixInt(h, static_cast<int>(cages.size()));
  for (const auto &cage : cages) {
    mixInt(h, static_cast<int>(cage.type));
    for (const int ringId : cage.rings) {
      mixInt(h, ringId);
    }
  }
  for (const auto kind : ringType) {
    mixInt(h, static_cast<int>(kind));
  }
  return h;
}

double secondsSince(Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

template <typename Fn>
double timeCall(Fn &&fn, std::uint64_t &digest) {
  asm volatile("" ::: "memory");
  const auto t0 = Clock::now();
  digest = fn();
  asm volatile("" : "+m"(digest)::"memory");
  return secondsSince(t0);
}

void report(const char *name, int rep, double seconds, std::uint64_t digest) {
  std::printf("case %s rep %d seconds %.9f checksum %016llx\n", name, rep,
              seconds, static_cast<unsigned long long>(digest));
}

struct Options {
  int cells = 8;
  int frame = 1;
  int type = 1;
  int warmup = 1;
  int reps = 1;
  double cutoff = 3.5;
  std::string traj;
  std::string which = "all";
};

bool parseInt(const char *text, int &out) {
  char *end = nullptr;
  const long v = std::strtol(text, &end, 10);
  if (end == text || *end != '\0') {
    return false;
  }
  out = static_cast<int>(v);
  return true;
}

bool parseDouble(const char *text, double &out) {
  char *end = nullptr;
  const double v = std::strtod(text, &end);
  if (end == text || *end != '\0') {
    return false;
  }
  out = v;
  return true;
}

void usage() {
  std::fprintf(stderr,
               "usage: seams_bench [--cells N] [--traj PATH] [--frame N]\n"
               "                   [--type N] [--cutoff R]\n"
               "                   [--case all|neigh|q3|q6|rings|ice]\n"
               "                   [--warmup K] [--reps K]\n");
}

} // namespace

int main(int argc, char **argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need = [&](const char *flag) -> const char * {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "missing value for %s\n", flag);
        std::exit(2);
      }
      return argv[++i];
    };
    if (arg == "--cells") {
      if (!parseInt(need("--cells"), opt.cells) || opt.cells < 1) {
        std::fprintf(stderr, "bad --cells\n");
        return 2;
      }
    } else if (arg == "--traj") {
      opt.traj = need("--traj");
    } else if (arg == "--frame") {
      if (!parseInt(need("--frame"), opt.frame) || opt.frame < 1) {
        std::fprintf(stderr, "bad --frame\n");
        return 2;
      }
    } else if (arg == "--type") {
      if (!parseInt(need("--type"), opt.type)) {
        std::fprintf(stderr, "bad --type\n");
        return 2;
      }
    } else if (arg == "--cutoff") {
      if (!parseDouble(need("--cutoff"), opt.cutoff) || !(opt.cutoff > 0.0)) {
        std::fprintf(stderr, "bad --cutoff\n");
        return 2;
      }
    } else if (arg == "--case") {
      opt.which = need("--case");
    } else if (arg == "--warmup") {
      if (!parseInt(need("--warmup"), opt.warmup) || opt.warmup < 0) {
        std::fprintf(stderr, "bad --warmup\n");
        return 2;
      }
    } else if (arg == "--reps") {
      if (!parseInt(need("--reps"), opt.reps) || opt.reps < 1) {
        std::fprintf(stderr, "bad --reps\n");
        return 2;
      }
    } else if (arg == "--help") {
      usage();
      return 0;
    } else {
      std::fprintf(stderr, "unknown argument %s\n", arg.c_str());
      usage();
      return 2;
    }
  }

#ifdef SEAMS_HAS_OPENMP
  omp_set_dynamic(0);
  if (const char *env = std::getenv("SEAMS_BENCH_THREADS")) {
    omp_set_num_threads(std::max(1, std::atoi(env)));
  }
#endif

  Cloud cloud;
  if (!opt.traj.empty()) {
    cloud = sinp::readLammpsTrjO(opt.traj, opt.frame, cloud, opt.type);
  } else {
    cloud = diamondIce(opt.cells);
  }
  if (cloud.nop <= 0) {
    std::fprintf(stderr, "empty configuration\n");
    return 1;
  }

  const bool all = opt.which == "all";
  const bool wantNeigh = all || opt.which == "neigh";
  const bool wantQ3 = all || opt.which == "q3";
  const bool wantQ6 = all || opt.which == "q6";
  const bool wantRings = all || opt.which == "rings";
  const bool wantIce = all || opt.which == "ice";
  if (!wantNeigh && !wantQ3 && !wantQ6 && !wantRings && !wantIce) {
    std::fprintf(stderr, "unknown case %s\n", opt.which.c_str());
    usage();
    return 2;
  }

  // Neighbour rows and rings are inputs of the later analyses.
  // Building them here keeps each timed case to its own work.
  Prepared prep;
  prep.cloud = std::move(cloud);
  if (wantQ3 || wantQ6 || wantRings || wantIce) {
    buildNeighbours(prep, opt.cutoff);
  }
  if (wantIce) {
    buildRings(prep);
  }

  std::printf("nop %d six %zu rings %zu\n", prep.cloud.nop, prep.six.size(),
              prep.rings.size());
  std::fflush(stdout);

  auto once = [&](const char *name, auto &&fn) {
    std::uint64_t digest = 0;
    for (int i = 0; i < opt.warmup; ++i) {
      digest = fn();
    }
    for (int rep = 0; rep < opt.reps; ++rep) {
      const double seconds = timeCall(fn, digest);
      report(name, rep, seconds, digest);
      std::fflush(stdout);
    }
    (void)digest;
  };

  if (wantNeigh) {
    once("neigh", [&]() { return runNeigh(prep.cloud, opt.cutoff); });
  }
  if (wantQ3) {
    once("q3", [&]() { return runQl(prep, 3); });
  }
  if (wantQ6) {
    once("q6", [&]() { return runQl(prep, 6); });
  }
  if (wantRings) {
    once("rings", [&]() { return runRings(prep); });
  }
  if (wantIce) {
    once("ice", [&]() { return runIce(prep); });
  }
  return 0;
}
