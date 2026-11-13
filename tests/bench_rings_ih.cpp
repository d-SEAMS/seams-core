// Time seams_rings_find on ideal ice Ih.
//
//   bench_rings_ih --session [--reps 8] [--warmup 1]
//   bench_rings_ih --molecules 4096 --threads 1 [--reps 1] [--warmup 0]
//
// The report names the core count. It does not name the host.

#include "ice_ih.hpp"
#include "seams.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

int coresOnline() {
  const long n = sysconf(_SC_NPROCESSORS_ONLN);
  return n > 0 ? static_cast<int>(n) : 1;
}

double millis(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

double medianOf(std::vector<double> samples) {
  if (samples.empty()) {
    return 0;
  }
  std::sort(samples.begin(), samples.end());
  const std::size_t n = samples.size();
  if (n % 2 == 1) {
    return samples[n / 2];
  }
  return 0.5 * (samples[n / 2 - 1] + samples[n / 2]);
}

struct Symmetry {
  int count = 0;
  int basal = 0;
  int prismatic = 0;
  int hexagon = 0;
  int prism = 0;
};

Symmetry labelsOf(seams_rings *rings) {
  Symmetry out;
  out.count = seams_rings_count(rings);
  std::vector<int> sym(static_cast<std::size_t>(out.count));
  if (seams_rings_labels(rings, nullptr, sym.data()) != SEAMS_OK) {
    return out;
  }
  for (int code : sym) {
    if (code == SEAMS_RING_HC_BASAL) {
      ++out.basal;
    } else if (code == SEAMS_RING_HC_PRISMATIC) {
      ++out.prismatic;
    } else if (code == SEAMS_RING_HEXAGON) {
      ++out.hexagon;
    } else if (code == SEAMS_RING_PRISM) {
      ++out.prism;
    }
  }
  return out;
}

int repsFor(int molecules) {
  int side = 0;
  for (int s = 1; s <= 64; ++s) {
    if (8 * s * s * s == molecules) {
      side = s;
      break;
    }
  }
  return side;
}

int timeOne(int molecules, int threads, int reps, int warmup, bool checkOracle) {
  const int side = repsFor(molecules);
  if (side == 0) {
    std::fprintf(stderr, "molecules %d is not 8*s^3\n", molecules);
    return 2;
  }
  if (seams_set_num_threads(threads) != SEAMS_OK) {
    std::fprintf(stderr, "threads: %s\n", seams_last_error());
    return 1;
  }
  const IceIhCell cell = iceIhOrtho(side, side, side);
  seams_system *sys = seams_system_create();
  if (sys == nullptr ||
      seams_system_set_atoms(sys, cell.xyz.data(), static_cast<std::size_t>(cell.n), cell.box9,
                             9, nullptr, nullptr) != SEAMS_OK) {
    std::fprintf(stderr, "frame: %s\n", seams_last_error());
    seams_system_destroy(sys);
    return 1;
  }
  Symmetry last;
  std::vector<double> samples;
  samples.reserve(static_cast<std::size_t>(reps));
  for (int i = 0; i < warmup + reps; ++i) {
    const auto t0 = Clock::now();
    seams_rings *rings = nullptr;
    const auto st = seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings);
    const auto t1 = Clock::now();
    if (st != SEAMS_OK) {
      std::fprintf(stderr, "rings: %s\n", seams_last_error());
      seams_system_destroy(sys);
      return 1;
    }
    last = labelsOf(rings);
    seams_rings_destroy(rings);
    if (i >= warmup) {
      samples.push_back(millis(t0, t1));
    }
  }
  seams_system_destroy(sys);
  if (checkOracle) {
    const bool ok4096 = molecules == 4096 && last.count == 8192 && last.basal == 2040 &&
                        last.prismatic == 5484 && last.hexagon == 668 && last.prism == 0;
    const bool ok32768 = molecules == 32768 && last.count == 65536 && last.basal == 16368 &&
                         last.prismatic == 46308 && last.hexagon == 2860 && last.prism == 0;
    const bool okOther = molecules != 4096 && molecules != 32768 && last.count == 2 * molecules &&
                         last.basal + last.prismatic + last.hexagon == last.count &&
                         last.prism == 0;
    if (!ok4096 && !ok32768 && !okOther) {
      std::fprintf(stderr, "census mismatch molecules %d count %d basal %d prismatic %d hexagon %d\n",
                   molecules, last.count, last.basal, last.prismatic, last.hexagon);
      return 3;
    }
  }
  std::printf("molecules %d threads %d reps %d median_ms %.3f count %d "
              "hc-basal %d hc-prismatic %d hexagon %d prism %d\n",
              molecules, threads, reps, medianOf(samples), last.count, last.basal, last.prismatic,
              last.hexagon, last.prism);
  for (double sample : samples) {
    std::printf("sample_ms %.3f\n", sample);
  }
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  bool session = false;
  int molecules = 4096;
  int threads = 1;
  int reps = 8;
  int warmup = 1;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need = [&](const char *flag) {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "missing value for %s\n", flag);
        std::exit(2);
      }
      return std::atoi(argv[++i]);
    };
    if (arg == "--session") {
      session = true;
    } else if (arg == "--molecules") {
      molecules = need("--molecules");
    } else if (arg == "--threads") {
      threads = need("--threads");
    } else if (arg == "--reps") {
      reps = need("--reps");
    } else if (arg == "--warmup") {
      warmup = need("--warmup");
    } else {
      std::fprintf(stderr, "unknown argument %s\n", arg.c_str());
      return 2;
    }
  }
  if (reps < 1) {
    std::fprintf(stderr, "reps must be positive\n");
    return 2;
  }
  std::printf("cores %d\n", coresOnline());
  if (!session) {
    return timeOne(molecules, threads, reps, warmup, true);
  }
  if (reps < 8) {
    std::fprintf(stderr, "a session records at least 8 repetitions\n");
    return 2;
  }
  const int sizes[] = {4096, 32768, 110592};
  const int widths[] = {1, 8};
  for (int size : sizes) {
    for (int width : widths) {
      const int rc = timeOne(size, width, reps, warmup, true);
      if (rc != 0) {
        return rc;
      }
    }
  }
  return 0;
}
