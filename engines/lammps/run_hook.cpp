// Copyright (c) d-SEAMS developers
// Drives the installed LAMMPS library. fix external calls the C API each step.

#include "library.h"
#include "seams.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdint>

struct Hook {
  void *lmp = nullptr;
  int calls = 0;
  int ice = 0;
  int atoms = 0;
  int stride = 1;
  double seconds = 0.0;
  std::vector<int64_t> steps;
};

void callback(void *ptr, int64_t step, int nlocal, int *, double **x, double **) {
  auto *hook = static_cast<Hook *>(ptr);
  std::vector<double> xyz(static_cast<std::size_t>(nlocal) * 3);
  for (int i = 0; i < nlocal; ++i) {
    xyz[static_cast<std::size_t>(3 * i)] = x[i][0];
    xyz[static_cast<std::size_t>(3 * i + 1)] = x[i][1];
    xyz[static_cast<std::size_t>(3 * i + 2)] = x[i][2];
  }
  double lo[3], hi[3], xy, yz, xz;
  int pflags[3], boxflag;
  lammps_extract_box(hook->lmp, lo, hi, &xy, &yz, &xz, pflags, &boxflag);
  const double box9[9] = {lo[0], lo[1], lo[2], hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2],
                          xy, xz, yz};
  std::vector<int> labels(static_cast<std::size_t>(nlocal));
  const auto t0 = std::chrono::steady_clock::now();
  seams_system *sys = seams_system_create();
  const auto set = seams_system_set_frame(sys, xyz.data(), nlocal, box9, nullptr, nullptr);
  int counts[9] = {};
  const auto code = seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                                     labels.data(), counts);
  seams_system_destroy(sys);
  const auto t1 = std::chrono::steady_clock::now();
  hook->seconds += std::chrono::duration<double>(t1 - t0).count();
  hook->calls += 1;
  hook->atoms = nlocal;
  hook->steps.push_back(step);
  if (set == SEAMS_OK && code == SEAMS_OK) {
    hook->ice += counts[SEAMS_LABEL_CUBIC];
  } else {
    hook->ice = -1;
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: run_hook in.dseams\n");
    return 2;
  }
  Hook hook;
  int narg = 1;
  char *extra[] = {argv[0]};
  hook.lmp = lammps_open_no_mpi(narg, extra, nullptr);
  if (hook.lmp == nullptr) {
    std::fprintf(stderr, "lammps_open failed\n");
    return 1;
  }
  std::ifstream in(argv[1]);
  if (!in) {
    std::fprintf(stderr, "cannot read %s\n", argv[1]);
    return 1;
  }
  std::string line;
  std::string run_line = "run 8";
  while (std::getline(in, line)) {
    const auto start = line.find_first_not_of(" \t");
    if (start == std::string::npos || line[start] == '#') {
      continue;
    }
    if (line.compare(start, 3, "run") == 0 &&
        (line.size() == start + 3 || line[start + 3] == ' ' || line[start + 3] == '\t')) {
      run_line = line.substr(start);
      continue;
    }
    const auto callback_at = line.find("pf/callback");
    if (callback_at != std::string::npos) {
      int ncall = 1;
      int napply = 1;
      if (std::sscanf(line.c_str() + callback_at, "pf/callback %d %d", &ncall, &napply) == 2 &&
          ncall > 0) {
        hook.stride = ncall;
      }
    }
    lammps_command(hook.lmp, line.c_str());
  }
  lammps_set_fix_external_callback(hook.lmp, "ext", callback, &hook);
  lammps_command(hook.lmp, run_line.c_str());
  lammps_close(hook.lmp);
  bool spaced = hook.steps.size() >= 2;
  for (std::size_t i = 1; i < hook.steps.size(); ++i) {
    if (hook.steps[i] - hook.steps[i - 1] != hook.stride) {
      spaced = false;
    }
  }
  if (!spaced || hook.ice != hook.calls * hook.atoms) {
    std::fprintf(stderr, "lammps hook calls=%d ice=%d atoms=%d stride=%d\n", hook.calls, hook.ice,
                 hook.atoms, hook.stride);
    return 1;
  }
  std::printf("lammps calls %d stride %d ice %d overhead_s %.6e\n", hook.calls, hook.stride,
              hook.ice, hook.seconds / hook.calls);
  return 0;
}
