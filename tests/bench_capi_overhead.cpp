/* Copyright (c) d-SEAMS developers
 * Cost of the C API boundary against the direct C++ calls it wraps.
 *
 * CHILL+ on a cutoff graph (published rule, 3.5 A) of an ice Ic supercell:
 *   capi  seams_chill_plus_n on a handle loaded with seams_system_set_atoms
 *   cpp   nneigh::orthoFullCsr + chill::chillPlusFromCsr on a Cloud
 * Both run serial and interleaved, one rep of each per round. The load
 * (set_atoms against building the Cloud) is timed apart from the analysis.
 * A trivial call (seams_system_natoms) gives the fixed cost per call.
 *
 * Usage: bench_capi_overhead [reps] [cells...]   default 11 reps, cells 8 24
 * Cells per edge k give 8 k^3 molecules: 8 -> 4096, 24 -> 110592.
 */
#include <seams.h>

#include <bop.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;
using Clock = std::chrono::steady_clock;

struct Lattice {
  std::vector<double> xyz;
  double box[9] = {};
  int n = 0;
};

Lattice ice_ic(int k) {
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  Lattice lat;
  lat.n = 8 * k * k * k;
  lat.xyz.reserve(static_cast<std::size_t>(lat.n) * 3);
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) {
      for (int l = 0; l < k; ++l) {
        for (const auto &f : fcc) {
          for (const auto &b : basis) {
            lat.xyz.push_back((i + f[0] + b[0]) * a);
            lat.xyz.push_back((j + f[1] + b[1]) * a);
            lat.xyz.push_back((l + f[2] + b[2]) * a);
          }
        }
      }
    }
  }
  const double edge = k * a;
  const double box[9] = {0, 0, 0, edge, edge, edge, 0, 0, 0};
  std::copy(box, box + 9, lat.box);
  return lat;
}

Cloud cloud_of(const Lattice &lat) {
  Cloud cloud;
  nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, lat.box[0], lat.box[1], lat.box[2],
                                lat.box[3], lat.box[4], lat.box[5], lat.box[6], lat.box[7],
                                lat.box[8]);
  cloud.nop = lat.n;
  cloud.currentFrame = 1;
  cloud.pts.resize(static_cast<std::size_t>(lat.n));
  for (int i = 0; i < lat.n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.x = lat.xyz[static_cast<std::size_t>(3 * i)];
    point.y = lat.xyz[static_cast<std::size_t>(3 * i + 1)];
    point.z = lat.xyz[static_cast<std::size_t>(3 * i + 2)];
    point.type = 1;
    point.atomID = i + 1;
    point.molID = point.atomID;
    point.iceType = molSys::atom_state_type::unclassified;
    cloud.idIndexMap.emplace(point.atomID, i);
  }
  return cloud;
}

double seconds(Clock::time_point t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t m = v.size() / 2;
  return v.size() % 2 ? v[m] : 0.5 * (v[m - 1] + v[m]);
}

double lowest(const std::vector<double> &v) { return *std::min_element(v.begin(), v.end()); }

void need(seams_status status, const char *what) {
  if (status != SEAMS_OK) {
    std::fprintf(stderr, "%s: %s (%s)\n", what, seams_status_string(status),
                 seams_last_error_message());
    std::exit(1);
  }
}

int cubic_of(const Cloud &cloud) {
  int found = 0;
  for (const auto &point : cloud.pts) {
    found += point.iceType == molSys::atom_state_type::cubic;
  }
  return found;
}

void trivial_calls(seams_system *sys) {
  const long calls = 20000000;
  size_t n = 0;
  size_t sum = 0;
  std::vector<double> per;
  for (int r = 0; r < 5; ++r) {
    const auto t0 = Clock::now();
    for (long i = 0; i < calls; ++i) {
      seams_system_natoms(sys, &n);
      sum += n;
    }
    per.push_back(seconds(t0) / static_cast<double>(calls) * 1e9);
  }
  std::printf("trivial seams_system_natoms median %.2f ns min %.2f ns per call (sum %zu)\n",
              median(per), lowest(per), sum % 7);
}

void run(int k, int reps) {
  const Lattice lat = ice_ic(k);
  const double cutoff = 3.5;
  seams_system *sys = seams_system_create();
  if (sys == nullptr) {
    std::exit(1);
  }
  std::vector<int> labels(static_cast<std::size_t>(lat.n));
  int counts9[9] = {};
  std::vector<int> offsets;
  std::vector<int> cols;
  std::vector<double> load_c, load_cpp, chill_c, chill_cpp;

  for (int r = -1; r < reps; ++r) {
    auto t0 = Clock::now();
    need(seams_system_set_atoms(sys, lat.xyz.data(), static_cast<size_t>(lat.n), lat.box, 9,
                                nullptr, nullptr),
         "set_atoms");
    const double lc = seconds(t0);

    t0 = Clock::now();
    need(seams_chill_plus_n(sys, cutoff, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                            labels.data(), labels.size(), counts9),
         "chill_plus_n");
    const double cc = seconds(t0);

    t0 = Clock::now();
    Cloud cloud = cloud_of(lat);
    const double lp = seconds(t0);

    t0 = Clock::now();
    if (!nneigh::orthoFullCsr(cloud, nullptr, cutoff, offsets, cols, false)) {
      std::fprintf(stderr, "orthoFullCsr refused the box\n");
      std::exit(1);
    }
    chill::chillPlusFromCsr(cloud, offsets, cols, false);
    const double cp = seconds(t0);

    if (counts9[SEAMS_LABEL_CUBIC] != cubic_of(cloud) || counts9[SEAMS_LABEL_CUBIC] != lat.n) {
      std::fprintf(stderr, "paths disagree: capi %d cpp %d of %d\n", counts9[SEAMS_LABEL_CUBIC],
                   cubic_of(cloud), lat.n);
      std::exit(1);
    }
    if (r >= 0) {
      load_c.push_back(lc);
      load_cpp.push_back(lp);
      chill_c.push_back(cc);
      chill_cpp.push_back(cp);
    }
  }

  const double mc = median(chill_c), mp = median(chill_cpp);
  std::printf("molecules %d reps %d\n", lat.n, reps);
  std::printf("  chill+ capi  median %.3f ms min %.3f ms\n", mc * 1e3, lowest(chill_c) * 1e3);
  std::printf("  chill+ cpp   median %.3f ms min %.3f ms\n", mp * 1e3, lowest(chill_cpp) * 1e3);
  std::printf("  chill+ overhead median %.3f ms (%+.2f%%), min-to-min %+.2f%%\n",
              (mc - mp) * 1e3, 100.0 * (mc - mp) / mp,
              100.0 * (lowest(chill_c) - lowest(chill_cpp)) / lowest(chill_cpp));
  std::printf("  load capi set_atoms median %.3f ms, cpp Cloud median %.3f ms\n",
              median(load_c) * 1e3, median(load_cpp) * 1e3);
  if (k == 8) {
    trivial_calls(sys);
  }
  seams_system_destroy(sys);
}

} // namespace

int main(int argc, char **argv) {
  const int reps = argc > 1 ? std::max(8, std::atoi(argv[1])) : 11;
  std::vector<int> cells;
  for (int i = 2; i < argc; ++i) {
    cells.push_back(std::atoi(argv[i]));
  }
  if (cells.empty()) {
    cells = {8, 24};
  }
  need(seams_set_num_threads(1), "threads");
#ifdef _OPENMP
  omp_set_num_threads(1);
#endif
  std::printf("seams %s, serial, cores visible %u\n", seams_version_string(),
              std::thread::hardware_concurrency());
  for (int k : cells) {
    run(k, reps);
  }
  return 0;
}
