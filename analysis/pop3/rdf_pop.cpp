#include "pop3.h"

#include <rdf2d.hpp>

#include <omp.h>

#include <chrono>
#include <cstring>
#include <dlfcn.h>
#include <iostream>
#include <random>
#include <string>
#include <vector>

/* In-plane RDF on the 8192-atom, 70 angstrom, cutoff-6 frame. The packed
 * cutoff grid is the path sampleRDF_AA takes here. Efficiencies come from
 * the GOMP_parallel interposer. The whole-call figure includes the serial
 * grid build, so it is reported beside the region efficiencies.
 */

namespace {

using clock_mono = std::chrono::steady_clock;

long long pair_count(const std::vector<int> &hist) {
  long long sum = 0;
  for (int bin : hist)
    sum += bin;
  return sum / 2;
}

molSys::PointCloud<molSys::Point<double>, double> make_frame() {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  constexpr int n = 8192;
  constexpr double box = 70.0;
  cloud.box = {box, box, box};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  cloud.nop = n;
  cloud.pts.resize(static_cast<std::size_t>(n));
  std::mt19937 rng(1);
  std::uniform_real_distribution<double> u(0.0, box);
  for (int i = 0; i < n; i++) {
    auto &pt = cloud.pts[static_cast<std::size_t>(i)];
    pt.type = 1;
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = u(rng);
    pt.y = u(rng);
    pt.z = u(rng);
    cloud.idIndexMap[i + 1] = i;
  }
  return cloud;
}

double seconds(clock_mono::duration d) {
  return std::chrono::duration<double>(d).count();
}

} // namespace

int main() {
  auto *reset = reinterpret_cast<void (*)()>(dlsym(RTLD_DEFAULT, "pop3_reset"));
  auto *arm = reinterpret_cast<void (*)(int)>(dlsym(RTLD_DEFAULT, "pop3_arm"));
  auto *probe = reinterpret_cast<void (*)()>(dlsym(RTLD_DEFAULT, "pop3_probe"));
  auto *snap = reinterpret_cast<int (*)(pop3_view *)>(
      dlsym(RTLD_DEFAULT, "pop3_snapshot"));
  if (!reset || !arm || !probe || !snap) {
    std::cerr << "LD_PRELOAD the pop3 GOMP_parallel interposer\n";
    return 2;
  }

  const auto cloud = make_frame();
  constexpr double cutoff = 6.0;
  constexpr double binwidth = 0.1;
  constexpr int nbin = 60;
  constexpr int reps = 200;
  const int widths[] = {1, 8};

  probe();
  std::cout << "pop3 rdf2d sampleRDF_AA\n";
  std::cout << "n=8192 box=70 cutoff=6 binwidth=0.1 nbin=60 seed=1 reps=200\n";
  std::cout << "useful_time=CLOCK_THREAD_CPUTIME_ID\n";
  std::cout << "region_runtime=CLOCK_MONOTONIC around GOMP_parallel\n";
  std::cout << "spin_wait counts as useful; communication efficiency is the slept gap\n";

  double whole[2] = {0.0, 0.0};
  pop3_view views[2];
  long long pairs[2] = {0, 0};

  for (int k = 0; k < 2; k++) {
    omp_set_num_threads(widths[k]);
    auto warm = rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);
    pairs[k] = pair_count(warm);
    reset();
    arm(1);
    const auto t0 = clock_mono::now();
    for (int r = 0; r < reps; r++)
      (void)rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);
    const auto t1 = clock_mono::now();
    arm(0);
    snap(&views[k]);
    whole[k] = seconds(t1 - t0) / static_cast<double>(reps);
  }

  std::cout.setf(std::ios::fixed);
  std::cout.precision(6);
  for (int k = 0; k < 2; k++) {
    const pop3_view &v = views[k];
    std::cout << "omp_threads=" << widths[k] << "\n";
    std::cout << "  pairs=" << pairs[k] << "\n";
    std::cout << "  regions=" << v.regions << "\n";
    std::cout << "  team=" << v.threads << "\n";
    std::cout << "  whole_call_s=" << whole[k] << "\n";
    std::cout << "  region_runtime_s=" << (v.runtime_s / reps) << "\n";
    std::cout << "  useful_mean_s=" << (v.useful_mean_s / reps) << "\n";
    std::cout << "  useful_max_s=" << (v.useful_max_s / reps) << "\n";
    std::cout << "  load_balance=" << v.load_balance << "\n";
    std::cout << "  communication_efficiency=" << v.communication_efficiency
              << "\n";
    std::cout << "  parallel_efficiency=" << v.parallel_efficiency << "\n";
    for (int t = 0; t < v.threads; t++)
      std::cout << "  useful_s[" << t << "]=" << (v.useful_s[t] / reps) << "\n";
    if (v.cycles > 0) {
      std::cout << "  ipc=" << v.ipc << "\n";
      std::cout << "  instructions=" << v.instructions << "\n";
      std::cout << "  cycles=" << v.cycles << "\n";
    } else {
      std::cout << "  ipc=unavailable errno=" << v.ipc_errno << " "
                << std::strerror(v.ipc_errno) << "\n";
    }
  }
  const double scaling = whole[0] / (8.0 * whole[1]);
  const double region_scale = views[0].runtime_s / (8.0 * views[1].runtime_s);
  const double computation =
      views[0].useful_mean_s / (8.0 * views[1].useful_mean_s);
  std::cout << "region_strong_scaling_efficiency=" << region_scale << "\n";
  std::cout << "computation_scalability=" << computation << "\n";
  std::cout << "strong_scaling_efficiency=" << scaling << "\n";
  if (pairs[0] != pairs[1]) {
    std::cerr << "pair counts differ across thread counts\n";
    return 1;
  }
  return 0;
}
