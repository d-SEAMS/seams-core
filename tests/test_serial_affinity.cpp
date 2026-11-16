// Runs C API calls under the cpu mask and thread width the driver
// (tests/test_serial_affinity.py) launched it with, and checks that the
// library keeps both. A sampler thread lists /proc/self/task during each
// call: the threads it sees, the cpus each may run on and the cpu each last
// ran on. After the call the process mask must still be the launch mask.
//
//   test_serial_affinity --cpus 0 --width 1 [--set-threads 1]
//
// --cpus is the mask taskset gave the process. --width is the thread count
// the process asked for; at width 1 a call may use the calling thread
// alone, and at any width every thread stays on --cpus. --set-threads
// calls seams_set_num_threads first.

#include <seams.h>

#include <dirent.h>
#include <pthread.h>
#include <sched.h>
#include <sys/syscall.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

std::set<int> parseCpus(const std::string &text) {
  std::set<int> cpus;
  std::stringstream in(text);
  std::string part;
  while (std::getline(in, part, ',')) {
    const auto dash = part.find('-');
    if (dash == std::string::npos) {
      cpus.insert(std::stoi(part));
    } else {
      for (int c = std::stoi(part.substr(0, dash)); c <= std::stoi(part.substr(dash + 1)); ++c) {
        cpus.insert(c);
      }
    }
  }
  return cpus;
}

std::string showCpus(const std::set<int> &cpus) {
  std::string out;
  for (const int c : cpus) {
    out += (out.empty() ? "" : ",") + std::to_string(c);
  }
  return out.empty() ? "none" : out;
}

std::set<int> maskOf(pid_t tid) {
  cpu_set_t mask;
  CPU_ZERO(&mask);
  std::set<int> cpus;
  if (sched_getaffinity(tid, sizeof(mask), &mask) != 0) {
    return cpus;
  }
  for (int c = 0; c < CPU_SETSIZE; ++c) {
    if (CPU_ISSET(c, &mask)) {
      cpus.insert(c);
    }
  }
  return cpus;
}

std::vector<pid_t> tasks() {
  std::vector<pid_t> ids;
  DIR *dir = opendir("/proc/self/task");
  if (dir == nullptr) {
    return ids;
  }
  while (const dirent *entry = readdir(dir)) {
    if (entry->d_name[0] != '.') {
      ids.push_back(static_cast<pid_t>(std::atoi(entry->d_name)));
    }
  }
  closedir(dir);
  return ids;
}

// Field 39 of /proc/<pid>/task/<tid>/stat, the cpu the thread last ran on.
int lastCpu(pid_t tid) {
  std::ifstream in("/proc/self/task/" + std::to_string(tid) + "/stat");
  std::string line;
  if (!std::getline(in, line)) {
    return -1;
  }
  const auto close = line.rfind(')');
  if (close == std::string::npos) {
    return -1;
  }
  std::stringstream rest(line.substr(close + 2));
  std::string field;
  for (int k = 3; k <= 39; ++k) {
    if (!(rest >> field)) {
      return -1;
    }
  }
  return std::atoi(field.c_str());
}

struct Seen {
  int peakTasks = 0;
  std::set<int> masks;
  std::set<int> ran;
};

// Lists the threads every 100 microseconds until stopped. The sampler
// itself is left out of the count.
Seen watch(const std::function<void()> &call) {
  Seen seen;
  std::atomic<bool> done{false};
  std::atomic<pid_t> self{0};
  std::thread sampler([&] {
    self = static_cast<pid_t>(syscall(SYS_gettid));
    while (!done.load()) {
      int count = 0;
      for (const pid_t tid : tasks()) {
        if (tid == self.load()) {
          continue;
        }
        ++count;
        for (const int c : maskOf(tid)) {
          seen.masks.insert(c);
        }
        const int cpu = lastCpu(tid);
        if (cpu >= 0) {
          seen.ran.insert(cpu);
        }
      }
      seen.peakTasks = std::max(seen.peakTasks, count);
      std::this_thread::sleep_for(std::chrono::microseconds(100));
    }
  });
  while (self.load() == 0) {
    std::this_thread::yield();
  }
  call();
  done = true;
  sampler.join();
  // pthread_join returns when the kernel clears the thread id, before the
  // task leaves /proc/self/task, and the next call's sampler would count
  // this one as a library thread.
  const std::string gone = "/proc/self/task/" + std::to_string(self.load());
  while (access(gone.c_str(), F_OK) == 0) {
    std::this_thread::sleep_for(std::chrono::microseconds(50));
  }
  return seen;
}

struct Frame {
  std::vector<double> xyz;
  double box[9] = {};
  int n = 0;
};

// Ice Ih, eight molecules per orthorhombic cell.
Frame iceIh(int ncell) {
  const double bond = 2.76;
  const double a = bond * std::sqrt(8.0 / 3.0);
  const double c = a * std::sqrt(8.0 / 3.0);
  const double b = a * std::sqrt(3.0);
  const double basis[4][3] = {{1.0 / 3.0, 2.0 / 3.0, 0.0},
                              {2.0 / 3.0, 1.0 / 3.0, 0.5},
                              {1.0 / 3.0, 2.0 / 3.0, 3.0 / 8.0},
                              {2.0 / 3.0, 1.0 / 3.0, 7.0 / 8.0}};
  Frame ice;
  const double lx = ncell * a;
  const double ly = ncell * b;
  const double lz = ncell * c;
  ice.box[3] = lx;
  ice.box[4] = ly;
  ice.box[5] = lz;
  ice.n = ncell * ncell * ncell * 8;
  ice.xyz.resize(static_cast<std::size_t>(ice.n) * 3);
  int atom = 0;
  for (int ix = 0; ix < ncell; ++ix) {
    for (int iy = 0; iy < ncell; ++iy) {
      for (int iz = 0; iz < ncell; ++iz) {
        for (const double shift : {0.0, 0.5}) {
          for (const auto &site : basis) {
            double x = std::fmod(shift * a + site[0] * a - site[1] * a / 2.0 + ix * a, lx);
            double y = std::fmod(shift * b + site[1] * b / 2.0 + iy * b, ly);
            const double z = (site[2] + iz) * c;
            x += x < 0.0 ? lx : 0.0;
            y += y < 0.0 ? ly : 0.0;
            ice.xyz[static_cast<std::size_t>(3 * atom)] = x;
            ice.xyz[static_cast<std::size_t>(3 * atom + 1)] = y;
            ice.xyz[static_cast<std::size_t>(3 * atom + 2)] = z;
            ++atom;
          }
        }
      }
    }
  }
  return ice;
}

int failures = 0;

void expect(bool ok, const std::string &what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL %s\n", what.c_str());
    ++failures;
  }
}

bool subset(const std::set<int> &inner, const std::set<int> &outer) {
  for (const int c : inner) {
    if (outer.count(c) == 0) {
      return false;
    }
  }
  return true;
}

} // namespace

int main(int argc, char **argv) {
  std::set<int> cpus;
  int width = 0;
  int setThreads = -1;
  for (int k = 1; k + 1 < argc; k += 2) {
    const std::string key = argv[k];
    if (key == "--cpus") {
      cpus = parseCpus(argv[k + 1]);
    } else if (key == "--width") {
      width = std::atoi(argv[k + 1]);
    } else if (key == "--set-threads") {
      setThreads = std::atoi(argv[k + 1]);
    }
  }
  if (cpus.empty() || width < 1) {
    std::fprintf(stderr, "usage: %s --cpus LIST --width N [--set-threads N]\n", argv[0]);
    return 2;
  }
  // The OpenMP runtime may already have bound this thread to its first
  // place, so the launch mask is the one the driver names.
  const std::set<int> start = maskOf(0);
  std::printf("launch cpus=%s start=%s width=%d set_threads=%d\n", showCpus(cpus).c_str(),
              showCpus(start).c_str(), width, setThreads);
  expect(subset(start, cpus), "start mask " + showCpus(start) + " is inside " + showCpus(cpus));

  if (setThreads >= 0) {
    expect(seams_set_num_threads(setThreads) == SEAMS_OK, "seams_set_num_threads");
  }
  const Frame ice = iceIh(16);
  seams_system *sys = seams_system_create();
  std::vector<int> types(static_cast<std::size_t>(ice.n), 1);
  std::vector<int> ids(static_cast<std::size_t>(ice.n));
  for (int i = 0; i < ice.n; ++i) {
    ids[static_cast<std::size_t>(i)] = i + 1;
  }
  expect(seams_system_set_frame(sys, ice.xyz.data(), ice.n, ice.box, types.data(), ids.data()) ==
             SEAMS_OK,
         "seams_system_set_frame");

  std::vector<int> labels(static_cast<std::size_t>(ice.n));
  std::vector<double> ql(static_cast<std::size_t>(ice.n));
  std::vector<double> qlbar(static_cast<std::size_t>(ice.n));
  int counts[9] = {};
  const std::vector<std::pair<const char *, std::function<int()>>> calls = {
      {"chill_plus_ovito_cutoff",
       [&] {
         return seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO,
                                 labels.data(), counts) == SEAMS_OK &&
                counts[SEAMS_LABEL_HEXAGONAL] == ice.n;
       }},
      {"chill_plus_published_mutual",
       [&] {
         return seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_PUBLISHED,
                                 labels.data(), counts) == SEAMS_OK &&
                counts[SEAMS_LABEL_HEXAGONAL] == ice.n;
       }},
      {"neighbours_cutoff",
       [&] {
         seams_neighbours *nb = nullptr;
         const bool ok = seams_neighbours_build(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK;
         seams_neighbours_destroy(nb);
         return ok;
       }},
      {"neighbours_knn",
       [&] {
         seams_neighbours *nb = nullptr;
         const bool ok = seams_neighbours_knn(sys, 4, 3.5, 1, 1, &nb) == SEAMS_OK;
         seams_neighbours_destroy(nb);
         return ok;
       }},
      {"q6_and_bar",
       [&] {
         return seams_steinhardt(sys, 3.5, 1, 6, ql.data(), ql.size(), qlbar.data(),
                                 qlbar.size()) == SEAMS_OK;
       }},
      {"rings",
       [&] {
         seams_rings *rings = nullptr;
         const bool ok = seams_rings_find(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK;
         seams_rings_destroy(rings);
         return ok;
       }},
  };
  for (int pass = 0; pass < 2; ++pass) {
    for (const auto &entry : calls) {
      const char *name = entry.first;
      const std::function<int()> &call = entry.second;
      bool ok = false;
      const Seen seen = watch([&] { ok = call() != 0; });
      const std::set<int> after = maskOf(0);
      std::printf("pass=%d call=%s ok=%d threads=%d masks=%s ran=%s after=%s\n", pass, name, ok,
                  seen.peakTasks, showCpus(seen.masks).c_str(), showCpus(seen.ran).c_str(),
                  showCpus(after).c_str());
      const std::string tag = std::string(name) + " pass " + std::to_string(pass);
      expect(ok, tag + " returned the ice labels");
      expect(after == cpus || (width == 1 && after == start),
             tag + ": process mask " + showCpus(after) + " after the call, launched on " +
                 showCpus(cpus));
      expect(subset(seen.masks, cpus),
             tag + ": a thread could run on " + showCpus(seen.masks) + ", launched on " +
                 showCpus(cpus));
      expect(subset(seen.ran, cpus),
             tag + ": a thread ran on " + showCpus(seen.ran) + ", launched on " + showCpus(cpus));
      if (width == 1) {
        expect(seen.peakTasks == 1,
               tag + ": " + std::to_string(seen.peakTasks) + " threads at width 1");
      }
    }
  }
  seams_system_destroy(sys);
  std::printf("cpus_effective=%zu\n", maskOf(0).size());
  if (failures != 0) {
    std::fprintf(stderr, "%d failures\n", failures);
    return 1;
  }
  std::printf("ok\n");
  return 0;
}
