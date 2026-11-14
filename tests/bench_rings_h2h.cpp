/*
** Ring search against a cutoff graph.
**
**   bench_rings_h2h [--repl] TRAJ frame atomType cutoff maxRing [reps warmup]
**
** One-shot prints one sample_ms line per timed call. Those samples are
** the ring search on a graph built before the timer.
**
** --repl builds the graph, prints a ready line, then answers each line:
**   go    ring search on the graph built before the timer
**   in    neighbour list, index, and ring search inside one timer
**   quit  stop
*/

#include <franzblau.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;

std::string cpuModel() {
  FILE *fh = std::fopen("/proc/cpuinfo", "r");
  if (fh == nullptr) {
    return "unknown";
  }
  char line[512];
  std::string model = "unknown";
  while (std::fgets(line, sizeof line, fh) != nullptr) {
    std::string text(line);
    const auto mark = text.find(':');
    if (text.rfind("model name", 0) == 0 && mark != std::string::npos) {
      model = text.substr(mark + 1);
      while (!model.empty() && (model.back() == '\n' || model.back() == ' ')) {
        model.pop_back();
      }
      while (!model.empty() && model.front() == ' ') {
        model.erase(model.begin());
      }
      break;
    }
  }
  std::fclose(fh);
  std::string clean;
  for (std::size_t i = 0; i < model.size();) {
    if (model.compare(i, 3, "(R)") == 0 || model.compare(i, 3, "(r)") == 0) {
      i += 3;
      continue;
    }
    clean.push_back(model[i]);
    i++;
  }
  std::string collapsed;
  bool space = false;
  for (char c : clean) {
    if (c == ' ') {
      space = true;
      continue;
    }
    if (space && !collapsed.empty()) {
      collapsed.push_back(' ');
    }
    space = false;
    collapsed.push_back(c);
  }
  return collapsed.empty() ? "unknown" : collapsed;
}

std::string histogram(const std::vector<std::vector<int>> &rings) {
  std::map<int, int> counts;
  for (const auto &ring : rings) {
    counts[static_cast<int>(ring.size())] += 1;
  }
  std::string text;
  for (const auto &item : counts) {
    if (!text.empty()) {
      text.push_back(' ');
    }
    text += std::to_string(item.first) + ":" + std::to_string(item.second);
  }
  if (text.empty()) {
    text = "none";
  }
  return text;
}

double millis(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}
} // namespace

int main(int argc, char **argv) {
  bool repl = false;
  int arg = 1;
  if (arg < argc && std::string(argv[arg]) == "--repl") {
    repl = true;
    arg++;
  }
  if (argc - arg < 5) {
    std::fprintf(stderr,
                 "usage: bench_rings_h2h [--repl] TRAJ frame atomType cutoff "
                 "maxRing [reps warmup]\n");
    return 2;
  }
  const std::string traj = argv[arg];
  const int frame = std::atoi(argv[arg + 1]);
  const int typeI = std::atoi(argv[arg + 2]);
  const double cutoff = std::atof(argv[arg + 3]);
  const int maxRing = std::atoi(argv[arg + 4]);
  const int reps = argc - arg > 5 ? std::atoi(argv[arg + 5]) : 8;
  const int warmup = argc - arg > 6 ? std::atoi(argv[arg + 6]) : 2;

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, frame, cloud, typeI);
  if (cloud.nop <= 0) {
    std::fprintf(stderr, "empty frame %d of %s\n", frame, traj.c_str());
    return 1;
  }
  const auto nList = nneigh::neighListO(cutoff, cloud, typeI);
  const auto idx = nneigh::neighbourListByIndex(cloud, nList);
  auto rings = primitive::ringNetwork(idx, maxRing);
  const std::string hist = histogram(rings);
  const auto expect = rings.size();

  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  if (repl) {
    std::printf("ready nop %d hist %s\n", cloud.nop, hist.c_str());
    char line[64];
    while (std::fgets(line, sizeof line, stdin) != nullptr) {
      if (line[0] == 'q') {
        break;
      }
      const bool inside = line[0] == 'i';
      if (line[0] != 'g' && !inside) {
        continue;
      }
      const auto t0 = Clock::now();
      if (inside) {
        // The published inside row is this whole call. Building the
        // list before t0 made that row a second copy of the prebuilt search.
        const auto built = nneigh::neighListO(cutoff, cloud, typeI);
        const auto builtIdx = nneigh::neighbourListByIndex(cloud, built);
        rings = primitive::ringNetwork(builtIdx, maxRing);
      } else {
        rings = primitive::ringNetwork(idx, maxRing);
      }
      const auto t1 = Clock::now();
      if (rings.size() != expect) {
        std::fprintf(stderr, "ring count changed\n");
        return 1;
      }
      std::printf("ms %.6f\n", millis(t0, t1));
    }
    return 0;
  }

  std::printf("nproc %ld\n", sysconf(_SC_NPROCESSORS_ONLN));
  std::printf("cpu_model %s\n", cpuModel().c_str());
  std::printf("nop %d\n", cloud.nop);
  std::printf("hist %s\n", hist.c_str());
  for (int i = 0; i < warmup + reps; i++) {
    const auto t0 = Clock::now();
    rings = primitive::ringNetwork(idx, maxRing);
    const auto t1 = Clock::now();
    if (rings.size() != expect) {
      std::fprintf(stderr, "ring count changed\n");
      return 1;
    }
    if (i >= warmup) {
      std::printf("sample_ms %.6f\n", millis(t0, t1));
    }
  }
  return 0;
}
