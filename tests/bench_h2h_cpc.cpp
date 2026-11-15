// Copyright (c) d-SEAMS developers
/*
** One serial call for the head-to-head.
**
**   bench_h2h_cpc neighbours|rings TRAJ CUTOFF MAX_RING CSV
**
** The frame is read before the timer. neighbours times neighListO.
** rings builds that list first and times primitive::ringNetwork.
** The CSV is the identity output: neighbour counts, or the per-molecule
** ring-size histogram. Thread count is whatever OpenMP was given; the
** campaign sets one thread before the process starts, and this program
** asks for one thread again.
*/

#include <chill_plus_scheme.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#if defined(SEAMS_HAS_OPENMP)
#include <omp.h>
#endif

namespace {

using Clock = std::chrono::steady_clock;

void pinOneThread() {
#if defined(SEAMS_HAS_OPENMP)
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif
}

int threadCount() {
#if defined(SEAMS_HAS_OPENMP)
  return omp_get_max_threads();
#else
  return 1;
#endif
}

std::vector<int> neighbourIndices(const std::vector<int> &row, int self) {
  std::vector<int> out;
  out.reserve(row.size());
  for (std::size_t k = 0; k < row.size(); ++k) {
    const int j = row[k];
    if (j == self && (k == 0 || j == row.front())) {
      continue;
    }
    if (j == self) {
      continue;
    }
    out.push_back(j);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 6) {
    std::fprintf(stderr,
                 "usage: bench_h2h_cpc neighbours|rings|chill-plus TRAJ CUTOFF MAX_RING CSV\n");
    return 2;
  }
  pinOneThread();
  const std::string mode = argv[1];
  const std::string traj = argv[2];
  const double cutoff = std::atof(argv[3]);
  const int maxRing = std::atoi(argv[4]);
  const std::string csvPath = argv[5];
  if ((mode != "neighbours" && mode != "rings" && mode != "chill-plus") ||
      cutoff <= 0.0 || maxRing < 3) {
    std::fprintf(stderr, "bad arguments\n");
    return 2;
  }

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrjO(traj, 1, cloud, 1);
  if (cloud.nop <= 0) {
    std::fprintf(stderr, "empty frame\n");
    return 1;
  }
  const int n = cloud.nop;

  std::vector<std::vector<int>> idx;
  double ms = 0.0;
  std::map<int, int> tally;
  std::string chillHist;
  std::ofstream csv(csvPath);
  if (!csv) {
    std::perror(csvPath.c_str());
    return 1;
  }

  if (mode == "chill-plus") {
    const auto tNeigh0 = Clock::now();
    const auto nList = nneigh::neighListO(cutoff, cloud, 1);
    const auto tNeigh1 = Clock::now();
    const auto tClass0 = Clock::now();
    chill::applyChillPlus(cloud, nList, chill::ChillPlusRule::Ovito, false);
    std::map<std::string, int> labels;
    for (const auto &pt : cloud.pts) {
      const char *name = "water";
      switch (pt.iceType) {
      case molSys::atom_state_type::cubic:
        name = "cubic";
        break;
      case molSys::atom_state_type::hexagonal:
        name = "hexagonal";
        break;
      case molSys::atom_state_type::interfacial:
        name = "interfacial";
        break;
      case molSys::atom_state_type::clathrate:
        name = "hydrate";
        break;
      case molSys::atom_state_type::interClathrate:
        name = "interfacial-hydrate";
        break;
      case molSys::atom_state_type::water:
        name = "other";
        break;
      default:
        name = "other";
        break;
      }
      labels[name] += 1;
    }
    const auto tClass1 = Clock::now();
    const double neighMs =
        std::chrono::duration<double, std::milli>(tNeigh1 - tNeigh0).count();
    const double classMs =
        std::chrono::duration<double, std::milli>(tClass1 - tClass0).count();
    ms = neighMs + classMs;
    csv << "label,count\n";
    std::printf("neigh_ms %.6f\nclass_ms %.6f\n", neighMs, classMs);
    for (const auto &kv : labels) {
      csv << kv.first << ',' << kv.second << '\n';
      if (!chillHist.empty()) {
        chillHist.push_back(' ');
      }
      chillHist += kv.first;
      chillHist.push_back(':');
      chillHist += std::to_string(kv.second);
    }
  } else if (mode == "neighbours") {
    std::vector<std::vector<int>> nList;
    const auto t0 = Clock::now();
    nList = nneigh::neighListO(cutoff, cloud, 1);
    idx = nneigh::neighbourListByIndex(cloud, nList);
    const auto t1 = Clock::now();
    ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    csv << "index,count,neighbours\n";
    if (static_cast<int>(idx.size()) != n) {
      std::fprintf(stderr, "neighbour rows %zu nop %d\n", idx.size(), n);
      return 1;
    }
    for (int i = 0; i < n; ++i) {
      const std::vector<int> nbrs = neighbourIndices(idx[static_cast<std::size_t>(i)], i);
      tally[static_cast<int>(nbrs.size())] += 1;
      csv << i << ',' << nbrs.size() << ',';
      for (std::size_t k = 0; k < nbrs.size(); ++k) {
        if (k) {
          csv << ' ';
        }
        csv << nbrs[k];
      }
      csv << '\n';
    }
  } else {
    const auto nList = nneigh::neighListO(cutoff, cloud, 1);
    idx = nneigh::neighbourListByIndex(cloud, nList);
    std::vector<std::vector<int>> rings;
    const auto t0 = Clock::now();
    rings = primitive::ringNetwork(idx, maxRing);
    std::vector<std::vector<int>> hist(static_cast<std::size_t>(n),
                                       std::vector<int>(static_cast<std::size_t>(maxRing + 1), 0));
    for (const auto &ring : rings) {
      const int size = static_cast<int>(ring.size());
      if (size < 3 || size > maxRing) {
        continue;
      }
      tally[size] += 1;
      for (const int atom : ring) {
        if (atom >= 0 && atom < n) {
          hist[static_cast<std::size_t>(atom)][static_cast<std::size_t>(size)] += 1;
        }
      }
    }
    const auto t1 = Clock::now();
    ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    csv << "index";
    for (int size = 3; size <= maxRing; ++size) {
      csv << ',' << size;
    }
    csv << '\n';
    for (int i = 0; i < n; ++i) {
      csv << i;
      for (int size = 3; size <= maxRing; ++size) {
        csv << ',' << hist[static_cast<std::size_t>(i)][static_cast<std::size_t>(size)];
      }
      csv << '\n';
    }
  }
  if (!csv) {
    std::fprintf(stderr, "csv write failed\n");
    return 1;
  }

  std::printf("status ok\n");
  std::printf("molecules %d\n", n);
  std::printf("time_ms %.6f\n", ms);
  std::printf("omp_threads %d\n", threadCount());
  if (!chillHist.empty()) {
    std::printf("histogram %s\n", chillHist.c_str());
  } else {
    std::printf("tally");
    for (const auto &kv : tally) {
      std::printf(" %d:%d", kv.first, kv.second);
    }
    std::printf("\n");
  }
  return 0;
}
