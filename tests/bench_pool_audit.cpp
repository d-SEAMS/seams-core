/*
** Time one kernel on ideal ice Ic. The crystal is built in memory.
** rings, q4-local, q6-local, q4-cutoff, q6-cutoff, chill, symmetry, hop.
*/
#include <bop.hpp>
#include <franzblau.hpp>
#include <local_symmetry.hpp>
#include <neighbours.hpp>

#include "ice_ic.hpp"

#include <chrono>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

double msSince(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - t0)
      .count();
}

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
  h ^= v + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
  return h;
}

std::uint64_t hashRings(const std::vector<std::vector<int>> &rings) {
  std::uint64_t h = rings.size();
  for (const auto &ring : rings) {
    h = mix(h, ring.size());
    for (int v : ring) {
      h = mix(h, static_cast<std::uint64_t>(v));
    }
  }
  return h;
}

std::uint64_t hashDoubles(const std::vector<double> &values) {
  std::uint64_t h = values.size();
  for (double value : values) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    h = mix(h, bits);
  }
  return h;
}

std::uint64_t hashInts(const std::vector<int> &values) {
  std::uint64_t h = values.size();
  for (int value : values) {
    h = mix(h, static_cast<std::uint64_t>(value));
  }
  return h;
}

std::vector<std::vector<int>> indexGraph(const Cloud &cloud) {
  auto listed = nneigh::neighListLinkcell(3.5, cloud, 1);
  if (listed.empty()) {
    listed = nneigh::neighListO(3.5, cloud, 1);
  }
  return nneigh::neighbourListByIndex(cloud, listed);
}

struct Sample {
  double ms = 0;
  std::uint64_t fingerprint = 0;
};

Sample once(const std::string &kernel, const Cloud &ice, const Cloud &moved,
            const std::vector<std::vector<int>> &graph,
            const std::vector<std::vector<int>> &graphMoved) {
  const auto t0 = std::chrono::steady_clock::now();
  std::uint64_t fingerprint = 0;
  if (kernel == "rings") {
    const auto &rings = primitive::ringNetwork(graph, 6);
    const double ms = msSince(t0);
    fingerprint = hashRings(rings);
    return {ms, fingerprint};
  }
  if (kernel == "q4-local") {
    const auto q = chill::steinhardtQl(ice, graph, 4);
    const double ms = msSince(t0);
    return {ms, hashDoubles(q.ql)};
  }
  if (kernel == "q6-local") {
    const auto q = chill::steinhardtQl(ice, graph, 6);
    const double ms = msSince(t0);
    return {ms, hashDoubles(q.ql)};
  }
  if (kernel == "q4-cutoff") {
    const auto q = chill::steinhardtQlCutoff(ice, 3.5, 4, 1, false, false);
    const double ms = msSince(t0);
    return {ms, hashDoubles(q.ql)};
  }
  if (kernel == "q6-cutoff") {
    const auto q = chill::steinhardtQlCutoff(ice, 3.5, 6, 1, false, false);
    const double ms = msSince(t0);
    return {ms, hashDoubles(q.ql)};
  }
  if (kernel == "chill") {
    Cloud work = ice;
    const auto tChill = std::chrono::steady_clock::now();
    const bool ok = chill::chillPlusPublishedCutoff(work, 3.5, 1);
    const double ms = msSince(tChill);
    std::uint64_t h = ok ? 1 : 0;
    for (int i = 0; i < work.nop; ++i) {
      h = mix(h, static_cast<std::uint64_t>(work.pts[static_cast<std::size_t>(i)].iceType));
    }
    return {ms, h};
  }
  if (kernel == "symmetry") {
    std::vector<int> kind;
    symm::classifyIceCloud(ice, kind);
    const double ms = msSince(t0);
    return {ms, hashInts(kind)};
  }
  if (kernel == "hop") {
    primitive::RingUpdater hop(6);
    hop.update(graph);
    const auto t1 = std::chrono::steady_clock::now();
    const auto &rings = hop.update(graphMoved);
    const double ms = msSince(t1);
    return {ms, hashRings(rings)};
  }
  return {};
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "usage: bench_pool_audit KERNEL MOLECULES WARMUP REPEATS\n";
    return 2;
  }
  const std::string kernel = argv[1];
  const int molecules = std::atoi(argv[2]);
  const int warmup = std::atoi(argv[3]);
  const int repeats = std::atoi(argv[4]);
  if (warmup < 0 || repeats < 1) {
    std::cerr << "bad warmup or repeats\n";
    return 2;
  }
  const Cloud ice = iceIc(molecules);
  if (ice.nop != molecules) {
    std::cerr << "unsupported molecule count " << molecules << "\n";
    return 2;
  }
  const Cloud moved = perturbIce(ice);
  const auto graph = indexGraph(ice);
  const auto graphMoved = indexGraph(moved);
  if (graph.empty()) {
    std::cerr << "empty neighbour list\n";
    return 1;
  }
  const Sample first = once(kernel, ice, moved, graph, graphMoved);
  if (first.fingerprint == 0 && kernel != "chill") {
    std::cerr << "empty " << kernel << "\n";
    return 1;
  }
  std::cout << "kernel " << kernel << "\n";
  std::cout << "molecules " << molecules << "\n";
  std::cout << "fingerprint " << first.fingerprint << "\n";
  std::cout.flush();
  for (int rep = 0; rep < warmup + repeats; ++rep) {
    const Sample got = once(kernel, ice, moved, graph, graphMoved);
    if (got.fingerprint != first.fingerprint) {
      std::cerr << "fingerprint changed during timing\n";
      return 1;
    }
    if (rep >= warmup) {
      std::cout << "ms " << got.ms << "\n";
      std::cout.flush();
    }
  }
  return 0;
}
