/*
** The pool-sized kernels match the one-thread result on ice Ic.
** Primitive rings, local and cutoff q4 and q6, CHILL+, the shell
** label, and one carried hop-ball update.
*/
#include <bop.hpp>
#include <franzblau.hpp>
#include <local_symmetry.hpp>
#include <neighbours.hpp>

#include "ice_ic.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <cstring>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

void setThreads(int n) {
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(n);
#else
  (void)n;
#endif
}

bool sameBits(const std::vector<double> &a, const std::vector<double> &b) {
  return a.size() == b.size() &&
         (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

std::vector<std::vector<int>> indexGraph(const Cloud &cloud) {
  auto listed = nneigh::neighListLinkcell(3.5, cloud, 1);
  if (listed.empty()) {
    listed = nneigh::neighListO(3.5, cloud, 1);
  }
  return nneigh::neighbourListByIndex(cloud, listed);
}

std::vector<int> iceLabels(const Cloud &cloud) {
  std::vector<int> labels(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    labels[static_cast<std::size_t>(i)] =
        static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
  }
  return labels;
}

std::vector<std::pair<int, std::uint64_t>> bondBits(const Cloud &cloud) {
  std::vector<std::pair<int, std::uint64_t>> out;
  out.reserve(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &bonds = cloud.pts[static_cast<std::size_t>(i)].c_ij;
    std::uint64_t bits = static_cast<std::uint64_t>(bonds.size());
    for (const auto &bond : bonds) {
      bits = bits * 1315423911u + static_cast<unsigned>(bond.classifier);
      std::uint64_t raw = 0;
      static_assert(sizeof(raw) == sizeof(bond.c_value));
      std::memcpy(&raw, &bond.c_value, sizeof(raw));
      bits ^= raw + 0x9e3779b97f4a7c15ull + (bits << 6) + (bits >> 2);
    }
    out.emplace_back(static_cast<int>(bonds.size()), bits);
  }
  return out;
}

} // namespace

TEST_CASE("pool-sized kernels match the one-thread result on ice Ic") {
  const Cloud ice = iceIc(4096);
  REQUIRE(ice.nop == 4096);
  const Cloud moved = perturbIce(ice);
  const auto graph = indexGraph(ice);
  const auto graphMoved = indexGraph(moved);
  REQUIRE(graph.size() == static_cast<std::size_t>(ice.nop));

  setThreads(1);
  const auto rings1 = primitive::ringNetwork(graph, 6);
  const auto q4Local1 = chill::steinhardtQl(ice, graph, 4);
  const auto q6Local1 = chill::steinhardtQl(ice, graph, 6);
  const auto q4Cut1 = chill::steinhardtQlCutoff(ice, 3.5, 4, 1, false, false);
  const auto q6Cut1 = chill::steinhardtQlCutoff(ice, 3.5, 6, 1, false, false);
  Cloud chill1 = ice;
  REQUIRE(chill::chillPlusPublishedCutoff(chill1, 3.5, 1));
  const auto labels1 = iceLabels(chill1);
  const auto bonds1 = bondBits(chill1);
  std::vector<int> sym1;
  symm::classifyIceCloud(ice, sym1);
  primitive::RingUpdater hop1(6);
  hop1.update(graph);
  const auto carried1 = hop1.update(graphMoved);

  setThreads(8);
  const auto rings8 = primitive::ringNetwork(graph, 6);
  const auto q4Local8 = chill::steinhardtQl(ice, graph, 4);
  const auto q6Local8 = chill::steinhardtQl(ice, graph, 6);
  const auto q4Cut8 = chill::steinhardtQlCutoff(ice, 3.5, 4, 1, false, false);
  const auto q6Cut8 = chill::steinhardtQlCutoff(ice, 3.5, 6, 1, false, false);
  Cloud chill8 = ice;
  REQUIRE(chill::chillPlusPublishedCutoff(chill8, 3.5, 1));
  const auto labels8 = iceLabels(chill8);
  const auto bonds8 = bondBits(chill8);
  std::vector<int> sym8;
  symm::classifyIceCloud(ice, sym8);
  primitive::RingUpdater hop8(6);
  hop8.update(graph);
  const auto carried8 = hop8.update(graphMoved);
  const auto rebuilt = primitive::ringNetwork(graphMoved, 6);

  REQUIRE(rings1 == rings8);
  REQUIRE(!rings1.empty());
  REQUIRE(sameBits(q4Local1.ql, q4Local8.ql));
  REQUIRE(sameBits(q6Local1.ql, q6Local8.ql));
  REQUIRE(sameBits(q4Cut1.ql, q4Cut8.ql));
  REQUIRE(sameBits(q6Cut1.ql, q6Cut8.ql));
  REQUIRE(labels1 == labels8);
  REQUIRE(bonds1 == bonds8);
  REQUIRE(sym1 == sym8);
  REQUIRE(carried1 == carried8);
  REQUIRE(carried1 == rebuilt);

  const auto full = [&](int threads) {
    setThreads(threads);
    auto listed = nneigh::neighListLinkcell(3.5, ice, 1);
    auto idx = nneigh::neighbourListByIndex(ice, listed);
    return primitive::ringNetwork(idx, 6);
  };
  const auto full1 = full(1);
  const auto full8 = full(8);
  REQUIRE(full1 == full8);
  REQUIRE(full1 == rings1);
}
