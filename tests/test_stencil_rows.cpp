/*
** The interior-cell stencil and the packed ring rows.
** Orthorhombic and triclinic cells with atoms on every face and corner
** keep the same directed pairs as the branchy walk and as brute force.
** Packed rows and per-atom vectors hash to the same rings at 1, 2, 4
** and 8 threads. The cells are small.
*/
#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#include "ice_ic.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <utility>
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

void place(Cloud &cloud, double x, double y, double z) {
  molSys::Point<double> point{};
  point.x = x;
  point.y = y;
  point.z = z;
  point.type = 1;
  point.atomID = cloud.nop + 1;
  cloud.idIndexMap[point.atomID] = cloud.nop;
  cloud.pts.push_back(point);
  cloud.nop += 1;
}

Cloud facesAndCorners(bool tilt) {
  Cloud cloud;
  const double length = 12.0;
  if (tilt) {
    cloud.box = {length, length, length, 1.2, -0.35, 0.55};
  } else {
    cloud.box = {length, length, length};
  }
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double edge[2] = {0.15, length - 0.15};
  for (double x : edge) {
    for (double y : edge) {
      for (double z : edge) {
        place(cloud, x, y, z);
      }
    }
  }
  // One atom just inside each face, so the periodic image is in range.
  place(cloud, 0.2, 6.0, 6.0);
  place(cloud, length - 0.2, 6.0, 6.0);
  place(cloud, 6.0, 0.2, 6.0);
  place(cloud, 6.0, length - 0.2, 6.0);
  place(cloud, 6.0, 6.0, 0.2);
  place(cloud, 6.0, 6.0, length - 0.2);
  // Interior of the centre cell. These pairs never cross a periodic face.
  place(cloud, 6.0, 6.0, 6.0);
  place(cloud, 7.1, 6.2, 6.1);
  place(cloud, 6.2, 7.2, 6.3);
  place(cloud, 6.4, 6.3, 7.3);
  return cloud;
}

std::vector<std::pair<int, int>> pairsOf(const std::vector<int> &offsets,
                                        const std::vector<int> &cols) {
  std::vector<std::pair<int, int>> pairs;
  const int n = static_cast<int>(offsets.size()) - 1;
  for (int i = 0; i < n; ++i) {
    for (int k = offsets[static_cast<std::size_t>(i)];
         k < offsets[static_cast<std::size_t>(i) + 1]; ++k) {
      pairs.emplace_back(i, cols[static_cast<std::size_t>(k)]);
    }
  }
  return pairs;
}

std::vector<std::pair<int, int>> brutePairs(const Cloud &cloud, double cutoff) {
  const double rc2 = cutoff * cutoff;
  std::vector<std::pair<int, int>> pairs;
  for (int i = 0; i < cloud.nop; ++i) {
    for (int j = 0; j < cloud.nop; ++j) {
      if (i == j) {
        continue;
      }
      if (gen::periodicDistSq(cloud, i, j) < rc2) {
        pairs.emplace_back(i, j);
      }
    }
  }
  return pairs;
}

std::uint64_t mix(std::uint64_t h, std::uint64_t v) {
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
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

} // namespace

TEST_CASE("interior stencil matches the branchy walk and brute force", "[stencil]") {
  setThreads(1);
  for (const bool tilt : {false, true}) {
    const Cloud cloud = facesAndCorners(tilt);
    std::vector<int> fastOff;
    std::vector<int> fastCol;
    std::vector<int> branchOff;
    std::vector<int> branchCol;
    REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, fastOff, fastCol, false, nullptr, nullptr,
                                 nullptr, false));
    REQUIRE(nneigh::orthoFullCsr(cloud, nullptr, 3.5, branchOff, branchCol, false, nullptr,
                                 nullptr, nullptr, true));
    const auto fast = pairsOf(fastOff, fastCol);
    const auto branchy = pairsOf(branchOff, branchCol);
    const auto brute = brutePairs(cloud, 3.5);
    REQUIRE(fast == branchy);
    REQUIRE(fast == brute);
    REQUIRE(fast.size() > 8);
  }
}

TEST_CASE("packed rows and per-atom vectors hash the same rings", "[rings]") {
  const auto ice = iceIc(216);
  const auto listed = nneigh::neighListO(3.5, ice, 1);
  const auto idx = nneigh::neighbourListByIndex(ice, listed);
  REQUIRE(!idx.empty());
  std::uint64_t first = 0;
  for (const int threads : {1, 2, 4, 8}) {
    setThreads(threads);
    const auto packed =
        hashRings(primitive::ringNetwork(idx, 6, primitive::RingRows::Packed));
    const auto vectors =
        hashRings(primitive::ringNetwork(idx, 6, primitive::RingRows::Vectors));
    REQUIRE(packed == vectors);
    REQUIRE(packed != 0);
    if (threads == 1) {
      first = packed;
    } else {
      REQUIRE(packed == first);
    }
  }
}
