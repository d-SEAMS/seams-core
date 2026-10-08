#include <catch2/catch_test_macros.hpp>

#include "domain_frames.hpp"

#include <domain.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <vector>

using Cloud = seams::domain::Cloud;

namespace {

// Shortest image over the wrapped difference and its 26 neighbours
double minImageSq(const gen::FracBox &b, const Cloud &c, int i, int j) {
  const auto &p = c.pts[static_cast<std::size_t>(i)];
  const auto &q = c.pts[static_cast<std::size_t>(j)];
  const auto d = gen::fracDelta(b, p.x, p.y, p.z, q.x, q.y, q.z);
  double best = 1e300;
  for (int na = -1; na <= 1; na++) {
    for (int nb = -1; nb <= 1; nb++) {
      for (int nc = -1; nc <= 1; nc++) {
        const double x = d[0] + na * b.lx + nb * b.xy + nc * b.xz;
        const double y = d[1] + nb * b.ly + nc * b.yz;
        const double z = d[2] + nc * b.lz;
        best = std::min(best, x * x + y * y + z * z);
      }
    }
  }
  return best;
}

} // namespace

TEST_CASE("Hilbert keys walk every cell by face steps", "[domain]") {
  for (int bits = 1; bits <= 7; bits++) {
    INFO("bits " << bits);
    const std::uint32_t side = 1u << bits;
    std::vector<std::array<int, 3>> at(side * side * side, {-1, -1, -1});
    bool permutation = true;
    for (std::uint32_t z = 0; z < side; z++) {
      for (std::uint32_t y = 0; y < side; y++) {
        for (std::uint32_t x = 0; x < side; x++) {
          const auto key = seams::domain::hilbertKey(x, y, z, bits);
          if (key >= at.size() || at[key][0] != -1) {
            permutation = false;
            continue;
          }
          at[key] = {static_cast<int>(x), static_cast<int>(y),
                     static_cast<int>(z)};
        }
      }
    }
    REQUIRE(permutation);
    bool faceSteps = true;
    for (std::size_t k = 1; k < at.size(); k++) {
      faceSteps = faceSteps && std::abs(at[k][0] - at[k - 1][0]) +
                                       std::abs(at[k][1] - at[k - 1][1]) +
                                       std::abs(at[k][2] - at[k - 1][2]) ==
                                   1;
    }
    REQUIRE(faceSteps);
  }
}

TEST_CASE("Shares partition the frame and hold every atom near an owned one",
          "[domain]") {
  const double halo = 4.0;
  for (const auto &frame : {tiltedFrame(800, 32.0, 0.0, 0.0, 0.0, 1),
                            tiltedFrame(800, 32.0, 6.0, -4.0, 3.0, 2)}) {
    const gen::FracBox b = gen::makeFracBox(frame);
    for (const int nRanks : {1, 2, 3, 5, 8, 64, 128}) {
      std::vector<int> owners(frame.pts.size(), 0);
      for (int rank = 0; rank < nRanks; rank++) {
        const auto share = seams::domain::decompose(frame, rank, nRanks, halo);
        REQUIRE(share.local.size() == share.owned.size());
        REQUIRE(std::is_sorted(share.local.begin(), share.local.end()));
        std::vector<char> local(frame.pts.size(), 0);
        for (std::size_t k = 0; k < share.local.size(); k++) {
          local[static_cast<std::size_t>(share.local[k])] = 1;
          owners[static_cast<std::size_t>(share.local[k])] += share.owned[k];
        }
        for (std::size_t k = 0; k < share.local.size(); k++) {
          if (!share.owned[k]) {
            continue;
          }
          for (int j = 0; j < frame.nop; j++) {
            if (minImageSq(b, frame, share.local[k], j) < halo * halo) {
              REQUIRE(local[static_cast<std::size_t>(j)]);
            }
          }
        }
        if (nRanks > 1) {
          REQUIRE(share.local.size() < frame.pts.size());
        }
      }
      for (const int o : owners) {
        REQUIRE(o == 1);
      }
    }
  }
}

TEST_CASE("An empty frame splits into empty shares", "[domain]") {
  const auto empty = tiltedFrame(0, 20.0, 0.0, 0.0, 0.0, 6);
  for (int rank = 0; rank < 4; rank++) {
    REQUIRE(seams::domain::decompose(empty, rank, 4, 7.0).local.empty());
    REQUIRE(seams::domain::rings(empty, 3.5, 6, rank, 4).empty());
  }
}

TEST_CASE("A source mask keeps exactly the rings those sources enumerate",
          "[domain]") {
  const auto frame = tiltedFrame(3000, 44.9, 0.0, 0.0, 0.0, 3);
  const auto nList = nneigh::getNewNeighbourListByIndex(frame, 3.6);
  const auto all = primitive::ringNetwork(nList, 6);
  std::vector<char> mask(nList.size(), 0);
  for (std::size_t v = 0; v < mask.size(); v += 3) {
    mask[v] = 1;
  }
  std::vector<std::vector<int>> expected;
  for (const auto &ring : all) {
    if (mask[static_cast<std::size_t>(ring[0])]) {
      expected.push_back(ring);
    }
  }
  REQUIRE(!expected.empty());
  REQUIRE(primitive::ringNetwork(nList, 6, mask) == expected);
}

TEST_CASE("Ranks' rings together are the frame's rings", "[domain]") {
  const double cutoff = 3.4;
  for (const auto &frame : {tiltedFrame(16000, 78.4, 0.0, 0.0, 0.0, 4),
                            tiltedFrame(16000, 78.4, 9.0, 5.0, -7.0, 5)}) {
    for (const int depth : {3, 4, 5, 6, 7, 8}) {
      const auto whole = sortedRowRings(frame, cutoff, depth);
      REQUIRE(!whole.empty());
      for (const int nRanks : {1, 2, 3, 5, 8}) {
        std::vector<std::vector<int>> joined;
        for (int rank = 0; rank < nRanks; rank++) {
          const auto mine =
              seams::domain::rings(frame, cutoff, depth, rank, nRanks);
          joined.insert(joined.end(), mine.begin(), mine.end());
        }
        std::stable_sort(joined.begin(), joined.end(),
                         [](const auto &a, const auto &b) { return a[0] < b[0]; });
        REQUIRE(joined == whole);
      }
    }
  }
}
