#include <catch2/catch_test_macros.hpp>

#include "domain_frames.hpp"

#include <domain.hpp>
#include <ring.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <neighbours.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <thread>
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

TEST_CASE("The partition halo is floor(k/2) bonds", "[domain]") {
  // formal/Formal/Halo.lean, primitive_ring_subset_ball.
  REQUIRE(seams::domain::kRingHaloHalf == 2);
  for (int k = 0; k <= 12; k++) {
    REQUIRE(seams::domain::ringHaloBonds(k) == k / 2);
    REQUIRE(seams::domain::ringHops(k) == seams::domain::ringHaloBonds(k));
    REQUIRE(seams::domain::ringRowHops(k) == seams::domain::ringHaloBonds(k));
  }
  // gatherWork takes this many cutoffs. 6 -> 3 bonds, 10 -> 5 bonds.
  REQUIRE(seams::domain::ringHaloBonds(6) == 3);
  REQUIRE(seams::domain::ringHaloBonds(7) == 3);
  REQUIRE(seams::domain::ringHaloBonds(8) == 4);
  REQUIRE(seams::domain::ringHaloBonds(10) == 5);
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

TEST_CASE("A summary of the ranks' rings is the frame's summary", "[domain]") {
  const double cutoff = 3.4;
  const auto frame = tiltedFrame(16000, 78.4, 9.0, 5.0, -7.0, 5);
  for (const int depth : {6, 8}) {
    const auto whole = sortedRowRings(frame, cutoff, depth);
    const auto want = seams::domain::summarize(whole, frame.nop, depth);
    // ring::assignPolygonType with every atom starting unclassified
    std::vector<int> types(static_cast<std::size_t>(frame.nop), 1);
    std::vector<int> nRingList;
    ring::assignPolygonType(whole, types, nRingList);
    REQUIRE(want.atomType == types);
    for (const int nRanks : {2, 3, 8}) {
      std::vector<std::vector<int>> joined;
      for (int rank = 0; rank < nRanks; rank++) {
        const auto mine = seams::domain::rings(frame, cutoff, depth, rank, nRanks);
        joined.insert(joined.end(), mine.begin(), mine.end());
      }
      const auto got = seams::domain::summarize(joined, frame.nop, depth);
      REQUIRE(got.sizeCount == want.sizeCount);
      REQUIRE(got.atomType == want.atomType);
    }
  }
}

TEST_CASE("a halo source is the atom share.local names", "[domain]") {
  Cloud cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double z[4] = {1.0, 15.0, 11.0, 15.0};
  for (int i = 0; i < 4; i++) {
    molSys::Point<double> p;
    p.x = 10.0;
    p.y = 10.0;
    p.z = z[i];
    p.atomID = i + 1;
    p.molID = p.atomID;
    p.type = 1;
    cloud.pts.push_back(p);
  }
  cloud.nop = static_cast<int>(cloud.pts.size());
  // local[k] is a frame index, and it is not k. Rank 0's block on this
  // box is the low half of z. Atom 2 sits just outside that block. Atom
  // 1 and atom 3 sit at the far side of the other block.
  seams::domain::Share share;
  share.local = {0, 2, 3};
  share.owned = {1, 0, 0};
  share.owner = {0, 1, 1};
  REQUIRE(share.local[1] != 1);
  REQUIRE(share.local[2] != 2);
  std::vector<char> source;
  const double inner = 3.0;
  REQUIRE(inner >= 0.0);
  seams::domain::markSources(cloud, share, inner, source);
  REQUIRE(source.size() == 3);
  REQUIRE(source[0] == 1);
  REQUIRE(source[1] == 1);
  REQUIRE(source[2] == 0);
}

namespace {

Cloud lattice(int n, double a) {
  Cloud cloud;
  const double L = n * a;
  cloud.box = {L, L, L};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  for (int ix = 0; ix < n; ix++) {
    for (int iy = 0; iy < n; iy++) {
      for (int iz = 0; iz < n; iz++) {
        molSys::Point<double> p;
        p.x = (ix + 0.5) * a;
        p.y = (iy + 0.5) * a;
        p.z = (iz + 0.5) * a;
        p.type = 1;
        p.atomID = cloud.nop + 1;
        p.molID = p.atomID;
        cloud.idIndexMap[p.atomID] = cloud.nop;
        cloud.pts.push_back(p);
        cloud.nop += 1;
      }
    }
  }
  return cloud;
}

} // namespace

TEST_CASE("overlapping shell and ring calls keep their own scratch", "[domain][threads]") {
  const auto cloudA = lattice(6, 3.2);
  const auto cloudB = lattice(6, 4.8);
  const double cutoff = 3.5;
  const std::vector<std::vector<int>> shellA =
      nneigh::reusedShellList(cutoff, cloudA, nullptr);
  const std::vector<std::vector<int>> shellB =
      nneigh::reusedShellList(cutoff, cloudB, nullptr);
  REQUIRE(shellA != shellB);
  const auto ringsA = seams::domain::rings(cloudA, cutoff, 4, 0, 1);
  const auto ringsB = seams::domain::rings(cloudB, cutoff, 4, 0, 1);
  REQUIRE(ringsA != ringsB);
  std::atomic<int> ready{0};
  std::atomic<int> bad{0};
  auto run = [&](const Cloud &cloud, const std::vector<std::vector<int>> &shell,
                 const std::vector<std::vector<int>> &rings) {
    ready.fetch_add(1);
    while (ready.load() < 2) {
    }
    for (int n = 0; n < 24; n++) {
      const auto &got = nneigh::reusedShellList(cutoff, cloud, nullptr);
      if (got != shell) {
        bad.fetch_add(1);
      }
      if (seams::domain::rings(cloud, cutoff, 4, 0, 1) != rings) {
        bad.fetch_add(1);
      }
    }
  };
  std::thread first(run, std::cref(cloudA), std::cref(shellA), std::cref(ringsA));
  std::thread second(run, std::cref(cloudB), std::cref(shellB), std::cref(ringsB));
  first.join();
  second.join();
  REQUIRE(bad.load() == 0);
}
