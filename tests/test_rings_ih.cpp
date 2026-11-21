// Copyright (c) d-SEAMS developers
// Census of seams_rings_find on ideal ice Ih.
//
// The symmetry counts follow from the lattice: each molecule sits in one
// puckered basal layer, whose hexagons number half its molecules, and the
// other 3N/2 of the 2N six-rings are prismatic faces of hexagonal cages.
// The fingerprints pin which rings carry which code, on the orthorhombic
// cell of scripts/gen_ideal_ice.py, cutoff 3.5 A, depth 6. Ring vertices
// are hashed after a rotation-invariant mix, then sorted, so the
// fingerprint does not depend on the order seams_rings_find returns the
// rings.

#include "ice_ih.hpp"
#include "seams.h"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace {

struct Census {
  int count = 0;
  int hist[8]{};
  int sym[12]{};
  std::uint64_t fingerprint = 0;
};

std::uint64_t fingerprintOf(seams_rings *rings) {
  const int n = seams_rings_count(rings);
  std::vector<int> sym(static_cast<std::size_t>(n));
  REQUIRE(seams_rings_labels(rings, nullptr, sym.data()) == SEAMS_OK);
  std::vector<std::uint64_t> keys(static_cast<std::size_t>(n));
  std::vector<int> atoms(16);
  for (int i = 0; i < n; ++i) {
    std::size_t nout = 0;
    REQUIRE(seams_rings_vertices(rings, i, atoms.data(), atoms.size(), &nout) == SEAMS_OK);
    std::uint64_t h = 1469598103934665603ULL;
    h ^= static_cast<std::uint64_t>(sym[static_cast<std::size_t>(i)]);
    h *= 1099511628211ULL;
    h ^= nout;
    h *= 1099511628211ULL;
    for (std::size_t k = 0; k < nout; ++k) {
      h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(atoms[k]));
      h *= 1099511628211ULL;
    }
    keys[static_cast<std::size_t>(i)] = h;
  }
  std::sort(keys.begin(), keys.end());
  std::uint64_t all = 1469598103934665603ULL;
  for (const auto key : keys) {
    all ^= key;
    all *= 1099511628211ULL;
  }
  return all;
}

Census run(int nx, int ny, int nz, int threads) {
  REQUIRE(seams_set_num_threads(threads) == SEAMS_OK);
  const IceIhCell cell = iceIhOrtho(nx, ny, nz);
  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  REQUIRE(seams_system_set_atoms(sys, cell.xyz.data(), static_cast<std::size_t>(cell.n),
                                 cell.box9, 9, nullptr, nullptr) == SEAMS_OK);
  seams_rings *rings = nullptr;
  REQUIRE(seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK);
  Census census;
  census.count = seams_rings_count(rings);
  REQUIRE(seams_rings_histogram(rings, census.hist, 8) == SEAMS_OK);
  std::vector<int> sym(static_cast<std::size_t>(census.count));
  REQUIRE(seams_rings_labels(rings, nullptr, sym.data()) == SEAMS_OK);
  for (int code : sym) {
    REQUIRE(code >= 0);
    REQUIRE(code < 12);
    ++census.sym[code];
  }
  census.fingerprint = fingerprintOf(rings);
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return census;
}

void expectOracle(const Census &census, int molecules, std::uint64_t fingerprint) {
  const int basal = molecules / 2;
  const int prismatic = 3 * molecules / 2;
  REQUIRE(census.count == 2 * molecules);
  REQUIRE(census.hist[6] == 2 * molecules);
  REQUIRE(census.hist[3] == 0);
  REQUIRE(census.hist[4] == 0);
  REQUIRE(census.hist[5] == 0);
  REQUIRE(census.sym[SEAMS_RING_HC_BASAL] == basal);
  REQUIRE(census.sym[SEAMS_RING_HC_PRISMATIC] == prismatic);
  REQUIRE(census.sym[SEAMS_RING_HEXAGON] == 0);
  REQUIRE(census.sym[SEAMS_RING_PRISM] == 0);
  REQUIRE(census.sym[SEAMS_RING_DEFORMED_PRISM] == 0);
  REQUIRE(census.sym[SEAMS_RING_MIXED_PRISM] == 0);
  REQUIRE(basal + prismatic == census.count);
  REQUIRE(census.fingerprint == fingerprint);
}

} // namespace

TEST_CASE("ice Ih ring census puts every six-ring on a hexagonal cage", "[rings_ih]") {
  const Census small = run(8, 8, 8, 1);
  expectOracle(small, 4096, 0x03e0b97b3f5999a7ULL);
  const Census mid = run(16, 16, 16, 1);
  expectOracle(mid, 32768, 0xe4a9479327642547ULL);
}

TEST_CASE("ice Ih ring census is the same at eight threads", "[rings_ih]") {
  const Census small = run(8, 8, 8, 8);
  expectOracle(small, 4096, 0x03e0b97b3f5999a7ULL);
  const Census mid = run(16, 16, 16, 8);
  expectOracle(mid, 32768, 0xe4a9479327642547ULL);
}
