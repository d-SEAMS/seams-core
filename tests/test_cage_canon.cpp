#include <catch2/catch_test_macros.hpp>

#include <cage_canon.hpp>

#include <vector>

TEST_CASE("nautyAvailable matches the compile-time flag", "[nauty]") {
#ifdef SEAMS_HAS_NAUTY
  REQUIRE(cage::nautyAvailable());
#else
  REQUIRE_FALSE(cage::nautyAvailable());
  REQUIRE(cage::canonicalCertificate({{0, 1, 2}}).empty());
  REQUIRE_FALSE(cage::isHexagonalPrism({{0, 1, 2, 3, 4, 5}}));
#endif
}

#ifdef SEAMS_HAS_NAUTY

static std::vector<std::vector<int>> hexPrismRings() {
  return {{0, 1, 2, 3, 4, 5},
          {6, 7, 8, 9, 10, 11},
          {0, 1, 7, 6},
          {1, 2, 8, 7},
          {2, 3, 9, 8},
          {3, 4, 10, 9},
          {4, 5, 11, 10},
          {5, 0, 6, 11}};
}

static std::vector<std::vector<int>> relabel(const std::vector<std::vector<int>> &rings,
                                             int n) {
  std::vector<std::vector<int>> out = rings;
  for (auto &ring : out) {
    for (int &a : ring) {
      a = n - 1 - a;
    }
  }
  return out;
}

TEST_CASE("hexagonal prism certificate is stable under relabel", "[nauty]") {
  const auto rings = hexPrismRings();
  REQUIRE(cage::isHexagonalPrism(rings));
  REQUIRE(cage::sameCertificate(rings, relabel(rings, 12)));
  REQUIRE_FALSE(cage::isHexagonalPrism({{0, 1, 2, 3, 4, 5}}));
}

#endif
