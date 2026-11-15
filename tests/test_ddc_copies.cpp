#include <catch2/catch_test_macros.hpp>

#include <topo_bulk.hpp>

#include <cstdint>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

extern "C" std::uint64_t seams_ddc_row_copies(void);

namespace {

std::vector<ring::EquatorChoice> cliques(int groups, int width) {
  std::vector<ring::EquatorChoice> rows(static_cast<std::size_t>(groups * width));
  for (int group = 0; group < groups; ++group) {
    for (int atom = 0; atom < width; ++atom) {
      auto &row = rows[static_cast<std::size_t>(group * width + atom)];
      row.n = width - 1;
      int slot = 0;
      for (int other = 0; other < width; ++other) {
        if (other == atom) {
          continue;
        }
        row.put(slot++, group * width + other);
      }
    }
  }
  return rows;
}

} // namespace

TEST_CASE("parallel equator rows count the same copies as one thread", "[ddc]") {
  const auto rows = cliques(40, 10);
  std::vector<int> order(rows.size());
  for (int i = 0; i < static_cast<int>(order.size()); ++i) {
    order[static_cast<std::size_t>(i)] = i;
  }
#ifdef SEAMS_HAS_OPENMP
  const int saved = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
  const auto serial = ring::maximumEquators(rows, order, nullptr);
  const auto copies1 = seams_ddc_row_copies();
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(8);
#endif
  REQUIRE(copies1 > 0);
  for (int trial = 0; trial < 40; ++trial) {
    const auto parallel = ring::maximumEquators(rows, order, nullptr);
    const auto copies8 = seams_ddc_row_copies();
    REQUIRE(parallel == serial);
    REQUIRE(copies8 == copies1);
  }
#ifdef SEAMS_HAS_OPENMP
  omp_set_num_threads(saved);
#endif
}
