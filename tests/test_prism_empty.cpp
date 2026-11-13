#include <catch2/catch_test_macros.hpp>

#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>

#include <vector>

TEST_CASE("empty basal rings are not a prism", "[prism]") {
  std::vector<std::vector<int>> rows{{0}};
  std::vector<int> basal1;
  std::vector<int> basal2;
  REQUIRE_FALSE(ring::basalPrismConditions(rows, basal1, basal2));
  REQUIRE_FALSE(prism3::basalPrismConditions(rows, basal1, basal2));
}
