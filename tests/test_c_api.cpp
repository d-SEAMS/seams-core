#include <catch2/catch_test_macros.hpp>
#include <seams_c_api.h>

#include <vector>

TEST_CASE("seams_chill_plus labels a four-atom tetrahedron without crashing",
          "[c_api]") {
  const double xyz[12] = {0, 0, 0, 2.7, 0, 0, 1.35, 2.34, 0, 1.35, 0.78, 2.21};
  const double box[6] = {0, 0, 0, 20, 20, 20};
  int labels[4] = {-1, -1, -1, -1};
  REQUIRE(seams_chill_plus(xyz, 4, box, labels) == 0);
  for (int i = 0; i < 4; i++) {
    REQUIRE(labels[i] >= 0);
    REQUIRE(labels[i] <= 8);
  }
}
