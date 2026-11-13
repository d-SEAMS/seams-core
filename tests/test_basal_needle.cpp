#include <catch2/catch_test_macros.hpp>

#include <seams_output.hpp>

#include <filesystem>
#include <vector>

TEST_CASE("unbonded basal rings keep a defined needle", "[output]") {
  std::vector<std::vector<int>> rings{{1, 2, 3, 4, 5, 6}, {7, 8, 9, 10, 11, 12}};
  std::vector<std::vector<int>> rows(12);
  for (int i = 0; i < 12; ++i) {
    rows[static_cast<std::size_t>(i)] = {i + 1};
  }
  const int wrote = sout::writeBasalRingsHex({0, 1}, 909001, rows, rings);
  REQUIRE(wrote == 1);
  std::filesystem::remove("../output/cages/hexBasalRings/basalRings909001.dat");
}
