#include <catch2/catch_test_macros.hpp>

#include <gpu_residency.hpp>

TEST_CASE("footprint scales with frames and atoms") {
  const auto one = gpu::estimateFootprint(4096, 1);
  const auto eleven = gpu::estimateFootprint(4096, 11);
  REQUIRE(one.totalBytes > 0);
  REQUIRE(eleven.totalBytes > one.totalBytes);
  REQUIRE(eleven.xyzBytes == 11 * one.xyzBytes);
}

TEST_CASE("zero atoms cannot reside") {
  gpu::DeviceInfo fake;
  fake.available = true;
  fake.freeBytes = 80ull * 1024ull * 1024ull * 1024ull;
  REQUIRE(gpu::maxResidentFrames(fake, 0) == 0);
}

TEST_CASE("an A100-sized budget holds the cubic trajectory") {
  gpu::DeviceInfo fake;
  fake.available = true;
  fake.freeBytes = 40ull * 1024ull * 1024ull * 1024ull;
  const int maxF = gpu::maxResidentFrames(fake, 4096);
  REQUIRE(maxF >= 11);
}

TEST_CASE("affiliation working set is counted in the footprint") {
  const auto foot = gpu::estimateFootprint(4096, 11, 16, 16);
  REQUIRE(foot.ringsBytes > foot.xyzBytes);
  REQUIRE(foot.labelBytes >= 11ull * 4096ull * 2ull * sizeof(int));
}

TEST_CASE("estimateFootprint budgets basal-pair device buffers") {
  const auto foot = gpu::estimateFootprint(4096, 11, 16, 16);
  const std::size_t maxRings = 4096ull * 16ull;
  const std::size_t maxPairs = maxRings * 8ull;
  const std::size_t pairBytes = 11ull * maxPairs * 2ull * sizeof(int);
  REQUIRE(foot.totalBytes >= pairBytes);
}

TEST_CASE("mutual four-nearest back-check uses the first four only") {
  // Host replica of the device contract: i keeps its first 4, and the
  // edge is mutual only if i is among j's first 4, not anywhere in kMax.
  constexpr int kMax = 16;
  int cols[2 * kMax];
  for (int t = 0; t < kMax; t++) {
    cols[t] = -1;
    cols[kMax + t] = -1;
  }
  cols[0] = 1;
  cols[1] = 2;
  cols[2] = 3;
  cols[3] = 4;
  cols[kMax + 0] = 5;
  cols[kMax + 1] = 6;
  cols[kMax + 2] = 7;
  cols[kMax + 3] = 8;
  cols[kMax + 4] = 0;
  auto amongFirst4 = [](const int *row, int who) {
    for (int t = 0; t < 4; t++) {
      if (row[t] == who) {
        return true;
      }
    }
    return false;
  };
  REQUIRE(amongFirst4(cols, 1));
  REQUIRE_FALSE(amongFirst4(cols + kMax, 0));
}
