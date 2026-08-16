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
