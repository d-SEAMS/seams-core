#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

// Helper to build a simple 2-atom cloud for distance tests
static molSys::PointCloud<molSys::Point<double>, double>
makeTwoAtomCloud(double x0, double y0, double z0, double x1, double y1,
                 double z1, double boxLen = 10.0) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {boxLen, boxLen, boxLen};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  molSys::Point<double> p0;
  p0.type = 1;
  p0.atomID = 0;
  p0.molID = 0;
  p0.x = x0;
  p0.y = y0;
  p0.z = z0;

  molSys::Point<double> p1;
  p1.type = 1;
  p1.atomID = 1;
  p1.molID = 1;
  p1.x = x1;
  p1.y = y1;
  p1.z = z1;

  cloud.pts.push_back(p0);
  cloud.pts.push_back(p1);
  cloud.nop = 2;
  cloud.idIndexMap[0] = 0;
  cloud.idIndexMap[1] = 1;
  return cloud;
}

// -- periodicDist tests --

TEST_CASE("periodicDist for atoms within the same image", "[generic]") {
  auto cloud = makeTwoAtomCloud(1.0, 0.0, 0.0, 4.0, 0.0, 0.0);
  double d = gen::periodicDist(cloud, 0, 1);
  REQUIRE_THAT(d, Catch::Matchers::WithinAbs(3.0, 1e-10));
}

TEST_CASE("periodicDist across periodic boundary", "[generic]") {
  auto cloud = makeTwoAtomCloud(0.5, 0.0, 0.0, 9.5, 0.0, 0.0);
  double d = gen::periodicDist(cloud, 0, 1);
  // Periodic distance should be 1.0
  REQUIRE_THAT(d, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("periodicDist for identical atoms is zero", "[generic]") {
  auto cloud = makeTwoAtomCloud(3.0, 4.0, 5.0, 3.0, 4.0, 5.0);
  double d = gen::periodicDist(cloud, 0, 1);
  REQUIRE_THAT(d, Catch::Matchers::WithinAbs(0.0, 1e-10));
}

// -- distance tests (no PBC) --

TEST_CASE("distance without PBC", "[generic]") {
  auto cloud = makeTwoAtomCloud(0.0, 0.0, 0.0, 3.0, 4.0, 0.0);
  double d = gen::distance(cloud, 0, 1);
  REQUIRE_THAT(d, Catch::Matchers::WithinAbs(5.0, 1e-10));
}

TEST_CASE("distance across box edge does NOT wrap", "[generic]") {
  auto cloud = makeTwoAtomCloud(0.5, 0.0, 0.0, 9.5, 0.0, 0.0);
  double d = gen::distance(cloud, 0, 1);
  // Without PBC, distance is 9.0
  REQUIRE_THAT(d, Catch::Matchers::WithinAbs(9.0, 1e-10));
}

// -- relDist tests --

TEST_CASE("relDist returns signed relative distances with PBC", "[generic]") {
  auto cloud = makeTwoAtomCloud(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  auto dr = gen::relDist(cloud, 0, 1);

  REQUIRE_THAT(dr[0], Catch::Matchers::WithinAbs(-3.0, 1e-10));
  REQUIRE_THAT(dr[1], Catch::Matchers::WithinAbs(-3.0, 1e-10));
  REQUIRE_THAT(dr[2], Catch::Matchers::WithinAbs(-3.0, 1e-10));
}

TEST_CASE("relDist wraps across boundary", "[generic]") {
  auto cloud = makeTwoAtomCloud(0.5, 0.0, 0.0, 9.5, 0.0, 0.0);
  auto dr = gen::relDist(cloud, 0, 1);

  // 0.5 - 9.5 = -9.0, which is < -5.0, so wrap: -9.0 + 10.0 = 1.0
  REQUIRE_THAT(dr[0], Catch::Matchers::WithinAbs(1.0, 1e-10));
}

// -- tokenizer tests --

TEST_CASE("tokenizer splits whitespace-delimited string", "[generic]") {
  auto tokens = gen::tokenizer("hello world foo");
  REQUIRE(tokens.size() == 3);
  REQUIRE(tokens[0] == "hello");
  REQUIRE(tokens[1] == "world");
  REQUIRE(tokens[2] == "foo");
}

TEST_CASE("tokenizer on empty string returns empty vector", "[generic]") {
  auto tokens = gen::tokenizer("");
  REQUIRE(tokens.empty());
}

TEST_CASE("tokenizerDouble parses doubles from string", "[generic]") {
  auto tokens = gen::tokenizerDouble("1.5 2.7 3.14");
  REQUIRE(tokens.size() == 3);
  REQUIRE_THAT(tokens[0], Catch::Matchers::WithinAbs(1.5, 1e-10));
  REQUIRE_THAT(tokens[2], Catch::Matchers::WithinAbs(3.14, 1e-10));
}

TEST_CASE("tokenizerInt parses ints from string", "[generic]") {
  auto tokens = gen::tokenizerInt("10 20 30");
  REQUIRE(tokens.size() == 3);
  REQUIRE(tokens[0] == 10);
  REQUIRE(tokens[2] == 30);
}

// -- calcMedian tests --

TEST_CASE("calcMedian of odd-length vector", "[generic]") {
  std::vector<double> v = {5.0, 1.0, 3.0};
  double m = gen::calcMedian(&v);
  REQUIRE_THAT(m, Catch::Matchers::WithinAbs(3.0, 1e-10));
}

TEST_CASE("calcMedian of even-length vector", "[generic]") {
  std::vector<double> v = {1.0, 2.0, 3.0, 4.0};
  double m = gen::calcMedian(&v);
  REQUIRE_THAT(m, Catch::Matchers::WithinAbs(2.5, 1e-10));
}

TEST_CASE("calcMedian of single element", "[generic]") {
  std::vector<double> v = {42.0};
  double m = gen::calcMedian(&v);
  REQUIRE_THAT(m, Catch::Matchers::WithinAbs(42.0, 1e-10));
}

// -- file_exists tests --

TEST_CASE("file_exists returns true for existing file", "[generic]") {
  // /dev/null always exists on unix
  REQUIRE(gen::file_exists("/dev/null"));
}

TEST_CASE("file_exists returns false for nonexistent file", "[generic]") {
  REQUIRE_FALSE(gen::file_exists("/tmp/this_file_should_not_exist_12345.xyz"));
}

// -- radDeg tests --

TEST_CASE("radDeg converts pi to 180 degrees", "[generic]") {
  double deg = gen::radDeg(gen::pi);
  REQUIRE_THAT(deg, Catch::Matchers::WithinAbs(180.0, 1e-10));
}

TEST_CASE("radDeg converts 0 to 0 degrees", "[generic]") {
  REQUIRE_THAT(gen::radDeg(0.0), Catch::Matchers::WithinAbs(0.0, 1e-10));
}

// -- avgVector tests --

TEST_CASE("avgVector divides each component by neigh count", "[generic]") {
  std::vector<std::complex<double>> v = {{4.0, 2.0}, {6.0, 0.0}, {0.0, 8.0}};
  auto avg = gen::avgVector(v, 1, 2); // l=1 -> 2*1+1=3 components, neigh=2

  REQUIRE_THAT(avg[0].real(), Catch::Matchers::WithinAbs(2.0, 1e-10));
  REQUIRE_THAT(avg[0].imag(), Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE_THAT(avg[1].real(), Catch::Matchers::WithinAbs(3.0, 1e-10));
  REQUIRE_THAT(avg[2].imag(), Catch::Matchers::WithinAbs(4.0, 1e-10));
}

TEST_CASE("avgVector with zero neighbours returns original", "[generic]") {
  std::vector<std::complex<double>> v = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  auto avg = gen::avgVector(v, 1, 0);

  REQUIRE_THAT(avg[0].real(), Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE_THAT(avg[2].imag(), Catch::Matchers::WithinAbs(6.0, 1e-10));
}

// -- compareByAtomID tests --

TEST_CASE("prettyPrintYoda does not crash", "[generic]") {
  auto cloud = makeTwoAtomCloud(1.0, 0.0, 0.0, 4.0, 0.0, 0.0);
  cloud.boxLow = {0.0, 0.0, 0.0};

  std::string outFile = std::filesystem::temp_directory_path().append("dseams_test_prettyprint.dat").string();
  int ret = gen::prettyPrintYoda(cloud, outFile);
  REQUIRE(ret == 0);

  std::filesystem::remove(outFile);
}

TEST_CASE("compareByAtomID sorts points by atomID", "[generic]") {
  molSys::Point<double> a, b;
  a.atomID = 5;
  b.atomID = 3;

  REQUIRE_FALSE(gen::compareByAtomID(a, b)); // 5 < 3 is false
  REQUIRE(gen::compareByAtomID(b, a));       // 3 < 5 is true
}

// -- getAverageWithoutOutliers regression tests --

TEST_CASE("getAverageWithoutOutliers does not crash when all values are outliers",
          "[generic]") {
  // Use an even-sized vector where all values are identical.
  // IQR=0, fences equal the value, so no outliers.
  std::vector<double> v = {1.0, 1.0, 1.0, 1.0};
  double avg = gen::getAverageWithoutOutliers(v);
  REQUIRE_THAT(avg, Catch::Matchers::WithinAbs(1.0, 1e-10));
}

TEST_CASE("getAverageWithoutOutliers handles odd length and empty input",
          "[generic]") {
  REQUIRE(gen::getAverageWithoutOutliers({}) == 0.0);
  REQUIRE_THAT(gen::getAverageWithoutOutliers({2.0}),
               Catch::Matchers::WithinAbs(2.0, 1e-10));
  REQUIRE_THAT(gen::getAverageWithoutOutliers({1.0, 2.0, 3.0}),
               Catch::Matchers::WithinAbs(2.0, 1e-10));
}

TEST_CASE("getAverageWithoutOutliers normal case with outlier excluded",
          "[generic]") {
  // {1, 2, 3, 4, 5, 100} -- 100 is an outlier
  // Sorted: {1,2,3,4,5,100}. Median=(3+4)/2=3.5
  // Lower half: {1,2,3}, Q1=2. Upper half: {4,5,100}, Q3=5
  // IQR=3, fences: 2-4.5=-2.5 and 5+4.5=9.5
  // So 100 excluded, avg of {1,2,3,4,5} = 3.0
  std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0, 100.0};
  double avg = gen::getAverageWithoutOutliers(v);
  REQUIRE_THAT(avg, Catch::Matchers::WithinAbs(3.0, 1e-10));
}

// -- eigenVecAngle regression tests --

TEST_CASE("eigenVecAngle does not produce NaN for nearly parallel vectors",
          "[generic]") {
  // Two nearly identical vectors -- dot/(|a|*|b|) may slightly exceed 1.0
  std::vector<double> a = {1.0, 0.0, 0.0};
  std::vector<double> b = {1.0, 1e-16, 0.0};
  double angle = gen::eigenVecAngle(a, b);
  REQUIRE_FALSE(std::isnan(angle));
  REQUIRE_THAT(angle, Catch::Matchers::WithinAbs(0.0, 1e-10));
}

TEST_CASE("eigenVecAngle does not produce NaN for antiparallel vectors",
          "[generic]") {
  std::vector<double> a = {1.0, 0.0, 0.0};
  std::vector<double> b = {-1.0, 0.0, 0.0};
  double angle = gen::eigenVecAngle(a, b);
  REQUIRE_FALSE(std::isnan(angle));
  REQUIRE_THAT(angle, Catch::Matchers::WithinAbs(gen::pi, 1e-10));
}
