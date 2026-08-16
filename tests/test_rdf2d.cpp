#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <rdf2d.hpp>

#include <cmath>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

// Helper: build a small PointCloud with atoms at known positions
static molSys::PointCloud<molSys::Point<double>, double>
makeRdfCloud(int nop, double boxLen = 10.0) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {boxLen, boxLen, boxLen};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  // Place atoms in a grid-like pattern within the box
  for (int i = 0; i < nop; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i;
    pt.x = 1.0 + 2.0 * (i % 4);
    pt.y = 1.0 + 2.0 * ((i / 4) % 4);
    pt.z = 0.5; // quasi-2D: thin slab in z
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = nop;
  return cloud;
}

// -- getSystemLengths tests --

TEST_CASE("getSystemLengths returns correct extents", "[rdf2d]") {
  auto cloud = makeRdfCloud(4);
  // Atoms at x={1,3,5,7}, y={1,1,1,1}, z={0.5,0.5,0.5,0.5}
  auto lengths = rdf2::getSystemLengths(cloud);

  REQUIRE(lengths.size() == 3);
  REQUIRE_THAT(lengths[0], Catch::Matchers::WithinAbs(6.0, 1e-10)); // x: 7-1
  REQUIRE_THAT(lengths[1], Catch::Matchers::WithinAbs(0.0, 1e-10)); // y: all 1.0
  REQUIRE_THAT(lengths[2], Catch::Matchers::WithinAbs(0.0, 1e-10)); // z: all 0.5
}

TEST_CASE("getSystemLengths with spread atoms", "[rdf2d]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  double coords[3][3] = {{1.0, 2.0, 3.0}, {5.0, 8.0, 1.0}, {3.0, 4.0, 7.0}};
  for (int i = 0; i < 3; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i;
    pt.molID = i;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i] = i;
  }
  cloud.nop = 3;

  auto lengths = rdf2::getSystemLengths(cloud);
  REQUIRE_THAT(lengths[0], Catch::Matchers::WithinAbs(4.0, 1e-10)); // x: 5-1
  REQUIRE_THAT(lengths[1], Catch::Matchers::WithinAbs(6.0, 1e-10)); // y: 8-2
  REQUIRE_THAT(lengths[2], Catch::Matchers::WithinAbs(6.0, 1e-10)); // z: 7-1
}

// -- getPlaneArea tests --

TEST_CASE("getPlaneArea returns product of two largest dimensions", "[rdf2d]") {
  std::vector<double> lengths = {2.0, 5.0, 3.0};
  double area = rdf2::getPlaneArea(lengths);
  // Sorted desc: {5,3,2}, area = 5*3 = 15
  REQUIRE_THAT(area, Catch::Matchers::WithinAbs(15.0, 1e-10));
}

TEST_CASE("getPlaneArea with equal dimensions", "[rdf2d]") {
  std::vector<double> lengths = {4.0, 4.0, 4.0};
  double area = rdf2::getPlaneArea(lengths);
  REQUIRE_THAT(area, Catch::Matchers::WithinAbs(16.0, 1e-10));
}

// -- sampleRDF_AA tests --

TEST_CASE("sampleRDF_AA produces histogram with correct bin count", "[rdf2d]") {
  // Two atoms separated by 2.0
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;

  molSys::Point<double> p0, p1;
  p0.type = 1; p0.atomID = 0; p0.molID = 0;
  p0.x = 0.0; p0.y = 0.0; p0.z = 0.0;
  p1.type = 1; p1.atomID = 1; p1.molID = 1;
  p1.x = 2.0; p1.y = 0.0; p1.z = 0.0;

  cloud.pts.push_back(p0);
  cloud.pts.push_back(p1);
  cloud.nop = 2;
  cloud.idIndexMap[0] = 0;
  cloud.idIndexMap[1] = 1;

  double cutoff = 5.0;
  double binwidth = 1.0;
  int nbin = static_cast<int>(cutoff / binwidth);

  auto hist = rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);

  REQUIRE(hist.size() == static_cast<size_t>(nbin));
  // Distance is 2.0, falls in bin 2 (ibin = int(2.0/1.0) = 2)
  REQUIRE(hist[2] == 2); // +2 for the pair (iatom and jatom)
  // Other bins should be 0
  REQUIRE(hist[0] == 0);
  REQUIRE(hist[1] == 0);
}

TEST_CASE("sampleRDF_AA histograms the tilt a-image pair", "[rdf2d]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  const double coords[2][3] = {{0.2, 0.1, 1.0}, {9.7, 0.1, 1.0}};
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = 1;
    pt.atomID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  REQUIRE_THAT(gen::periodicDistSq(cloud, 0, 1),
               Catch::Matchers::WithinAbs(0.25, 1e-9));
  const double cutoff = 1.0;
  const double binwidth = 0.1;
  const int nbin = 10;
  auto hist = rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);
  REQUIRE(hist.size() == static_cast<std::size_t>(nbin));
  REQUIRE(hist[5] == 2);
}

// -- normalizeRDF tests --

TEST_CASE("normalizeRDF produces non-negative values", "[rdf2d]") {
  int nopA = 10;
  std::vector<double> rdfValues(5, 0.0);
  std::vector<int> histogram = {0, 4, 8, 2, 0};
  double binwidth = 1.0;
  int nbin = 5;
  std::vector<double> volumeLengths = {10.0, 10.0, 2.0};
  int nIter = 1;

  int ret =
      rdf2::normalizeRDF(nopA, rdfValues, histogram, binwidth, nbin,
                         volumeLengths, nIter);

  REQUIRE(ret == 0);
  for (int i = 0; i < nbin; i++) {
    REQUIRE(rdfValues[i] >= 0.0);
  }
}

// -- rdf2Danalysis_AA integration test --

TEST_CASE("rdf2Danalysis_AA runs without error for single frame", "[rdf2d]") {
  auto cloud = makeRdfCloud(8, 20.0);
  cloud.currentFrame = 1;

  std::string tmpPath = fs::temp_directory_path().append("dseams_test_rdf2d/").string();
  fs::create_directories(tmpPath);

  std::vector<double> rdfValues;
  double cutoff = 5.0;
  double binwidth = 0.5;

  // Single frame: firstFrame == finalFrame == currentFrame
  int ret = rdf2::rdf2Danalysis_AA(tmpPath, rdfValues, cloud, cutoff, binwidth,
                                    1, 1);

  REQUIRE(ret == 0);
  REQUIRE_FALSE(rdfValues.empty());

  // Check that RDF file was written
  REQUIRE(fs::exists(tmpPath + "topoMonolayer/rdf.dat"));

  // Cleanup
  std::error_code _ec_; fs::remove_all(tmpPath, _ec_);
}
