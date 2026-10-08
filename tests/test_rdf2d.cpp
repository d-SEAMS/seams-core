#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <random>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf2d.hpp>

#include <algorithm>
#include <array>
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

// Every pair once, at gen::periodicDist, binned as sampleRDF_AA bins.
static std::vector<int>
referenceHistogram(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                   double cutoff, double binwidth, int nbin) {
  std::vector<int> reference(static_cast<std::size_t>(nbin), 0);
  const gen::FracBox frame = gen::makeFracBox(cloud);
  for (int i = 0; i < cloud.nop; i++) {
    for (int j = i + 1; j < cloud.nop; j++) {
      const double r = std::sqrt(gen::periodicDistSq(frame, cloud, i, j));
      if (r > cutoff) {
        continue;
      }
      const int b = std::min(static_cast<int>(r / binwidth), nbin - 1);
      reference[static_cast<std::size_t>(b)] += 2;
    }
  }
  return reference;
}

static void addPoint(molSys::PointCloud<molSys::Point<double>, double> &cloud,
                     double x, double y, double z) {
  molSys::Point<double> pt;
  pt.type = 1;
  pt.atomID = cloud.nop + 1;
  pt.x = x;
  pt.y = y;
  pt.z = z;
  cloud.pts.push_back(pt);
  cloud.idIndexMap[pt.atomID] = cloud.nop;
  cloud.nop++;
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

TEST_CASE("getSystemLengths hex-prism uses dump H not bound spans",
          "[rdf2d]") {
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
  auto lengths = rdf2::getSystemLengths(cloud);
  REQUIRE(lengths.size() == 3);
  REQUIRE_THAT(lengths[0], Catch::Matchers::WithinAbs(10.0, 1e-9));
  REQUIRE_THAT(lengths[1], Catch::Matchers::WithinAbs(8.660254037844386, 1e-9));
  REQUIRE_THAT(lengths[2], Catch::Matchers::WithinAbs(10.0, 1e-9));
  const double vol = lengths[0] * lengths[1] * lengths[2];
  const double dumpVol = nneigh::dumpVolume(cloud);
  const double planeArea = lengths[0] * lengths[1];
  const double bound = 15.0 * 8.660254037844386 * 10.0;
  REQUIRE_THAT(vol, Catch::Matchers::WithinAbs(dumpVol, 1e-9));
  REQUIRE_THAT(planeArea,
               Catch::Matchers::WithinAbs(10.0 * 8.660254037844386, 1e-9));
  REQUIRE(vol < bound - 1.0);
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

TEST_CASE("sampleRDF_AA at half an edge is MIC-once", "[rdf2d]") {
  // Atoms sit half a box apart. Both the direct pair and the wrapped
  // image have length L/2. The histogram keeps one of them.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  const double coords[2][3] = {{0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}};
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
  REQUIRE_THAT(gen::periodicDist(cloud, 0, 1),
               Catch::Matchers::WithinAbs(5.0, 1e-12));
  const double cutoff = 5.0;
  const double binwidth = 1.0;
  const int nbin = 5;
  auto hist = rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);
  int total = 0;
  for (int bin : hist) {
    total += bin;
  }
  REQUIRE(total == 2);
}

TEST_CASE("sampleRDF_AA cutoff past L/2 matches MIC-once", "[rdf2d]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {5.0, 5.0, 5.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  molSys::Point<double> p0, p1;
  p0.type = 1;
  p0.atomID = 0;
  p0.molID = 0;
  p0.x = 0.0;
  p0.y = 0.0;
  p0.z = 0.0;
  p1.type = 1;
  p1.atomID = 1;
  p1.molID = 1;
  p1.x = 2.0;
  p1.y = 0.0;
  p1.z = 0.0;
  cloud.pts.push_back(p0);
  cloud.pts.push_back(p1);
  cloud.nop = 2;
  cloud.idIndexMap[0] = 0;
  cloud.idIndexMap[1] = 1;
  REQUIRE_THAT(gen::periodicDist(cloud, 0, 1),
               Catch::Matchers::WithinAbs(2.0, 1e-12));
  const double cutoff = 4.0;
  const double binwidth = 0.5;
  const int nbin = 8;
  auto hist = rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);
  REQUIRE(hist.size() == static_cast<std::size_t>(nbin));
  REQUIRE(hist[4] == 2);
  int extra = 0;
  for (int i = 0; i < nbin; i++) {
    if (i != 4) {
      extra += hist[static_cast<std::size_t>(i)];
    }
  }
  REQUIRE(extra == 0);
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
  const double r = gen::periodicDist(cloud, 0, 1);
  REQUIRE_THAT(r * r, Catch::Matchers::WithinAbs(0.25, 1e-9));
  const double cutoff = 1.0;
  const double binwidth = 0.1;
  const int nbin = 10;
  const int expectBin = static_cast<int>(r / binwidth);
  auto hist = rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);
  REQUIRE(hist.size() == static_cast<std::size_t>(nbin));
  REQUIRE(expectBin >= 0);
  REQUIRE(expectBin < nbin);
  REQUIRE(hist[static_cast<std::size_t>(expectBin)] == 2);
}

TEST_CASE("sampleRDF_AA tilt between edge and span is MIC-once", "[rdf2d]") {
  // xy = 60 and x span 61 recover lx = 1. The span half-minimum is 6;
  // the edge half-minimum is 0.5. Cutoff 2 sits between them.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {61.0, 12.0, 50.0, 60.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  const double coords[2][3] = {{0.1, 1.0, 1.0}, {0.9, 1.0, 1.0}};
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
  const double r = gen::periodicDist(cloud, 0, 1);
  REQUIRE_THAT(r, Catch::Matchers::WithinAbs(0.2, 1e-9));
  double lengths[3];
  nneigh::dumpCellLengths(cloud.box, cloud.boxLow, lengths);
  REQUIRE_THAT(lengths[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
  const double cutoff = 2.0;
  const double binwidth = 0.1;
  const int nbin = 20;
  const int expectBin = static_cast<int>(r / binwidth);
  auto hist = rdf2::sampleRDF_AA(cloud, cutoff, binwidth, nbin);
  REQUIRE(hist.size() == static_cast<std::size_t>(nbin));
  REQUIRE(expectBin >= 0);
  REQUIRE(expectBin < nbin);
  int total = 0;
  for (int i = 0; i < nbin; i++) {
    total += hist[static_cast<std::size_t>(i)];
  }
  // One MIC pair, stored once for each order.
  REQUIRE(total == 2);
  REQUIRE(hist[static_cast<std::size_t>(expectBin)] == 2);
}

TEST_CASE("sampleRDF_AA cell pairs match the direct minimum image",
          "[rdf2d]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {20.0, 20.0, 20.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 0;
  std::mt19937 rng(20);
  std::uniform_real_distribution<double> u(0.0, 20.0);
  for (int i = 0; i < 64; i++) {
    const double x = u(rng);
    const double y = u(rng);
    const double z = u(rng);
    addPoint(cloud, x, y, z);
  }
  const auto reference = referenceHistogram(cloud, 3.0, 0.1, 30);
  REQUIRE(rdf2::sampleRDF_AA(cloud, 3.0, 0.1, 30) == reference);
}

TEST_CASE("sampleRDF_AA cell pairs match a sheared minimum image",
          "[rdf2d]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {30.0, 20.0, 25.0, 4.0, 0.5, -0.3};
  cloud.boxLow = {1.0, -2.0, 0.5};
  cloud.nop = 0;
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> ux(1.0, 31.0);
  std::uniform_real_distribution<double> uy(-2.0, 18.0);
  std::uniform_real_distribution<double> uz(0.5, 25.5);
  for (int i = 0; i < 80; i++) {
    const double x = ux(rng);
    const double y = uy(rng);
    const double z = uz(rng);
    addPoint(cloud, x, y, z);
  }
  const auto reference = referenceHistogram(cloud, 3.0, 0.1, 30);
  REQUIRE(reference != std::vector<int>(30, 0));
  REQUIRE(rdf2::sampleRDF_AA(cloud, 3.0, 0.1, 30) == reference);
}

TEST_CASE("sampleRDF_AA keeps a tilted pair the fractional wrap pushes out",
          "[rdf2d]") {
  // yz = ly/2: the width across b is 10/sqrt(1.25), so a 4.9 cutoff passes
  // half of it while staying below half of every H diagonal.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {10.0, 15.0, 10.0, 0.0, 0.0, 5.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 0;
  addPoint(cloud, 5.0, 5.0, 2.0);
  addPoint(cloud, 5.0, 1.4, 5.0);
  REQUIRE_THAT(gen::periodicDist(cloud, 0, 1),
               Catch::Matchers::WithinAbs(std::sqrt(3.6 * 3.6 + 9.0), 1e-12));
  const auto reference = referenceHistogram(cloud, 4.9, 0.1, 49);
  REQUIRE(rdf2::sampleRDF_AA(cloud, 4.9, 0.1, 49) == reference);
}

TEST_CASE("sampleRDF_AA cell pairs size cells on the face separation",
          "[rdf2d]") {
  // yz = ly/2 in a 42 cell: four rows along the b diagonal are 9.4 apart
  // across the tilt, so a pair inside a 10 cutoff can sit two rows apart.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {42.0, 63.0, 42.0, 0.0, 0.0, 21.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 0;
  addPoint(cloud, 21.0, 42.0 * 0.249 + 21.0 * 0.6, 42.0 * 0.6);
  addPoint(cloud, 21.0, 42.0 * 0.51 + 21.0 * 0.4956, 42.0 * 0.4956);
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < 200; i++) {
    const double sx = u(rng);
    const double sy = u(rng);
    const double sz = u(rng);
    addPoint(cloud, 42.0 * sx, 42.0 * sy + 21.0 * sz, 42.0 * sz);
  }
  REQUIRE(gen::periodicDist(cloud, 0, 1) < 10.0);
  const auto reference = referenceHistogram(cloud, 10.0, 0.25, 40);
  REQUIRE(rdf2::sampleRDF_AA(cloud, 10.0, 0.25, 40) == reference);
}

TEST_CASE("sampleRDF_AA direct loop past the ball matches the minimum image",
          "[rdf2d]") {
  // The cutoff passes half the narrowest face separation, so the direct loop
  // runs; under tilt the pairs past the Smith ball take the Euclidean image.
  for (const auto &tilt : {std::array<double, 3>{0.0, 0.0, 0.0},
                           std::array<double, 3>{6.0, -4.0, 5.0}}) {
    const double L = 20.0;
    const double xy = tilt[0], xz = tilt[1], yz = tilt[2];
    const double xmin = std::min({0.0, xy, xz, xy + xz});
    const double xmax = std::max({0.0, xy, xz, xy + xz});
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud.box = {L + xmax - xmin, L + std::max(0.0, yz) - std::min(0.0, yz), L,
                 xy, xz, yz};
    cloud.boxLow = {xmin, std::min(0.0, yz), 0.0};
    if (xy == 0.0 && xz == 0.0 && yz == 0.0) {
      cloud.box.resize(3);
    }
    cloud.nop = 0;
    std::mt19937 rng(13);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int i = 0; i < 600; i++) {
      const double sx = u(rng);
      const double sy = u(rng);
      const double sz = u(rng);
      addPoint(cloud, L * sx + xy * sy + xz * sz, L * sy + yz * sz, L * sz);
    }
    const auto reference = referenceHistogram(cloud, 13.0, 0.1, 130);
    REQUIRE(rdf2::sampleRDF_AA(cloud, 13.0, 0.1, 130) == reference);
  }
}

TEST_CASE("sampleRDF_AA cell pairs match a threaded frame tilted on every axis",
          "[rdf2d]") {
  // Enough atoms for the parallel build and the threaded walk.
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {50.0, 45.0, 40.0, 6.0, -4.0, 5.0};
  cloud.boxLow = {-4.0, 0.0, 0.0};
  cloud.nop = 0;
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < 4500; i++) {
    const double sx = u(rng);
    const double sy = u(rng);
    const double sz = u(rng);
    addPoint(cloud, 40.0 * sx + 6.0 * sy - 4.0 * sz, 40.0 * sy + 5.0 * sz,
             40.0 * sz);
  }
  const auto reference = referenceHistogram(cloud, 5.0, 0.05, 100);
  REQUIRE(rdf2::sampleRDF_AA(cloud, 5.0, 0.05, 100) == reference);
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

TEST_CASE("normalizeRDF hex-prism density is dumpVolume not span product",
          "[rdf2d]") {
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
  auto lengths = rdf2::getSystemLengths(cloud);
  const double dumpVol = nneigh::dumpVolume(cloud);
  REQUIRE_THAT(lengths[0] * lengths[1] * lengths[2],
               Catch::Matchers::WithinAbs(dumpVol, 1e-9));
  REQUIRE_THAT(lengths[0] * lengths[1],
               Catch::Matchers::WithinAbs(10.0 * 8.660254037844386, 1e-9));

  const int nopA = 2;
  const int nbin = 5;
  const double binwidth = 1.0;
  const int nIter = 1;
  std::vector<int> histogram = {2, 0, 0, 0, 0};
  std::vector<double> rdfDump(static_cast<std::size_t>(nbin), 0.0);
  std::vector<double> rdfSpan(static_cast<std::size_t>(nbin), 0.0);
  REQUIRE(rdf2::normalizeRDF(nopA, rdfDump, histogram, binwidth, nbin, lengths,
                             nIter) == 0);
  REQUIRE(rdf2::normalizeRDF(nopA, rdfSpan, histogram, binwidth, nbin,
                             {15.0, 8.660254037844386, 10.0}, nIter) == 0);
  REQUIRE(rdfDump[0] > 0.0);
  REQUIRE(rdfSpan[0] > 0.0);
  REQUIRE(rdfDump[0] < rdfSpan[0]);
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
