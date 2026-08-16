#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf.hpp>

#include <array>
#include <vector>

namespace {

molSys::PointCloud<molSys::Point<double>, double>
twoTypeCloud(const std::vector<double> &box,
             const double coords[2][3]) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = box;
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = 2;
  cloud.currentFrame = 1;
  const int types[2] = {1, 2};
  for (int i = 0; i < 2; i++) {
    molSys::Point<double> pt;
    pt.type = types[i];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[i][0];
    pt.y = coords[i][1];
    pt.z = coords[i][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  return cloud;
}

molSys::PointCloud<molSys::Point<double>, double>
typedCloud(const std::vector<double> &box,
           const std::vector<int> &types,
           const std::vector<std::array<double, 3>> &coords) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = box;
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.nop = static_cast<int>(coords.size());
  cloud.currentFrame = 1;
  for (int i = 0; i < cloud.nop; ++i) {
    molSys::Point<double> pt;
    pt.type = types[static_cast<std::size_t>(i)];
    pt.atomID = i + 1;
    pt.molID = i + 1;
    pt.x = coords[static_cast<std::size_t>(i)][0];
    pt.y = coords[static_cast<std::size_t>(i)][1];
    pt.z = coords[static_cast<std::size_t>(i)][2];
    cloud.pts.push_back(pt);
    cloud.idIndexMap[i + 1] = i;
  }
  return cloud;
}

double meanTypeIDegree(const std::vector<std::vector<int>> &nList,
                       const molSys::PointCloud<molSys::Point<double>, double> &cloud,
                       int typeI) {
  int nI = 0;
  int deg = 0;
  for (int i = 0; i < cloud.nop; ++i) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != typeI) {
      continue;
    }
    ++nI;
    if (i < static_cast<int>(nList.size()) && nList[static_cast<std::size_t>(i)].size() > 1) {
      deg += static_cast<int>(nList[static_cast<std::size_t>(i)].size()) - 1;
    }
  }
  return nI > 0 ? static_cast<double>(deg) / static_cast<double>(nI) : 0.0;
}

} // namespace

TEST_CASE("dumpVolume uses det H not bound spans", "[rdf]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  const double vol = nneigh::dumpVolume(cloud);
  const double lxLyLz = 10.0 * 8.660254037844386 * 10.0;
  const double bound = 15.0 * 8.660254037844386 * 10.0;
  REQUIRE_THAT(vol, Catch::Matchers::WithinAbs(lxLyLz, 1e-9));
  REQUIRE(vol < bound - 1.0);
}

TEST_CASE("partial RDF first bin sees the sheared a-image pair", "[rdf]") {
  const double coords[2][3] = {{0.2, 0.1, 1.0}, {9.7, 0.1, 1.0}};
  auto cloud =
      twoTypeCloud({15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0}, coords);
  REQUIRE_THAT(gen::periodicDistSq(cloud, 0, 1),
               Catch::Matchers::WithinAbs(0.25, 1e-9));
  const auto gr = rdf::partialRdf(cloud, 1, 2, 1.0, 1);
  REQUIRE(gr.count.size() == 1);
  REQUIRE(gr.count[0] == 1);
  REQUIRE(gr.nI == 1);
  REQUIRE(gr.nJ == 1);
  REQUIRE(gr.g[0] > 0.0);
}

TEST_CASE("partial RDF like-type does not count the unlike neighbour",
          "[rdf]") {
  const double coords[2][3] = {{0.0, 0.0, 0.0}, {5.0, 0.0, 0.0}};
  auto cloud = twoTypeCloud({10.0, 10.0, 10.0}, coords);
  const auto like = rdf::partialRdf(cloud, 1, 1, 6.0, 6);
  int total = 0;
  for (int c : like.count) {
    total += c;
  }
  REQUIRE(total == 0);
  const auto unlike = rdf::partialRdf(cloud, 1, 2, 6.0, 6);
  REQUIRE(unlike.count[5] == 1);
}

TEST_CASE("runningCN at rmax is one for a single I-J pair", "[rdf]") {
  const double coords[2][3] = {{1.0, 1.0, 1.0}, {3.0, 1.0, 1.0}};
  auto cloud = twoTypeCloud({10.0, 10.0, 10.0}, coords);
  const auto gr = rdf::partialRdf(cloud, 1, 2, 5.0, 10);
  REQUIRE(gr.nI == 1);
  REQUIRE(gr.nJ == 1);
  REQUIRE(gr.binwidth == 0.5);
  REQUIRE(gr.count[4] == 1);
  const auto cn = rdf::runningCN(gr);
  REQUIRE(cn.size() == gr.g.size());
  REQUIRE_THAT(cn[3], Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(cn.back(), Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(rdf::coordinationNumber(gr, 5.0),
               Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(rdf::coordinationNumber(gr, 1.5),
               Catch::Matchers::WithinAbs(0.0, 1e-12));
}

TEST_CASE("neighList degree matches site-site CN past every pair and before none",
          "[rdf]") {
  auto cloud = typedCloud(
      {20.0, 20.0, 20.0}, {1, 1, 2, 2},
      {{{0.0, 0.0, 0.0}, {10.0, 10.0, 10.0}, {2.0, 0.0, 0.0}, {10.0, 12.0, 10.0}}});
  const auto gr = rdf::partialRdf(cloud, 1, 2, 10.0, 10);
  const double rhoJ = static_cast<double>(gr.nJ) / gr.volume;

  const auto past = nneigh::neighList(3.0, cloud, 1, 2);
  const double degPast = meanTypeIDegree(past, cloud, 1);
  REQUIRE_THAT(degPast, Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(rdf::coordinationNumber(gr, 3.0, rhoJ),
               Catch::Matchers::WithinAbs(degPast, 1e-9));

  const auto none = nneigh::neighList(1.0, cloud, 1, 2);
  const double degNone = meanTypeIDegree(none, cloud, 1);
  REQUIRE_THAT(degNone, Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(rdf::coordinationNumber(gr, 1.0, rhoJ),
               Catch::Matchers::WithinAbs(degNone, 1e-9));
}

TEST_CASE("firstMinimumBin finds the valley after the first of two peaks",
          "[rdf]") {
  rdf::PartialRdf h;
  h.binwidth = 0.5;
  h.g = {0.1, 0.4, 2.0, 1.1, 0.3, 0.5, 1.6, 0.8};
  REQUIRE(rdf::firstMinimumBin(h) == 4);

  rdf::PartialRdf empty;
  REQUIRE(rdf::firstMinimumBin(empty) == -1);

  rdf::PartialRdf rising;
  rising.g = {0.1, 0.2, 0.4, 0.8, 1.2};
  REQUIRE(rdf::firstMinimumBin(rising) == -1);
}
