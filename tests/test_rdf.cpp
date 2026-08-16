#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <rdf.hpp>

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
