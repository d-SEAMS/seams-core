#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <density.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <site.hpp>

#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

void addAtom(Cloud &cloud, int atomID, int type, double x, double y, double z) {
  molSys::Point<double> pt;
  pt.atomID = atomID;
  pt.molID = atomID;
  pt.type = type;
  pt.x = x;
  pt.y = y;
  pt.z = z;
  cloud.idIndexMap[atomID] = static_cast<int>(cloud.pts.size());
  cloud.pts.push_back(pt);
  cloud.nop = static_cast<int>(cloud.pts.size());
}

Cloud uniformZCloud(const std::vector<double> &box, int n, int type) {
  Cloud cloud;
  cloud.box = box;
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  const double lz = box.size() > 2 ? box[2] : 0.0;
  for (int i = 0; i < n; ++i) {
    const double z = (static_cast<double>(i) + 0.5) * lz / static_cast<double>(n);
    addAtom(cloud, i + 1, type, 1.0, 1.0, z);
  }
  return cloud;
}

double integrateRho(const site::DensityZ &d, double area, double dz) {
  double acc = 0.0;
  for (double rho : d.rho) {
    acc += rho * area * dz;
  }
  return acc;
}

} // namespace

TEST_CASE("uniform z histogram integrates to N on an ortho box", "[density]") {
  auto cloud = uniformZCloud({10.0, 8.0, 20.0}, 10, 1);
  const int nbin = 20;
  const auto d = site::densityZ(cloud, 0, nbin, 2);
  REQUIRE(d.z.size() == static_cast<std::size_t>(nbin));
  REQUIRE(d.rho.size() == static_cast<std::size_t>(nbin));
  REQUIRE(d.type == 0);
  const double area = 10.0 * 8.0;
  const double dz = 20.0 / static_cast<double>(nbin);
  REQUIRE_THAT(integrateRho(d, area, dz), Catch::Matchers::WithinAbs(10.0, 1e-9));
  REQUIRE_THAT(d.z.front(), Catch::Matchers::WithinAbs(0.5 * dz, 1e-12));
}

TEST_CASE("tilted box uses det H face area not bound spans", "[density]") {
  auto cloud =
      uniformZCloud({15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0}, 10, 1);
  const double vol = nneigh::dumpVolume(cloud);
  double lengths[3] = {0.0, 0.0, 0.0};
  nneigh::dumpCellLengths(cloud.box, cloud.boxLow, lengths);
  const double area = vol / lengths[2];
  const double boundArea = 15.0 * 8.660254037844386;
  REQUIRE(area < boundArea - 1.0);
  REQUIRE_THAT(area, Catch::Matchers::WithinAbs(lengths[0] * lengths[1], 1e-9));

  const int nbin = 10;
  const auto d = site::densityZ(cloud, 0, nbin, 2);
  const double dz = 10.0 / static_cast<double>(nbin);
  REQUIRE_THAT(integrateRho(d, area, dz), Catch::Matchers::WithinAbs(10.0, 1e-9));
  REQUIRE_THAT(integrateRho(d, boundArea, dz),
               Catch::Matchers::WithinAbs(10.0 * boundArea / area, 1e-9));
}

TEST_CASE("type 0 is every atom and type I filters", "[density]") {
  Cloud cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  addAtom(cloud, 1, 1, 1.0, 1.0, 1.0);
  addAtom(cloud, 2, 1, 1.0, 1.0, 3.0);
  addAtom(cloud, 3, 2, 1.0, 1.0, 7.0);

  const auto all = site::densityZ(cloud, 0, 10, 2);
  REQUIRE_THAT(integrateRho(all, 100.0, 1.0), Catch::Matchers::WithinAbs(3.0, 1e-9));

  const auto t1 = site::densityZ(cloud, 1, 10, 2);
  REQUIRE(t1.type == 1);
  REQUIRE_THAT(integrateRho(t1, 100.0, 1.0), Catch::Matchers::WithinAbs(2.0, 1e-9));

  const auto t2 = site::densityZ(cloud, 2, 10, 2);
  REQUIRE(t2.type == 2);
  REQUIRE_THAT(integrateRho(t2, 100.0, 1.0), Catch::Matchers::WithinAbs(1.0, 1e-9));
}

TEST_CASE("kind histogram uses indicesOf", "[density]") {
  Cloud cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  addAtom(cloud, 1, 1, 1.0, 1.0, 2.0);
  addAtom(cloud, 2, 2, 1.0, 1.0, 8.0);
  site::Table table;
  table.typeToKind[1] = site::Kind::cationHead;
  table.typeToKind[2] = site::Kind::anion;

  const auto heads = site::densityZ(cloud, table, site::Kind::cationHead, 10, 2);
  REQUIRE_THAT(integrateRho(heads, 100.0, 1.0),
               Catch::Matchers::WithinAbs(1.0, 1e-9));
}
