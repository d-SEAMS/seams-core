#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <site.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

Cloud makeCloud() {
  Cloud cloud;
  cloud.box = {10.0, 10.0, 10.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  return cloud;
}

void addAtom(Cloud &cloud, int atomID, int molID, int type, double x, double y,
             double z) {
  molSys::Point<double> pt;
  pt.atomID = atomID;
  pt.molID = molID;
  pt.type = type;
  pt.x = x;
  pt.y = y;
  pt.z = z;
  cloud.idIndexMap[atomID] = static_cast<int>(cloud.pts.size());
  cloud.pts.push_back(pt);
  cloud.nop = static_cast<int>(cloud.pts.size());
}

} // namespace

TEST_CASE("of override wins over type", "[site]") {
  site::Table table;
  table.typeToKind[1] = site::Kind::cationHead;
  table.typeToKind[2] = site::Kind::anion;
  table.atomOverride[17] = site::Kind::tail;

  molSys::Point<double> head;
  head.type = 1;
  head.atomID = 3;
  REQUIRE(table.of(head) == site::Kind::cationHead);
  REQUIRE(table.ofType(1) == site::Kind::cationHead);

  molSys::Point<double> overridden;
  overridden.type = 1;
  overridden.atomID = 17;
  REQUIRE(table.of(overridden) == site::Kind::tail);
  REQUIRE(table.ofType(1) == site::Kind::cationHead);

  molSys::Point<double> unknown;
  unknown.type = 9;
  unknown.atomID = 99;
  REQUIRE(table.of(unknown) == site::Kind::unspecified);
  REQUIRE(table.ofType(9) == site::Kind::unspecified);
}

TEST_CASE("type 1 is not chemistry", "[site]") {
  site::Table table;
  REQUIRE(table.family == site::Family::waterIce);
  REQUIRE(table.ofType(1) == site::Kind::unspecified);

  molSys::Point<double> pt;
  pt.type = 1;
  pt.atomID = 1;
  REQUIRE(table.of(pt) == site::Kind::unspecified);

  table.family = site::Family::ionicLiquid;
  REQUIRE(table.of(pt) == site::Kind::unspecified);

  table.typeToKind[1] = site::Kind::waterO;
  REQUIRE(table.of(pt) == site::Kind::waterO);

  table.typeToKind[1] = site::Kind::cationHead;
  REQUIRE(table.of(pt) == site::Kind::cationHead);
}

TEST_CASE("indicesOf polar is the cationHead anion polar union", "[site]") {
  site::Table table;
  table.typeToKind[1] = site::Kind::cationHead;
  table.typeToKind[2] = site::Kind::anion;
  table.typeToKind[3] = site::Kind::tail;
  table.typeToKind[4] = site::Kind::polar;
  table.typeToKind[5] = site::Kind::apolar;

  auto cloud = makeCloud();
  addAtom(cloud, 1, 1, 1, 1.0, 0.0, 0.0);
  addAtom(cloud, 2, 2, 2, 2.0, 0.0, 0.0);
  addAtom(cloud, 3, 1, 3, 3.0, 0.0, 0.0);
  addAtom(cloud, 4, 1, 3, 4.0, 0.0, 0.0);
  addAtom(cloud, 5, 3, 4, 5.0, 0.0, 0.0);
  addAtom(cloud, 6, 4, 5, 6.0, 0.0, 0.0);

  const auto polar = site::indicesOf(cloud, table, site::Kind::polar);
  REQUIRE(polar.size() == 3);
  REQUIRE(polar[0] == 0);
  REQUIRE(polar[1] == 1);
  REQUIRE(polar[2] == 4);

  const auto apolar = site::indicesOf(cloud, table, site::Kind::apolar);
  REQUIRE(apolar.size() == 3);
  REQUIRE(apolar[0] == 2);
  REQUIRE(apolar[1] == 3);
  REQUIRE(apolar[2] == 5);

  REQUIRE(site::indicesOf(cloud, table, site::Kind::cationHead).size() == 1);
  REQUIRE(site::indicesOf(cloud, table, site::Kind::tail).size() == 2);
}

TEST_CASE("ionCloud two mols cation heads plus anion 1", "[site]") {
  site::Table table;
  table.typeToKind[1] = site::Kind::cationHead;
  table.typeToKind[2] = site::Kind::anion;
  table.typeToKind[3] = site::Kind::tail;
  table.typeToKind[4] = site::Kind::donorH;

  auto cloud = makeCloud();
  // Cation: two head atoms split across the periodic seam, plus a tail
  // and a donorH that must not enter the COM.
  addAtom(cloud, 10, 7, 1, 0.5, 1.0, 2.0);
  addAtom(cloud, 11, 7, 1, 9.5, 1.0, 2.0);
  addAtom(cloud, 12, 7, 3, 5.0, 4.0, 5.0);
  addAtom(cloud, 13, 7, 4, 0.6, 1.1, 2.1);
  // Monatomic anion.
  addAtom(cloud, 20, 8, 2, 4.0, 5.0, 6.0);

  const auto ions = site::ionCloud(cloud, table);
  REQUIRE(ions.nop == 2);
  REQUIRE(ions.box == cloud.box);
  REQUIRE(ions.boxLow == cloud.boxLow);
  REQUIRE(ions.currentFrame == cloud.currentFrame);

  REQUIRE(ions.pts[0].type == 1);
  REQUIRE(ions.pts[0].molID == 7);
  REQUIRE(ions.pts[1].type == 2);
  REQUIRE(ions.pts[1].molID == 8);
  REQUIRE_THAT(ions.pts[1].x, Catch::Matchers::WithinAbs(4.0, 1e-12));
  REQUIRE_THAT(ions.pts[1].y, Catch::Matchers::WithinAbs(5.0, 1e-12));
  REQUIRE_THAT(ions.pts[1].z, Catch::Matchers::WithinAbs(6.0, 1e-12));

  // First head at 0.5; second unwraps across L=10 to -0.5; COM x = 0.
  // Tail and donorH are not averaged in.
  REQUIRE_THAT(ions.pts[0].x, Catch::Matchers::WithinAbs(0.0, 1e-12));
  REQUIRE_THAT(ions.pts[0].y, Catch::Matchers::WithinAbs(1.0, 1e-12));
  REQUIRE_THAT(ions.pts[0].z, Catch::Matchers::WithinAbs(2.0, 1e-12));

  const auto dr = gen::relDist(cloud, 1, 0);
  const double handX = 0.5 * (cloud.pts[0].x + cloud.pts[0].x + dr[0]);
  const double handY = 0.5 * (cloud.pts[0].y + cloud.pts[0].y + dr[1]);
  const double handZ = 0.5 * (cloud.pts[0].z + cloud.pts[0].z + dr[2]);
  REQUIRE_THAT(ions.pts[0].x, Catch::Matchers::WithinAbs(handX, 1e-12));
  REQUIRE_THAT(ions.pts[0].y, Catch::Matchers::WithinAbs(handY, 1e-12));
  REQUIRE_THAT(ions.pts[0].z, Catch::Matchers::WithinAbs(handZ, 1e-12));
}

TEST_CASE("ionCloud sheared box COM uses relDist not bound-span wrap",
          "[site]") {
  site::Table table;
  table.typeToKind[1] = site::Kind::cationHead;

  Cloud cloud;
  cloud.box = {15.0, 8.660254037844386, 10.0, 5.0, 0.0, 0.0};
  cloud.boxLow = {0.0, 0.0, 0.0};
  cloud.currentFrame = 1;
  addAtom(cloud, 1, 1, 1, 0.5, 0.5, 1.0);
  addAtom(cloud, 2, 1, 1, 1.0, 8.0, 1.0);

  const auto ions = site::ionCloud(cloud, table);
  REQUIRE(ions.nop == 1);
  REQUIRE(ions.pts[0].type == 1);
  REQUIRE(ions.box == cloud.box);
  REQUIRE(ions.boxLow == cloud.boxLow);

  const auto dr = gen::relDist(cloud, 1, 0);
  const double expectX = 0.5 * (cloud.pts[0].x + cloud.pts[0].x + dr[0]);
  const double expectY = 0.5 * (cloud.pts[0].y + cloud.pts[0].y + dr[1]);
  const double expectZ = 0.5 * (cloud.pts[0].z + cloud.pts[0].z + dr[2]);
  REQUIRE_THAT(ions.pts[0].x, Catch::Matchers::WithinAbs(expectX, 1e-10));
  REQUIRE_THAT(ions.pts[0].y, Catch::Matchers::WithinAbs(expectY, 1e-10));
  REQUIRE_THAT(ions.pts[0].z, Catch::Matchers::WithinAbs(expectZ, 1e-12));

  std::array<double, 3> span = {cloud.pts[1].x - cloud.pts[0].x,
                                cloud.pts[1].y - cloud.pts[0].y,
                                cloud.pts[1].z - cloud.pts[0].z};
  for (int k = 0; k < 3; k++) {
    span[k] -= cloud.box[static_cast<std::size_t>(k)] *
               std::round(span[k] / cloud.box[static_cast<std::size_t>(k)]);
  }
  const double spanX = 0.5 * (cloud.pts[0].x + cloud.pts[0].x + span[0]);
  const double spanY = 0.5 * (cloud.pts[0].y + cloud.pts[0].y + span[1]);
  REQUIRE(std::abs(ions.pts[0].x - spanX) > 1.0);
  REQUIRE(std::abs(expectX - spanX) > 1.0);
  (void)spanY;
}

TEST_CASE("ionCloud monatomic anion is a copy", "[site]") {
  site::Table table;
  table.typeToKind[2] = site::Kind::anion;

  auto cloud = makeCloud();
  addAtom(cloud, 20, 8, 2, 4.0, 5.0, 6.0);
  cloud.pts[0].inSlice = false;

  const auto ions = site::ionCloud(cloud, table);
  REQUIRE(ions.nop == 1);
  REQUIRE(ions.pts[0].type == 2);
  REQUIRE(ions.pts[0].atomID == 20);
  REQUIRE(ions.pts[0].molID == 8);
  REQUIRE_THAT(ions.pts[0].x, Catch::Matchers::WithinAbs(4.0, 1e-12));
  REQUIRE_THAT(ions.pts[0].y, Catch::Matchers::WithinAbs(5.0, 1e-12));
  REQUIRE_THAT(ions.pts[0].z, Catch::Matchers::WithinAbs(6.0, 1e-12));
  REQUIRE(ions.pts[0].inSlice == false);
  REQUIRE(ions.pts[0].iceType == cloud.pts[0].iceType);
}

TEST_CASE("parseSiteSpec reads type=kind and optional family", "[site]") {
  const auto table =
      site::parseSiteSpec("1=cationHead, 2=anion, 3=tail, family=ionicLiquid");
  REQUIRE(table.family == site::Family::ionicLiquid);
  REQUIRE(table.ofType(1) == site::Kind::cationHead);
  REQUIRE(table.ofType(2) == site::Kind::anion);
  REQUIRE(table.ofType(3) == site::Kind::tail);
  REQUIRE(table.ofType(4) == site::Kind::unspecified);
  REQUIRE(site::lammpsTypeOfKind(table, site::Kind::anion) == 2);
}

TEST_CASE("parseSiteSpec rejects unknown kind names", "[site]") {
  REQUIRE_THROWS_AS(site::parseSiteSpec("1=oxygen"), std::invalid_argument);
}

TEST_CASE("lammpsTypeOfKind errors when the kind is not unique", "[site]") {
  site::Table table;
  table.typeToKind[1] = site::Kind::cationHead;
  table.typeToKind[2] = site::Kind::anion;
  table.typeToKind[3] = site::Kind::cationHead;
  REQUIRE(site::lammpsTypeOfKind(table, site::Kind::anion) == 2);
  REQUIRE_THROWS_AS(site::lammpsTypeOfKind(table, site::Kind::cationHead),
                    std::runtime_error);
  REQUIRE_THROWS_AS(site::lammpsTypeOfKind(table, site::Kind::tail),
                    std::runtime_error);
}
