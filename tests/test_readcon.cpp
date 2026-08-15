#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <seams_input.hpp>

// The fixture con/tiny_multi_cuh2.con (from the readcon-core test suite)
// holds two frames of two Cu and two H atoms in a 15.3456 x 21.702 x 100
// cell, spec version 2, atom IDs 0 to 3.

TEST_CASE("readCon reads the first frame of a multi-frame con file",
          "[readcon]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readCon("con/tiny_multi_cuh2.con", 1, cloud);

  REQUIRE(cloud.nop == 4);
  REQUIRE(cloud.pts.size() == 4);
  REQUIRE(cloud.idIndexMap.size() == 4);
  REQUIRE(cloud.currentFrame == 1);

  REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(15.3456, 1e-6));
  REQUIRE_THAT(cloud.box[1], Catch::Matchers::WithinAbs(21.702, 1e-6));
  REQUIRE_THAT(cloud.box[2], Catch::Matchers::WithinAbs(100.0, 1e-6));

  // Types carry the atomic number: Cu = 29, H = 1
  REQUIRE(cloud.pts[0].type == 29);
  REQUIRE(cloud.pts[1].type == 29);
  REQUIRE(cloud.pts[2].type == 1);
  REQUIRE(cloud.pts[3].type == 1);

  // First Cu and first H coordinates from the fixture
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(0.6394, 1e-6));
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(6.9753, 1e-6));
  REQUIRE_THAT(cloud.pts[2].x, Catch::Matchers::WithinAbs(8.6823, 1e-6));
  REQUIRE_THAT(cloud.pts[2].z, Catch::Matchers::WithinAbs(11.733, 1e-6));
}

TEST_CASE("readCon selects the requested frame", "[readcon]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readCon("con/tiny_multi_cuh2.con", 2, cloud);

  REQUIRE(cloud.nop == 4);
  REQUIRE(cloud.currentFrame == 2);
  // The H atoms move between the frames; the Cu z relaxes by 1e-4
  REQUIRE_THAT(cloud.pts[2].x, Catch::Matchers::WithinAbs(8.8549, 1e-6));
  REQUIRE_THAT(cloud.pts[2].z, Catch::Matchers::WithinAbs(11.165, 1e-6));
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(6.9752, 1e-6));
}

TEST_CASE("readCon on a missing frame returns an empty cloud", "[readcon]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readCon("con/tiny_multi_cuh2.con", 99, cloud);
  REQUIRE(cloud.nop == 0);
  REQUIRE(cloud.pts.empty());
}
