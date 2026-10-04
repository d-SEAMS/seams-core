#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numbers>

#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
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

  REQUIRE(cloud.box.size() == 3);
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

TEST_CASE("readCon stores a sheared CON cell as the minimage dump box",
          "[readcon]") {
  const double b = 21.702;
  const double xy = b * 0.5;
  const double yy = b * std::numbers::sqrt3 / 2.0;
  const auto path =
      std::filesystem::temp_directory_path() / "seams-sheared-cuh2.con";
  {
    std::ofstream out(path);
    out << "Random Number Seed\n"
        << "{\"con_spec_version\":2}\n"
        << "15.345600\t" << b << "\t100.000000\n"
        << "90.000000\t90.000000\t60.000000\n"
        << "0 0\n"
        << "218 0 1\n"
        << "2\n"
        << "2 2\n"
        << "63.546000 1.007930\n"
        << "Cu\n"
        << "Coordinates of Component 1\n"
        << "0.000000 0.000000 0.000000 1 0\n"
        << std::setprecision(17) << xy << " " << yy << " 0.000000 1 1\n"
        << "H\n"
        << "Coordinates of Component 2\n"
        << "8.6823 9.947 11.733 0 2\n"
        << "7.9421 9.947 11.733 0 3\n";
  }

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readCon(path.string(), 1, cloud);
  REQUIRE(cloud.nop == 4);
  REQUIRE(cloud.box.size() == 6);
  REQUIRE_THAT(cloud.box[3], Catch::Matchers::WithinAbs(xy, 1e-8));
  REQUIRE_THAT(cloud.box[4], Catch::Matchers::WithinAbs(0.0, 1e-8));
  REQUIRE_THAT(cloud.box[5], Catch::Matchers::WithinAbs(0.0, 1e-8));

  double lengths[3];
  nneigh::dumpCellLengths(cloud.box, cloud.boxLow, lengths);
  REQUIRE_THAT(lengths[0], Catch::Matchers::WithinAbs(15.3456, 1e-6));
  REQUIRE_THAT(lengths[1], Catch::Matchers::WithinAbs(yy, 1e-6));
  REQUIRE_THAT(lengths[2], Catch::Matchers::WithinAbs(100.0, 1e-6));
  REQUIRE_THAT(gen::periodicDist(cloud, 0, 1),
               Catch::Matchers::WithinAbs(0.0, 1e-6));

#ifdef SEAMS_HAS_LINKCELL
  const lc_cell cell = nneigh::lammpsBoxToLcCell(cloud.box, cloud.boxLow);
  REQUIRE_THAT(cell.ax, Catch::Matchers::WithinAbs(15.3456, 1e-6));
  REQUIRE_THAT(cell.bx, Catch::Matchers::WithinAbs(xy, 1e-8));
  REQUIRE_THAT(cell.by, Catch::Matchers::WithinAbs(yy, 1e-6));
  REQUIRE_THAT(cell.cz, Catch::Matchers::WithinAbs(100.0, 1e-6));
#endif
  std::filesystem::remove(path);
}

#ifdef SEAMS_HAS_READCON_DB
TEST_CASE("readConCorpus selects a campaign hit and decodes it with readcon",
          "[readcon]") {
  const auto root =
      std::filesystem::temp_directory_path() / "seams-readcon-corpus";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  const int n = sinp::ingestConCorpus(root.string(), 1, "con/tiny_multi_cuh2.con");
  REQUIRE(n == 2);

  molSys::PointCloud<molSys::Point<double>, double> direct;
  molSys::PointCloud<molSys::Point<double>, double> fromDb;
  direct = sinp::readCon("con/tiny_multi_cuh2.con", 1, direct);
  fromDb = sinp::readConCorpus(root.string(), "Cu", "Cu:2|H:2", 1, fromDb);
  REQUIRE(fromDb.nop == direct.nop);
  REQUIRE(fromDb.box == direct.box);
  REQUIRE_THAT(fromDb.pts[0].x, Catch::Matchers::WithinAbs(direct.pts[0].x, 1e-6));
  REQUIRE(fromDb.pts[0].type == 29);
  REQUIRE_THAT(fromDb.pts[2].z, Catch::Matchers::WithinAbs(direct.pts[2].z, 1e-6));

  fromDb = sinp::readConCorpus(root.string(), "Cu", "Cu:2|H:2", 2, fromDb);
  direct = sinp::readCon("con/tiny_multi_cuh2.con", 2, direct);
  REQUIRE_THAT(fromDb.pts[2].x, Catch::Matchers::WithinAbs(direct.pts[2].x, 1e-6));

  fromDb = sinp::readConCorpus(root.string(), "Xe", "", 1, fromDb);
  REQUIRE(fromDb.nop == 0);
  fromDb = sinp::readConCorpus(root.string(), "", "Cu:9", 1, fromDb);
  REQUIRE(fromDb.nop == 0);
  std::filesystem::remove_all(root);
}

TEST_CASE("readConCorpus fans a sharded campaign root out by traj id",
          "[readcon]") {
  const auto root =
      std::filesystem::temp_directory_path() / "seams-readcon-shards";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  {
    std::ofstream manifest(root / "shards.json");
    manifest << "{\"n_shards\":2,\"version\":1}\n";
  }
  const int n = sinp::ingestConCorpus(root.string(), 1, "con/tiny_multi_cuh2.con");
  REQUIRE(n == 2);
  REQUIRE(std::filesystem::is_regular_file(root / "shard_0001" / "data.mdb"));
  REQUIRE_FALSE(std::filesystem::is_regular_file(root / "data.mdb"));

  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readConCorpus(root.string(), "Cu", "", 1, cloud);
  REQUIRE(cloud.nop == 4);
  REQUIRE(cloud.pts[0].type == 29);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(0.6394, 1e-6));
  std::filesystem::remove_all(root);
}
#endif
