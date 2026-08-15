#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <seams_input.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// -- readXYZ tests --

TEST_CASE("readXYZ throws for nonexistent file", "[seams_input]") {
  REQUIRE_THROWS_AS(sinp::readXYZ("/tmp/nonexistent_12345.xyz"),
                    std::runtime_error);
}

TEST_CASE("readXYZ creates a synthetic file and reads it back",
          "[seams_input]") {
  std::string tmpFile = fs::temp_directory_path().append("dseams_test_readxyz.xyz").string();
  {
    std::ofstream f(tmpFile);
    f << "3\n";
    f << "test comment\n";
    f << "O 1.0 2.0 3.0\n";
    f << "O 4.0 5.0 6.0\n";
    f << "O 7.0 8.0 9.0\n";
  }

  auto cloud = sinp::readXYZ(tmpFile);

  REQUIRE(cloud.nop == 3);
  REQUIRE(cloud.pts.size() == 3);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE_THAT(cloud.pts[1].y, Catch::Matchers::WithinAbs(5.0, 1e-10));
  REQUIRE_THAT(cloud.pts[2].z, Catch::Matchers::WithinAbs(9.0, 1e-10));
  REQUIRE(cloud.box.size() == 3);
  REQUIRE(cloud.boxLow.size() == 3);
  REQUIRE_THAT(cloud.boxLow[0], Catch::Matchers::WithinAbs(1.0, 1e-10));
  REQUIRE_THAT(cloud.boxLow[1], Catch::Matchers::WithinAbs(2.0, 1e-10));
  REQUIRE_THAT(cloud.boxLow[2], Catch::Matchers::WithinAbs(3.0, 1e-10));
  REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(6.0, 1e-10));
  REQUIRE_THAT(cloud.box[1], Catch::Matchers::WithinAbs(6.0, 1e-10));
  REQUIRE_THAT(cloud.box[2], Catch::Matchers::WithinAbs(6.0, 1e-10));

  fs::remove(tmpFile);
}

// -- readLammpsTrj tests --

TEST_CASE("readLammpsTrj reads frame 1 from trajectory", "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;

  auto cloud = sinp::readLammpsTrj(
      "traj/exampleTraj.lammpstrj", 1, yCloud);

  REQUIRE(cloud.nop > 0);
  REQUIRE(cloud.pts.size() == static_cast<size_t>(cloud.nop));
  REQUIRE(cloud.box.size() >= 3);
  REQUIRE(cloud.currentFrame == 1);
}

TEST_CASE("readLammpsTrj handles nonexistent file gracefully",
          "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;

  auto cloud = sinp::readLammpsTrj("/tmp/nonexistent.lammpstrj", 1, yCloud);

  // Should return an empty cloud (no crash)
  REQUIRE(cloud.pts.empty());
}

TEST_CASE("readLammpsTrj handles invalid frame number", "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;

  auto cloud = sinp::readLammpsTrj(
      "traj/exampleTraj.lammpstrj", 999, yCloud);

  // Frame 999 doesn't exist; cloud should be empty
  REQUIRE(cloud.pts.empty());
}

TEST_CASE("readLammpsTrj with slice filtering", "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {25.0, 25.0, 25.0};

  auto cloud = sinp::readLammpsTrj(
      "traj/exampleTraj.lammpstrj", 1, yCloud, true, lo, hi);

  REQUIRE(cloud.nop > 0);
  // Some atoms should have inSlice=true, some false
  bool hasInSlice = false;
  bool hasOutSlice = false;
  for (const auto &pt : cloud.pts) {
    if (pt.inSlice) hasInSlice = true;
    else hasOutSlice = true;
  }
  // At least some atoms should be in the slice
  REQUIRE(hasInSlice);
}

// -- readLammpsTrjO tests --

TEST_CASE("readLammpsTrjO reads only oxygen atoms", "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;

  auto cloud = sinp::readLammpsTrjO(
      "traj/exampleTraj.lammpstrj", 1, yCloud, 1);

  REQUIRE(cloud.nop > 0);
  // All atoms should be of typeO
  for (const auto &pt : cloud.pts) {
    REQUIRE(pt.type == 1);
  }
}

TEST_CASE("readLammpsTrjO handles nonexistent file", "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;

  auto cloud =
      sinp::readLammpsTrjO("/tmp/nonexistent.lammpstrj", 1, yCloud, 1);

  REQUIRE(cloud.pts.empty());
}

// -- readLammpsTrjreduced tests --

TEST_CASE("readLammpsTrjreduced reads filtered atoms", "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;

  auto cloud = sinp::readLammpsTrjreduced(
      "traj/exampleTraj.lammpstrj", 1, yCloud, 1);

  REQUIRE(cloud.nop > 0);
  for (const auto &pt : cloud.pts) {
    REQUIRE(pt.type == 1);
  }
}

TEST_CASE("readLammpsTrjreduced with slice and type filter", "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {25.0, 25.0, 25.0};

  auto cloud = sinp::readLammpsTrjreduced(
      "traj/exampleTraj.lammpstrj", 1, yCloud, 1, true, lo, hi);

  // All returned atoms should be of the desired type and in the slice
  for (const auto &pt : cloud.pts) {
    REQUIRE(pt.type == 1);
  }
}

// -- readBonds tests --

TEST_CASE("readBonds returns empty for nonexistent file", "[seams_input]") {
  auto bonds = sinp::readBonds("/tmp/nonexistent_bonds.dat");
  REQUIRE(bonds.empty());
}

TEST_CASE("readBonds reads synthetic bond file", "[seams_input]") {
  std::string tmpFile = fs::temp_directory_path().append("dseams_test_bonds.dat").string();
  {
    std::ofstream f(tmpFile);
    f << "3 Bonds\n";
    f << "1 2 3\n";
    f << "4 5 6\n";
    f << "7 8 9\n";
  }

  auto bonds = sinp::readBonds(tmpFile);

  REQUIRE(bonds.size() == 3);
  REQUIRE(bonds[0].size() == 3);
  REQUIRE(bonds[0][0] == 1);
  REQUIRE(bonds[0][1] == 2);
  REQUIRE(bonds[0][2] == 3);
  REQUIRE(bonds[2][2] == 9);

  fs::remove(tmpFile);
}

// -- atomInSlice tests --

TEST_CASE("atomInSlice returns true for atom inside slice", "[seams_input]") {
  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {5.0, 5.0, 5.0};

  REQUIRE(sinp::atomInSlice(2.0, 3.0, 1.0, lo, hi) == true);
}

TEST_CASE("atomInSlice returns false for atom outside slice", "[seams_input]") {
  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {5.0, 5.0, 5.0};

  REQUIRE(sinp::atomInSlice(6.0, 3.0, 1.0, lo, hi) == false);
}

TEST_CASE("atomInSlice with equal lo/hi ignores that dimension",
          "[seams_input]") {
  // When coordLow[i] == coordHigh[i], that dimension is not filtered
  std::array<double, 3> lo = {0.0, 0.0, 0.0};
  std::array<double, 3> hi = {0.0, 0.0, 0.0};

  // Any coordinate should pass when all dims are "ignored"
  REQUIRE(sinp::atomInSlice(100.0, 200.0, 300.0, lo, hi) == true);
}

TEST_CASE("atomInSlice boundary check", "[seams_input]") {
  std::array<double, 3> lo = {1.0, 1.0, 1.0};
  std::array<double, 3> hi = {5.0, 5.0, 5.0};

  // On the boundary should be included (>=, <=)
  REQUIRE(sinp::atomInSlice(1.0, 1.0, 1.0, lo, hi) == true);
  REQUIRE(sinp::atomInSlice(5.0, 5.0, 5.0, lo, hi) == true);
}


#ifdef SEAMS_HAS_CHEMFILES
TEST_CASE("readChemfiles matches the hand-rolled LAMMPS reader", "[seams_input]") {
  // Same IDs, types and coordinates on both repo dump fixtures, including
  // the type filter. The chemfiles path loses the box origin (its cell has
  // no origin concept), so boxLow is not compared.
  struct Fixture {
    const char *path;
    int filter;
  } cases[] = {{"traj/mW_cubic.lammpstrj", 1},
               {"traj/exampleTraj.lammpstrj", 2}};

  for (const auto &fx : cases) {
    molSys::PointCloud<molSys::Point<double>, double> hand, chem;
    hand = sinp::readLammpsTrjO(fx.path, 1, hand, fx.filter);
    chem = sinp::readChemfiles(fx.path, 1, chem, fx.filter);

    REQUIRE(hand.nop > 0);
    REQUIRE(chem.nop == hand.nop);
    for (int k = 0; k < 3; k++) {
      REQUIRE_THAT(chem.box[k],
                   Catch::Matchers::WithinAbs(hand.box[k], 1e-9));
    }
    for (int i = 0; i < hand.nop; i++) {
      REQUIRE(chem.pts[i].atomID == hand.pts[i].atomID);
      REQUIRE(chem.pts[i].type == hand.pts[i].type);
      REQUIRE_THAT(chem.pts[i].x,
                   Catch::Matchers::WithinAbs(hand.pts[i].x, 1e-9));
      REQUIRE_THAT(chem.pts[i].y,
                   Catch::Matchers::WithinAbs(hand.pts[i].y, 1e-9));
      REQUIRE_THAT(chem.pts[i].z,
                   Catch::Matchers::WithinAbs(hand.pts[i].z, 1e-9));
    }
  }
}

TEST_CASE("readChemfiles reports unreadable files without terminating",
          "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readChemfiles("no/such/file.lammpstrj", 1, cloud, -1);
  REQUIRE(cloud.nop == 0);
}
#endif // SEAMS_HAS_CHEMFILES
