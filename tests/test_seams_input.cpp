#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/resource.h>
#include <unistd.h>

namespace fs = std::filesystem;

// -- readXYZ tests --

TEST_CASE("readXYZ throws for nonexistent file", "[seams_input]") {
  REQUIRE_THROWS_AS(sinp::readXYZ("missing/nonexistent_12345.xyz"),
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

TEST_CASE("readXYZ gives each element its own type", "[seams_input]") {
  std::string tmpFile = fs::temp_directory_path().append("dseams_test_readxyz_sin.xyz").string();
  {
    std::ofstream f(tmpFile);
    f << "3\n\n";
    f << "Si 0.0 0.0 0.0\n";
    f << "N 1.7 0.0 0.0\n";
    f << "N 2.8 0.0 0.0\n";
  }
  auto cloud = sinp::readXYZ(tmpFile);
  REQUIRE(cloud.nop == 3);
  REQUIRE(cloud.pts[0].type == 1);
  REQUIRE(cloud.pts[1].type == 2);
  REQUIRE(cloud.pts[2].type == 2);
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

  auto cloud = sinp::readLammpsTrj("missing/nonexistent.lammpstrj", 1, yCloud);

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
      sinp::readLammpsTrjO("missing/nonexistent.lammpstrj", 1, yCloud, 1);

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

// -- LAMMPS dump live session (ReaderNative position + offset table) --

namespace {
std::string writeTinyDump() {
  auto path =
      fs::temp_directory_path() / "dseams_test_lammps_index.lammpstrj";
  std::ofstream f(path);
  for (int frame = 0; frame < 3; ++frame) {
    f << "ITEM: TIMESTEP\n" << (100 * frame) << "\n";
    f << "ITEM: NUMBER OF ATOMS\n2\n";
    f << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n";
    f << "ITEM: ATOMS id type x y z\n";
    f << "1 1 " << frame + 0.25 << " 0 0\n";
    f << "2 1 " << frame + 0.75 << " 0 0\n";
  }
  return path.string();
}

std::string writeUnwrappedDump() {
  auto path =
      fs::temp_directory_path() / "dseams_test_lammps_xu.lammpstrj";
  std::ofstream f(path);
  f << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n";
  f << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n";
  f << "ITEM: ATOMS id type xu yu zu\n";
  f << "7 2 11.5 12.25 -3.0\n";
  return path.string();
}

std::string writeEmptyMiddleDump() {
  auto path =
      fs::temp_directory_path() / "dseams_test_lammps_empty_frame.lammpstrj";
  std::ofstream f(path);
  const int counts[3] = {2, 0, 2};
  for (int frame = 0; frame < 3; ++frame) {
    f << "ITEM: TIMESTEP\n" << frame << "\n";
    f << "ITEM: NUMBER OF ATOMS\n" << counts[frame] << "\n";
    f << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n";
    f << "ITEM: ATOMS id type x y z\n";
    for (int i = 0; i < counts[frame]; ++i) {
      f << (i + 1) << " 1 " << (i + 0.25) << " 0 0\n";
    }
  }
  return path.string();
}
} // namespace

TEST_CASE("nLammpsFrames counts ITEM: TIMESTEP markers", "[seams_input]") {
  REQUIRE(sinp::nLammpsFrames("missing/nonexistent.lammpstrj") == 0);
  REQUIRE(sinp::nLammpsFrames("traj/mW_cubic.lammpstrj") == 11);
  const auto tiny = writeTinyDump();
  sinp::dropLammpsDumpIndex(tiny);
  REQUIRE(sinp::nLammpsFrames(tiny) == 3);
  fs::remove(tiny);
  sinp::dropLammpsDumpIndex(tiny);
}

TEST_CASE("readLammpsTrjreduced seeks the last frame by index",
          "[seams_input]") {
  const auto tiny = writeTinyDump();
  sinp::dropLammpsDumpIndex(tiny);
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  auto last = sinp::readLammpsTrjreduced(tiny, 3, yCloud, 1);
  REQUIRE(last.nop == 2);
  REQUIRE(last.currentFrame == 3);
  REQUIRE_THAT(last.pts[0].x, Catch::Matchers::WithinAbs(2.25, 1e-10));
  REQUIRE_THAT(last.pts[1].x, Catch::Matchers::WithinAbs(2.75, 1e-10));
  auto first = sinp::readLammpsTrjreduced(tiny, 1, yCloud, 1);
  REQUIRE_THAT(first.pts[0].x, Catch::Matchers::WithinAbs(0.25, 1e-10));
  auto mid = sinp::readLammpsTrjreduced(tiny, 2, yCloud, 1);
  REQUIRE_THAT(mid.pts[0].x, Catch::Matchers::WithinAbs(1.25, 1e-10));
  fs::remove(tiny);
  sinp::dropLammpsDumpIndex(tiny);
}

TEST_CASE("sequential load_frame walks keep the live dump position",
          "[seams_input]") {
  const auto tiny = writeTinyDump();
  sinp::dropLammpsDumpIndex(tiny);
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  for (int frame = 1; frame <= 3; ++frame) {
    auto cloud = sinp::readLammpsTrj(tiny, frame, yCloud);
    REQUIRE(cloud.nop == 2);
    REQUIRE(cloud.currentFrame == frame);
    REQUIRE_THAT(cloud.pts[0].x,
                 Catch::Matchers::WithinAbs(frame - 0.75, 1e-10));
  }
  auto again = sinp::readLammpsTrj(tiny, 1, yCloud);
  REQUIRE_THAT(again.pts[0].x, Catch::Matchers::WithinAbs(0.25, 1e-10));
  fs::remove(tiny);
  sinp::dropLammpsDumpIndex(tiny);
}

TEST_CASE("indexed last-frame read of mW_cubic matches a full scan",
          "[seams_input]") {
  REQUIRE(sinp::nLammpsFrames("traj/mW_cubic.lammpstrj") == 11);
  molSys::PointCloud<molSys::Point<double>, double> a, b;
  auto last = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 11, a, 1);
  auto first = sinp::readLammpsTrjO("traj/mW_cubic.lammpstrj", 1, b, 1);
  REQUIRE(last.nop == first.nop);
  REQUIRE(last.nop == 4096);
  REQUIRE(last.currentFrame == 11);
  REQUIRE(first.currentFrame == 1);
}

TEST_CASE("forEachLammpsFrame matches sequential reads", "[seams_input]") {
  const auto tiny = writeTinyDump();
  sinp::dropLammpsDumpIndex(tiny);
  std::vector<double> xs(3, -1.0);
  std::vector<int> nops(3, 0);
  sinp::forEachLammpsFrame(
      tiny, 1, 3, 1,
      [&](int frame, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
        const auto i = static_cast<std::size_t>(frame - 1);
        nops[i] = cloud.nop;
        if (!cloud.pts.empty()) {
          xs[i] = cloud.pts[0].x;
        }
      },
      2);
  REQUIRE(nops[0] == 2);
  REQUIRE(nops[1] == 2);
  REQUIRE(nops[2] == 2);
  REQUIRE_THAT(xs[0], Catch::Matchers::WithinAbs(0.25, 1e-10));
  REQUIRE_THAT(xs[1], Catch::Matchers::WithinAbs(1.25, 1e-10));
  REQUIRE_THAT(xs[2], Catch::Matchers::WithinAbs(2.25, 1e-10));
  fs::remove(tiny);
  sinp::dropLammpsDumpIndex(tiny);
}

TEST_CASE("forEachLammpsFrame deals frames to parts in serpentine order",
          "[seams_input]") {
  // Frames 2 to 11 in rounds of three, odd rounds reversed: {2, 3, 4},
  // {7, 6, 5}, {8, 9, 10}, {-, -, 11}
  const std::vector<std::vector<int>> expected = {
      {2, 7, 8}, {3, 6, 9}, {4, 5, 10, 11}, {}};
  for (int part = 0; part < 4; part++) {
    std::vector<int> seen(12, 0);
    sinp::forEachLammpsFrame(
        "traj/mW_cubic.lammpstrj", 2, 11, 1,
        [&](int frame, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
          seen[static_cast<std::size_t>(frame)] =
              cloud.currentFrame == frame ? 1 : -1;
        },
        2, part, 3);
    std::vector<int> got;
    for (int frame = 0; frame < 12; frame++) {
      if (seen[static_cast<std::size_t>(frame)] != 0) {
        REQUIRE(seen[static_cast<std::size_t>(frame)] == 1);
        got.push_back(frame);
      }
    }
    REQUIRE(got == expected[static_cast<std::size_t>(part)]);
  }
}

TEST_CASE("empty NUMBER OF ATOMS frames are valid snapshots",
          "[seams_input]") {
  const auto dump = writeEmptyMiddleDump();
  sinp::dropLammpsDumpIndex(dump);
  REQUIRE(sinp::nLammpsFrames(dump) == 3);

  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  auto first = sinp::readLammpsTrjO(dump, 1, yCloud, 1);
  REQUIRE(first.nop == 2);
  REQUIRE(first.box.size() >= 3);

  auto empty = sinp::readLammpsTrjO(dump, 2, yCloud, 1);
  REQUIRE(empty.nop == 0);
  REQUIRE(empty.pts.empty());
  REQUIRE(empty.currentFrame == 2);
  REQUIRE(empty.box.size() >= 3);

  auto last = sinp::readLammpsTrjO(dump, 3, yCloud, 1);
  REQUIRE(last.nop == 2);
  REQUIRE(last.currentFrame == 3);

  std::vector<int> nops(3, -1);
  sinp::forEachLammpsFrame(
      dump, 1, 3, 1,
      [&](int frame, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
        nops[static_cast<std::size_t>(frame - 1)] = cloud.nop;
      },
      1);
  REQUIRE(nops[0] == 2);
  REQUIRE(nops[1] == 0);
  REQUIRE(nops[2] == 2);

  fs::remove(dump);
  sinp::dropLammpsDumpIndex(dump);
}

TEST_CASE("readLammpsTrj binds xu yu zu when x y z are absent",
          "[seams_input]") {
  const auto path = writeUnwrappedDump();
  sinp::dropLammpsDumpIndex(path);
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  auto cloud = sinp::readLammpsTrj(path, 1, yCloud);
  REQUIRE(cloud.nop == 1);
  REQUIRE(cloud.pts[0].atomID == 7);
  REQUIRE(cloud.pts[0].type == 2);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(11.5, 1e-10));
  REQUIRE_THAT(cloud.pts[0].y, Catch::Matchers::WithinAbs(12.25, 1e-10));
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(-3.0, 1e-10));
  fs::remove(path);
  sinp::dropLammpsDumpIndex(path);
}

TEST_CASE("readLammpsTrjO parses scientific-notation dump coordinates",
          "[seams_input]") {
  auto path = fs::temp_directory_path() / "dseams_test_lammps_sci.lammpstrj";
  {
    std::ofstream f(path);
    f << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n";
    f << "ITEM: BOX BOUNDS pp pp pp\n0 1.0e1\n0 1e1\n0 10\n";
    f << "ITEM: ATOMS id type x y z\n";
    f << "1 1 1.25e0 2.5E-1 -3.0e0\n";
    f << "2 1 4.2e1 0 1.5e-1\n";
  }
  sinp::dropLammpsDumpIndex(path.string());
  molSys::PointCloud<molSys::Point<double>, double> yCloud;
  auto cloud = sinp::readLammpsTrjO(path.string(), 1, yCloud, 1);
  REQUIRE(cloud.nop == 2);
  REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(10.0, 1e-10));
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(1.25, 1e-10));
  REQUIRE_THAT(cloud.pts[0].y, Catch::Matchers::WithinAbs(0.25, 1e-10));
  REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(-3.0, 1e-10));
  REQUIRE_THAT(cloud.pts[1].x, Catch::Matchers::WithinAbs(42.0, 1e-10));
  REQUIRE_THAT(cloud.pts[1].z, Catch::Matchers::WithinAbs(0.15, 1e-10));
  fs::remove(path);
  sinp::dropLammpsDumpIndex(path.string());
}

// -- readBonds tests --

TEST_CASE("readBonds returns empty for nonexistent file", "[seams_input]") {
  auto bonds = sinp::readBonds("missing/nonexistent_bonds.dat");
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

TEST_CASE("readLammpsTrj assigns ids when the dump has no id column",
          "[seams_input]") {
  const auto path =
      fs::temp_directory_path() / "dseams_test_noid.lammpstrj";
  {
    std::ofstream f(path);
    f << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\n"
      << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      << "ITEM: ATOMS type x y z\n"
      << "1 1 1 1\n1 4 1 1\n1 1 5 1\n";
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path.string(), 1, cloud);
  REQUIRE(cloud.nop == 3);
  REQUIRE(cloud.idIndexMap.size() == 3);
  REQUIRE(cloud.pts[0].atomID == 1);
  REQUIRE(cloud.pts[1].atomID == 2);
  REQUIRE(cloud.pts[2].atomID == 3);
  fs::remove(path);
  sinp::dropLammpsDumpIndex(path.string());
}

TEST_CASE("readLammpsTrj keeps a mol column that sits in column 0",
          "[seams_input]") {
  const auto path =
      fs::temp_directory_path() / "dseams_test_molfirst.lammpstrj";
  {
    std::ofstream f(path);
    f << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
      << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      << "ITEM: ATOMS mol id type x y z\n"
      << "10 1 1 1 1 1\n20 2 1 3 1 1\n";
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path.string(), 1, cloud);
  REQUIRE(cloud.pts[0].atomID == 1);
  REQUIRE(cloud.pts[0].molID == 10);
  REQUIRE(cloud.pts[1].atomID == 2);
  REQUIRE(cloud.pts[1].molID == 20);
  fs::remove(path);
  sinp::dropLammpsDumpIndex(path.string());
}

TEST_CASE("readLammpsTrj sets nop from the atoms when the count overflows",
          "[seams_input]") {
  const auto path =
      fs::temp_directory_path() / "dseams_test_huge.lammpstrj";
  {
    std::ofstream f(path);
    f << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n999999999999\n"
      << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      << "ITEM: ATOMS id type x y z\n1 1 1 1 1\n2 1 2 1 1\n";
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path.string(), 1, cloud);
  REQUIRE(cloud.nop == 2);
  REQUIRE(cloud.pts.size() == 2);
  fs::remove(path);
  sinp::dropLammpsDumpIndex(path.string());
}

TEST_CASE("readXYZ sets nop to the atoms that were parsed", "[seams_input]") {
  const auto path = fs::temp_directory_path() / "dseams_test_short.xyz";
  {
    std::ofstream f(path);
    f << "5\ncomment\nO 0 0 0\nO 1 0 0\n";
  }
  auto cloud = sinp::readXYZ(path.string());
  REQUIRE(cloud.nop == 2);
  REQUIRE(cloud.pts.size() == 2);
  fs::remove(path);
}

TEST_CASE("readXYZ skips a coordinate that is not a number", "[seams_input]") {
  const auto path = fs::temp_directory_path() / "dseams_test_bad.xyz";
  {
    std::ofstream f(path);
    f << "2\ncomment\nO 0 0 0\nO foo bar baz\n";
  }
  auto cloud = sinp::readXYZ(path.string());
  REQUIRE(cloud.nop == 1);
  REQUIRE(cloud.pts.size() == 1);
  fs::remove(path);
}

TEST_CASE("readBonds drops a constant trailing 0/1 flag", "[seams_input]") {
  const auto path = fs::temp_directory_path() / "dseams_test_ringsflag.dat";
  {
    std::ofstream f(path);
    f << "2 Bonds\n272 214 906 1361 388 1\n388 1361 1042 1548 237 1\n";
  }
  auto bonds = sinp::readBonds(path.string());
  REQUIRE(bonds.size() == 2);
  REQUIRE(bonds[0].size() == 5);
  REQUIRE(bonds[0].back() == 388);
  fs::remove(path);
}

TEST_CASE("lengths and angles of 90 degrees stay an orthorhombic dump box",
          "[seams_input]") {
  std::vector<double> box;
  std::vector<double> boxLow;
  REQUIRE(nneigh::setBoxFromLengthsAngles(box, boxLow, 10, 12, 14, 90, 90, 90));
  REQUIRE(box.size() == 3);
  REQUIRE_THAT(box[0], Catch::Matchers::WithinAbs(10, 1e-9));
  REQUIRE_THAT(box[1], Catch::Matchers::WithinAbs(12, 1e-9));
  REQUIRE_THAT(box[2], Catch::Matchers::WithinAbs(14, 1e-9));
  const double vol = nneigh::dumpVolume(box, boxLow);
  REQUIRE_THAT(vol, Catch::Matchers::WithinAbs(10 * 12 * 14, 1e-6));
}

TEST_CASE("a 60 degree gamma uses the cell volume not the bound-span product",
          "[seams_input]") {
  std::vector<double> box;
  std::vector<double> boxLow;
  REQUIRE(nneigh::setBoxFromLengthsAngles(box, boxLow, 10, 10, 10, 90, 90, 60));
  REQUIRE(box.size() == 6);
  const double vol = nneigh::dumpVolume(box, boxLow);
  const double expect = 1000.0 * std::sin(60.0 * 3.14159265358979323846 / 180.0);
  REQUIRE_THAT(vol, Catch::Matchers::WithinAbs(expect, 1e-6));
  REQUIRE(box[0] * box[1] * box[2] > vol + 1.0);
}

TEST_CASE("readChemfiles reports unreadable files without terminating",
          "[seams_input]") {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readChemfiles("no/such/file.lammpstrj", 1, cloud, -1);
  REQUIRE(cloud.nop == 0);
}
#endif // SEAMS_HAS_CHEMFILES

namespace {

void writeDumpLines(const fs::path &path, int declared, int lines) {
  std::ofstream out(path);
  out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n" << declared
      << "\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      << "ITEM: ATOMS id type x y z\n";
  for (int i = 0; i < lines; ++i) {
    out << i << " 1 0 0 0\n";
  }
}

} // namespace

// The session from the first read keeps its descriptor. Drop it before
// the limit, so the next read has to open and that open fails.
TEST_CASE("a second open fails when the process is out of descriptors",
          "[seams_input]") {
  const auto dir = fs::temp_directory_path() / "dseams_nofile_descriptor";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "nofile.lammpstrj";
  writeDumpLines(dump, 220000, 220000);
  const auto path = dump.string();
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path, 1, cloud);
  REQUIRE(cloud.nop == 220000);
  sinp::dropLammpsDumpIndex(path);
  struct rlimit saved {};
  REQUIRE(::getrlimit(RLIMIT_NOFILE, &saved) == 0);
  std::vector<int> held;
  int maxfd = 0;
  while (held.size() < 128) {
    const int fd = ::open("/dev/null", O_RDONLY);
    if (fd < 0) {
      break;
    }
    if (fd > maxfd) {
      maxfd = fd;
    }
    held.push_back(fd);
  }
  {
    struct Restore {
      struct rlimit saved;
      std::vector<int> *held;
      bool armed = false;
      ~Restore() {
        if (armed) {
          ::setrlimit(RLIMIT_NOFILE, &saved);
        }
        if (held != nullptr) {
          for (const int fd : *held) {
            ::close(fd);
          }
        }
      }
    } restore{saved, &held, false};
    REQUIRE(!held.empty());
    struct rlimit tight = saved;
    tight.rlim_cur = static_cast<rlim_t>(maxfd + 1);
    REQUIRE(::setrlimit(RLIMIT_NOFILE, &tight) == 0);
    restore.armed = true;
    cloud = sinp::readLammpsTrj(path, 1, cloud);
    REQUIRE(cloud.nop == 0);
    sinp::dropLammpsDumpIndex(path);
  }
  fs::remove_all(dir);
}
