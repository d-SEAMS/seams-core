// Copyright (c) d-SEAMS developers
#include "fresh_temp.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <thread>
#include <vector>

#include <csignal>
#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

namespace fs = std::filesystem;

// -- readXYZ tests --

TEST_CASE("readXYZ throws for nonexistent file", "[seams_input]") {
  REQUIRE_THROWS_AS(sinp::readXYZ("missing/nonexistent_12345.xyz"),
                    std::runtime_error);
}

TEST_CASE("readXYZ creates a synthetic file and reads it back",
          "[seams_input]") {
  std::string tmpFile = dseams_test::freshFile("dseams_test_readxyz.xyz").string();
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
  std::string tmpFile = dseams_test::freshFile("dseams_test_readxyz_sin.xyz").string();
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
      dseams_test::freshFile("dseams_test_lammps_index.lammpstrj");
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
      dseams_test::freshFile("dseams_test_lammps_xu.lammpstrj");
  std::ofstream f(path);
  f << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n";
  f << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n";
  f << "ITEM: ATOMS id type xu yu zu\n";
  f << "7 2 11.5 12.25 -3.0\n";
  return path.string();
}

std::string writeEmptyMiddleDump() {
  auto path =
      dseams_test::freshFile("dseams_test_lammps_empty_frame.lammpstrj");
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
  auto path = dseams_test::freshFile("dseams_test_lammps_sci.lammpstrj");
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
  std::string tmpFile = dseams_test::freshFile("dseams_test_bonds.dat").string();
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
  REQUIRE(sinp::atomInSlice(1.0, 6.0, 1.0, lo, hi) == false);
  REQUIRE(sinp::atomInSlice(1.0, 1.0, 6.0, lo, hi) == false);
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
      dseams_test::freshFile("dseams_test_noid.lammpstrj");
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
      dseams_test::freshFile("dseams_test_molfirst.lammpstrj");
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
      dseams_test::freshFile("dseams_test_huge.lammpstrj");
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

TEST_CASE("readXYZ rejects a short frame", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_test_short.xyz");
  {
    std::ofstream f(path);
    f << "5\ncomment\nO 0 0 0\nO 1 0 0\n";
  }
  try {
    sinp::readXYZ(path.string());
    FAIL("a short XYZ must not come back as a frame");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()).find("truncated") != std::string::npos);
  }
  fs::remove(path);
}

TEST_CASE("readXYZ rejects a coordinate that is not a number", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_test_bad.xyz");
  {
    std::ofstream f(path);
    f << "2\ncomment\nO 0 0 0\nO foo bar baz\n";
  }
  try {
    sinp::readXYZ(path.string());
    FAIL("a bad coordinate must not come back as a frame");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()).find("not a number") != std::string::npos);
  }
  fs::remove(path);
}

TEST_CASE("readBonds drops a constant trailing 0/1 flag", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_test_ringsflag.dat");
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

std::string oneDump(const std::string &atoms, int declared, const std::string &tail = "") {
  return "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n" + std::to_string(declared) +
         "\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\nITEM: ATOMS id type x y z\n" +
         atoms + tail;
}

void writeText(const fs::path &path, const std::string &text) {
  std::ofstream out(path);
  out << text;
}

void writeDumpLines(const fs::path &path, int declared, int lines) {
  std::ofstream out(path);
  out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n" << declared
      << "\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      << "ITEM: ATOMS id type x y z\n";
  for (int i = 0; i < lines; ++i) {
    out << i << " 1 0 0 0\n";
  }
}

std::vector<int> fdsForPath(const std::string &path) {
  std::vector<int> ids;
  for (const auto &ent : fs::directory_iterator("/proc/self/fd")) {
    std::error_code ec;
    const auto link = fs::read_symlink(ent.path(), ec);
    if (ec) {
      continue;
    }
    if (link.string().rfind(path, 0) == 0) {
      ids.push_back(std::stoi(ent.path().filename().string()));
    }
  }
  std::sort(ids.begin(), ids.end());
  return ids;
}

std::string conText(const char *cell, const char *angles, const char *ids) {
  return std::string("Random Number Seed\n{\"con_spec_version\":2}\n") + cell +
         "\n" + angles +
         "\n0 0\n218 0 1\n1\n2\n63.546000\nCu\nCoordinates of Component 1\n" +
         ids;
}

} // namespace

TEST_CASE("dump and bond readers take the rejected inputs", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  molSys::PointCloud<molSys::Point<double>, double> cloud;

  const auto dup = dir / "dup.lammpstrj";
  writeText(dup, oneDump("1 1 0 0 0\n1 1 1 0 0\n", 2));
  cloud = sinp::readLammpsTrj(dup.string(), 1, cloud);
  REQUIRE(cloud.nop == 2);

  const auto junk = dir / "junk.lammpstrj";
  writeText(junk, oneDump("1 1 abc 0 0\n2 1 1 0 0\n", 2));
  cloud = sinp::readLammpsTrj(junk.string(), 1, cloud);

  const auto scaled = dir / "scaled.lammpstrj";
  writeText(scaled,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type xu ys zs\n1 1 0.1 0.2 0.3\n");
  cloud = sinp::readLammpsTrj(scaled.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);

  int seen = 0;
  sinp::forEachLammpsFrame(dup.string(), 0, 0, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; });
  REQUIRE(seen == 1);
  sinp::forEachLammpsFrame(dup.string(), 4, 1, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; });
  REQUIRE(seen == 1);

  const auto empty = dir / "empty.lammpstrj";
  writeText(empty, "");
  REQUIRE(sinp::nLammpsFrames(empty.string()) == 0);
  const auto prose = dir / "prose.lammpstrj";
  writeText(prose, "hello\n");
  REQUIRE(sinp::nLammpsFrames(prose.string()) == 0);
  cloud = sinp::readLammpsTrj(prose.string(), 0, cloud);
  cloud = sinp::readLammpsTrj(prose.string(), 1, cloud);

  const auto tail = dir / "tail.lammpstrj";
  writeText(tail, oneDump("1 1 0 0 0\n", 1, "not a timestep\n"));
  cloud = sinp::readLammpsTrj(tail.string(), 1, cloud);
  cloud = sinp::readLammpsTrj(tail.string(), 2, cloud);

  const auto grow = dir / "grow.lammpstrj";
  writeText(grow, oneDump("1 1 0 0 0\n", 1));
  cloud = sinp::readLammpsTrj(grow.string(), 1, cloud);
  {
    std::ofstream extra(grow, std::ios::app);
    extra << "\n";
  }
  cloud = sinp::readLammpsTrj(grow.string(), 1, cloud);

  const auto flipped = dir / "flip.lammpstrj";
  writeText(flipped, oneDump("1 1 0 0 0\n", 1));
  cloud = sinp::readLammpsTrj(flipped.string(), 1, cloud);
  const auto stamp = fs::last_write_time(flipped);
  {
    std::fstream edit(flipped, std::ios::in | std::ios::out | std::ios::binary);
    edit.seekp(0);
    edit.write("JTEM: TIMESTEP", 14);
  }
  fs::last_write_time(flipped, stamp);
  cloud = sinp::readLammpsTrj(flipped.string(), 1, cloud);

  const auto asDir = dir / "subdir";
  fs::create_directory(asDir);
  cloud = sinp::readLammpsTrj(asDir.string(), 1, cloud);
  REQUIRE(cloud.nop == 0);

  const auto bonds = dir / "bonds.dat";
  writeText(bonds, "3 Bonds\n1 2 abc\n1 2 3 0\n4 5 6 1\n");
  const auto parsed = sinp::readBonds(bonds.string());
  REQUIRE(parsed.size() == 2);

  const auto xyz = dir / "odd.xyz";
  writeText(xyz, "nope\ncomment\n");
  REQUIRE_THROWS_AS(sinp::readXYZ(xyz.string()), std::runtime_error);
  writeText(xyz, "-1\ncomment\n");
  REQUIRE_THROWS_AS(sinp::readXYZ(xyz.string()), std::runtime_error);
  writeText(xyz, "3\ncomment\n\n1 0 0 0\nO 1 0 0\n");
  REQUIRE_THROWS_AS(sinp::readXYZ(xyz.string()), std::runtime_error);
  const auto tri = dir / "tri.lammpstrj";
  writeText(tri,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS xy xz yz pp pp pp\n0 10 0.5\n0 10 0.2\n0 10 0.1\n"
            "ITEM: ATOMS id type x y z\n1 1 1 2 3\n");
  cloud = sinp::readLammpsTrj(tri.string(), 1, cloud);
  REQUIRE(cloud.box.size() >= 6);

  sinp::dropLammpsDumpIndex(dup.string());
  fs::remove_all(dir);
}

#ifdef SEAMS_HAS_READCON
TEST_CASE("con cells and the corpus errors", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const char *atoms = "   1.0    2.0    3.0 1    0\n   4.0    5.0    6.0 1    1\n";
  const auto good = dir / "good.con";
  writeText(good, conText("10.0\t10.0\t10.0", "90.0\t90.0\t90.0", atoms));
  cloud = sinp::readCon(good.string(), 0, cloud);
  REQUIRE(cloud.nop == 0);
  cloud = sinp::readCon(good.string(), 1, cloud);
  REQUIRE(cloud.nop == 2);
  cloud = sinp::readCon(good.string(), 9, cloud);
  REQUIRE(cloud.nop == 0);
  cloud = sinp::readCon((dir / "missing.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 0);
  writeText(dir / "zero.con", conText("0.0\t0.0\t0.0", "90.0\t90.0\t90.0", atoms));
  cloud = sinp::readCon((dir / "zero.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 2);
  writeText(dir / "bzero.con", conText("10.0\t0.0\t10.0", "90.0\t90.0\t90.0", atoms));
  cloud = sinp::readCon((dir / "bzero.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 2);
  writeText(dir / "czero.con", conText("10.0\t10.0\t0.0", "90.0\t90.0\t90.0", atoms));
  cloud = sinp::readCon((dir / "czero.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 2);
  writeText(dir / "empty.con",
            "Random Number Seed\n{\"con_spec_version\":2}\n10.0\t10.0\t10.0\n"
            "90.0\t90.0\t90.0\n0 0\n218 0 1\n1\n0\n63.546000\nCu\n"
            "Coordinates of Component 1\n");
  cloud = sinp::readCon((dir / "empty.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 0);
  writeText(dir / "one.con",
            "Random Number Seed\n{\"con_spec_version\":2}\n10.0\t10.0\t10.0\n"
            "90.0\t90.0\t90.0\n0 0\n218 0 1\n1\n1\n63.546000\nCu\n"
            "Coordinates of Component 1\n   1.0    2.0    3.0 1    0\n");
  cloud = sinp::readCon((dir / "one.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  REQUIRE(cloud.pts[0].x == 1.0);
  writeText(dir / "flat.con", conText("10.0\t10.0\t10.0", "90.0\t90.0\t0.0", atoms));
  cloud = sinp::readCon((dir / "flat.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 2);
  writeText(dir / "shear.con", conText("10.0\t10.0\t10.0", "90.0\t0.0\t90.0", atoms));
  cloud = sinp::readCon((dir / "shear.con").string(), 1, cloud);
  REQUIRE(cloud.nop == 2);
  writeText(dir / "tilt.con", conText("10.0\t10.0\t10.0", "80.0\t75.0\t70.0", atoms));
  cloud = sinp::readCon((dir / "tilt.con").string(), 1, cloud);
  REQUIRE(cloud.box.size() >= 6);

#ifdef SEAMS_HAS_READCON_DB
  const auto corpus = dir / "corpus";
  fs::create_directories(corpus);
  REQUIRE(sinp::ingestConCorpus(corpus.string(), 1, good.string()) > 0);
  REQUIRE(sinp::ingestConCorpus(corpus.string(), 2, good.string()) > 0);
  cloud = sinp::readConCorpus(corpus.string(), "", "", 1, cloud);
  REQUIRE(cloud.nop == 2);
  cloud = sinp::readConCorpus(corpus.string(), "Cu", "", 2, cloud);
  cloud = sinp::readConCorpus(corpus.string(), "", "", 9, cloud);
  REQUIRE(cloud.nop == 0);
  REQUIRE(sinp::ingestConCorpus(corpus.string(), 3, (dir / "missing.con").string()) < 0);
  const auto blocker = dir / "not-a-dir";
  writeText(blocker, "x");
  REQUIRE(sinp::ingestConCorpus(blocker.string(), 1, good.string()) < 0);
  cloud = sinp::readConCorpus(blocker.string(), "", "", 1, cloud);
  cloud = sinp::readConCorpus((dir / "no-such-corpus").string(), "", "", 1, cloud);
  cloud = sinp::readConCorpus(corpus.string(), "", "", 0, cloud);
  cloud = sinp::readConCorpus(corpus.string(), "Cu", "", 1, cloud);
  writeText(corpus / "shards.json", "{}\n");
  cloud = sinp::readConCorpus(corpus.string(), "", "", 1, cloud);
  writeText(corpus / "shards.json", "{\"n_shards\"}\n");
  cloud = sinp::readConCorpus(corpus.string(), "", "", 1, cloud);
  writeText(corpus / "shards.json", "{\"n_shards\": 0}\n");
  cloud = sinp::readConCorpus(corpus.string(), "", "", 1, cloud);
  writeText(corpus / "shards.json", "{\"n_shards\": 9000}\n");
  cloud = sinp::readConCorpus(corpus.string(), "", "", 1, cloud);
  writeText(corpus / "shards.json", "{\"n_shards\": 1}\n");
  fs::create_directories(corpus / "shard_0000");
  writeText(corpus / "shard_0000" / "data.mdb", "");
  cloud = sinp::readConCorpus(corpus.string(), "", "H2", 4, cloud);
#endif
  fs::remove_all(dir);
}
#endif

TEST_CASE("a dump without a timestep and a raced missing frame",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  const auto blank = dir / "blank.lammpstrj";
  writeText(blank, "this file has no timestep header\n");
  cloud = sinp::readLammpsTrj(blank.string(), 1, cloud);
  REQUIRE(cloud.nop == 0);
  cloud = sinp::readLammpsTrj((dir / "as-dir").string(), 1, cloud);
  fs::create_directories(dir / "as-dir");
  cloud = sinp::readLammpsTrj((dir / "as-dir").string(), 1, cloud);
  cloud = sinp::readLammpsTrj("/proc/self/mem", 1, cloud);
  const auto empty = dir / "empty.lammpstrj";
  writeText(empty, "");
  cloud = sinp::readLammpsTrj(empty.string(), 1, cloud);

  const auto dump = dir / "race.lammpstrj";
  std::string body =
      "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3000\n"
      "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      "ITEM: ATOMS id type x y z\n";
  body.reserve(body.size() + 3000 * 16);
  for (int i = 0; i < 3000; ++i) {
    body += std::to_string(i);
    body += " 1 0 0 0\n";
  }
  writeText(dump, body);
  std::atomic<bool> go{true};
  std::thread worker([&] {
    // A raced dump may throw. The exception has to stay in this thread.
    try {
      molSys::PointCloud<molSys::Point<double>, double> local;
      while (go.load()) {
        local = sinp::readLammpsTrj(dump.string(), 1, local);
      }
    } catch (const std::exception &) {
    }
  });
  try {
    for (int i = 0; i < 200; ++i) {
      cloud = sinp::readLammpsTrj(dump.string(), 0, cloud);
      cloud = sinp::readLammpsTrj(dump.string(), 9, cloud);
    }
    fs::remove(dump);
    for (int i = 0; i < 50; ++i) {
      cloud = sinp::readLammpsTrj(dump.string(), 1, cloud);
    }
  } catch (const std::exception &) {
  }
  go.store(false);
  worker.join();
  REQUIRE(cloud.nop == 0);
  const auto flicker = dir / "flicker.lammpstrj";
  writeText(flicker, body);
  std::atomic<bool> spin{true};
  std::thread remover([&] {
    try {
      while (spin.load()) {
        fs::remove(flicker);
        writeText(flicker, body);
      }
    } catch (const std::exception &) {
    }
  });
  try {
    for (int i = 0; i < 2000; ++i) {
      cloud = sinp::readLammpsTrj(flicker.string(), 1, cloud);
    }
  } catch (const std::exception &) {
  }
  spin.store(false);
  remover.join();
  std::atomic<bool> spinFrames{true};
  std::thread frameRemover([&] {
    try {
      while (spinFrames.load()) {
        fs::remove(flicker);
        writeText(flicker, body);
      }
    } catch (const std::exception &) {
    }
  });
  for (int i = 0; i < 2000; ++i) {
    sinp::dropLammpsDumpIndex(flicker.string());
    REQUIRE(sinp::nLammpsFrames(flicker.string()) >= 0);
  }
  spinFrames.store(false);
  frameRemover.join();
  REQUIRE(sinp::nLammpsFrames("/proc/kcore") == 0);
  REQUIRE(sinp::nLammpsFrames("/proc/self/pagemap") == 0);
  std::atomic<bool> hold{true};
  std::thread holder([&] {
    try {
      molSys::PointCloud<molSys::Point<double>, double> local;
      while (hold.load()) {
        local = sinp::readLammpsTrj(dump.string(), 1, local);
      }
    } catch (const std::exception &) {
    }
  });
  writeText(dump, body);
  try {
    for (int i = 0; i < 4000; ++i) {
      if (i == 50) {
        fs::remove(dump);
      }
      cloud = sinp::readLammpsTrj(dump.string(), 9, cloud);
    }
  } catch (const std::exception &) {
  }
  hold.store(false);
  holder.join();
  fs::remove_all(dir);
}

TEST_CASE("a signal during the index is retried and a short read stops",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  const auto dump = dir / "alarm.lammpstrj";
  writeDumpLines(dump, 80000, 80000);
  const auto path = dump.string();
  struct Handler {
    static void onAlarm(int) {}
  };
  const auto previous = std::signal(SIGALRM, &Handler::onAlarm);
  struct itimerval tick {};
  tick.it_interval.tv_usec = 200;
  tick.it_value.tv_usec = 200;
  REQUIRE(::setitimer(ITIMER_REAL, &tick, nullptr) == 0);
  std::atomic<bool> go{true};
  std::thread shrinker([&] {
    while (go.load()) {
      ::truncate(path.c_str(), 0);
      writeDumpLines(dump, 40, 40);
    }
  });
  int frames = -1;
  for (int i = 0; i < 40; ++i) {
    frames = sinp::nLammpsFrames(path);
    REQUIRE(frames >= 0);
    REQUIRE(frames <= 1);
  }
  go.store(false);
  shrinker.join();
  tick = {};
  REQUIRE(::setitimer(ITIMER_REAL, &tick, nullptr) == 0);
  std::signal(SIGALRM, previous);
  sinp::dropLammpsDumpIndex(path);
  fs::remove_all(dir);
}

TEST_CASE("a file with no timestep is not a dump", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  const auto dump = dir / "note.lammpstrj";
  writeText(dump, "not a dump\n");
  REQUIRE(sinp::nLammpsFrames(dump.string()) == 0);
  sinp::dropLammpsDumpIndex(dump.string());
  fs::remove_all(dir);
}

TEST_CASE("a dump truncated to nothing reports no frames", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  const auto dump = dir / "tail.lammpstrj";
  writeDumpLines(dump, 4, 4);
  const auto path = dump.string();
  REQUIRE(sinp::nLammpsFrames(path) == 1);
  REQUIRE(::truncate(path.c_str(), 0) == 0);
  REQUIRE(sinp::nLammpsFrames(path) == 0);
  sinp::dropLammpsDumpIndex(path);
  fs::remove_all(dir);
}

TEST_CASE("an empty dump reports no frames", "[seams_input]") {
  const auto path = dseams_test::freshFile("empty.lammpstrj");
  {
    std::ofstream out(path);
  }
  REQUIRE(sinp::nLammpsFrames(path.string()) == 0);
}

TEST_CASE("replacing a dump during the index does not fault",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  const auto dump = dir / "replace.lammpstrj";
  std::string body =
      "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n40\n"
      "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      "ITEM: ATOMS id type x y z\n";
  for (int i = 0; i < 40; ++i) {
    body += std::to_string(i);
    body += " 1 0 0 0\n";
  }
  writeText(dump, body);
  std::atomic<bool> go{true};
  std::thread indexer([&] {
    while (go.load()) {
      sinp::dropLammpsDumpIndex(dump.string());
      REQUIRE(sinp::nLammpsFrames(dump.string()) >= 0);
      molSys::PointCloud<molSys::Point<double>, double> local;
      local = sinp::readLammpsTrj(dump.string(), 1, local);
    }
  });
  for (int i = 0; i < 300; ++i) {
    fs::remove(dump);
    writeText(dump, body);
  }
  go.store(false);
  indexer.join();
  sinp::dropLammpsDumpIndex(dump.string());
  const int frames = sinp::nLammpsFrames(dump.string());
  REQUIRE(frames >= 0);
  fs::remove_all(dir);
}

TEST_CASE("a short dump unlinked during the read is not finished",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "short.lammpstrj";
  writeDumpLines(dump, 200000, 80000);
  const auto path = dump.string();
  std::atomic<bool> saw{false};
  std::thread unlinker([&] {
    for (int i = 0; i < 2000000 && !saw.load(); ++i) {
      if (fdsForPath(path).size() >= 1) {
        fs::remove(dump);
        saw.store(true);
        return;
      }
    }
  });
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path, 1, cloud);
  unlinker.join();
  REQUIRE(saw);
  REQUIRE(cloud.nop == 80000);
  sinp::dropLammpsDumpIndex(path);
  fs::remove_all(dir);
}

TEST_CASE("a finished short dump names the missing atoms", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "short.lammpstrj";
  writeDumpLines(dump, 4, 1);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  REQUIRE_THROWS_WITH(
      sinp::readLammpsTrj(dump.string(), 1, cloud),
      "LAMMPS dump frame is truncated: expected 4 atoms, read 1.");
  sinp::dropLammpsDumpIndex(dump.string());
  fs::remove_all(dir);
}

TEST_CASE("an atom row past the field cap is still one row", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "wide.lammpstrj";
  std::string row = "1 1";
  for (int i = 0; i < 40; ++i) {
    row += " 0";
  }
  writeText(dump,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type x y z\n" +
                row + "\n");
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(dump.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  REQUIRE(cloud.pts.size() == 1);
  sinp::dropLammpsDumpIndex(dump.string());
  fs::remove_all(dir);
}

TEST_CASE("readers survive a dump rewritten under them", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "stress.lammpstrj";
  writeDumpLines(dump, 2000, 2000);
  const auto path = dump.string();
  std::atomic<bool> go{true};
  std::atomic<int> reads{0};
  std::vector<std::thread> readers;
  for (int n = 0; n < 4; ++n) {
    readers.emplace_back([&] {
      try {
        molSys::PointCloud<molSys::Point<double>, double> local;
        while (go.load()) {
          local = sinp::readLammpsTrj(path, 1, local);
          reads.fetch_add(1);
        }
      } catch (const std::exception &) {
        reads.fetch_add(1);
      }
    });
  }
  const auto until = std::chrono::steady_clock::now() + std::chrono::milliseconds(200);
  while (std::chrono::steady_clock::now() < until) {
    writeDumpLines(dump, 2000, 2000);
  }
  go.store(false);
  for (auto &th : readers) {
    th.join();
  }
  REQUIRE(reads.load() > 0);
  const auto finished = dir / "finished.lammpstrj";
  writeDumpLines(finished, 3, 1);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  REQUIRE_THROWS_WITH(
      sinp::readLammpsTrj(finished.string(), 1, cloud),
      "LAMMPS dump frame is truncated: expected 3 atoms, read 1.");
  sinp::dropLammpsDumpIndex(path);
  sinp::dropLammpsDumpIndex(finished.string());
  fs::remove_all(dir);
}

TEST_CASE("a second reader closes a clone that is no longer a timestep",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "clone.lammpstrj";
  writeDumpLines(dump, 180000, 180000);
  const auto path = dump.string();
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path, 1, cloud);
  REQUIRE(cloud.nop == 180000);
  struct stat st {};
  REQUIRE(::stat(path.c_str(), &st) == 0);
  std::atomic<bool> go{true};
  std::thread holder([&] {
    try {
      molSys::PointCloud<molSys::Point<double>, double> local;
      while (go.load()) {
        local = sinp::readLammpsTrj(path, 1, local);
      }
    } catch (const std::exception &) {
    }
  });
  std::this_thread::sleep_for(std::chrono::milliseconds(30));
  const int w = ::open(path.c_str(), O_WRONLY);
  REQUIRE(w >= 0);
  const char junk[] = "NOT A TIMESTEP!\n";
  REQUIRE(::pwrite(w, junk, sizeof junk - 1, 0) ==
          static_cast<ssize_t>(sizeof junk - 1));
  ::close(w);
  struct timespec ts[2];
  ts[0] = st.st_atim;
  ts[1] = st.st_mtim;
  REQUIRE(::utimensat(AT_FDCWD, path.c_str(), ts, 0) == 0);
  cloud = sinp::readLammpsTrj(path, 1, cloud);
  REQUIRE(cloud.nop == 0);
  go.store(false);
  holder.join();
  sinp::dropLammpsDumpIndex(path);
  fs::remove_all(dir);
}

TEST_CASE("a dump atom row without coordinates is not an atom",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto noType = dir / "notype.lammpstrj";
  writeText(noType,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id y z\n"
            "1 0 0\n");
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(noType.string(), 1, cloud);
  REQUIRE(cloud.nop == 0);
  const auto noX = dir / "nox.lammpstrj";
  writeText(noX,
            "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
            "ITEM: ATOMS id type z\n"
            "1 1 0\n");
  cloud = sinp::readLammpsTrj(noX.string(), 1, cloud);
  REQUIRE(cloud.nop == 0);
  sinp::dropLammpsDumpIndex(noType.string());
  sinp::dropLammpsDumpIndex(noX.string());
  fs::remove_all(dir);
}

TEST_CASE("a dump that grows during the read keeps the atoms it has",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "grow.lammpstrj";
  writeDumpLines(dump, 80000, 20000);
  const auto path = dump.string();
  std::atomic<bool> reading{true};
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  std::thread reader([&] {
    try {
      cloud = sinp::readLammpsTrj(path, 1, cloud);
    } catch (const std::exception &) {
    }
    reading.store(false);
  });
  bool grew = false;
  for (int i = 0; i < 2000000 && reading.load(); ++i) {
    if (fdsForPath(path).size() >= 1) {
      std::ofstream out(path, std::ios::app);
      for (int n = 0; n < 20; ++n) {
        out << "1 1 0 0 0\n";
      }
      grew = true;
      break;
    }
  }
  reader.join();
  REQUIRE(grew);
  REQUIRE(cloud.nop > 0);
  sinp::dropLammpsDumpIndex(path);
  fs::remove_all(dir);
}

TEST_CASE("a blocked reader rejects a frame the index cannot hold",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "block.lammpstrj";
  writeDumpLines(dump, 120000, 120000);
  const auto path = dump.string();
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path, 1, cloud);
  REQUIRE(cloud.nop == 120000);
  std::atomic<bool> go{true};
  std::thread holder([&] {
    try {
      molSys::PointCloud<molSys::Point<double>, double> local;
      while (go.load()) {
        local = sinp::readLammpsTrj(path, 1, local);
      }
    } catch (const std::exception &) {
    }
  });
  const auto until = std::chrono::steady_clock::now() + std::chrono::milliseconds(200);
  while (std::chrono::steady_clock::now() < until) {
    cloud = sinp::readLammpsTrj(path, 0, cloud);
    cloud = sinp::readLammpsTrj(path, 50, cloud);
  }
  REQUIRE(cloud.nop == 0);
  go.store(false);
  holder.join();
  sinp::dropLammpsDumpIndex(path);
  fs::remove_all(dir);
}

TEST_CASE("a short dump that grows during the wait is kept",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  for (int attempt = 0; attempt < 8; ++attempt) {
    const auto dump = dir / ("w" + std::to_string(attempt) + ".lammpstrj");
    writeDumpLines(dump, 40, 6);
    const auto path = dump.string();
    std::atomic<bool> reading{true};
    std::thread reader([&] {
      try {
        molSys::PointCloud<molSys::Point<double>, double> local;
        local = sinp::readLammpsTrj(path, 1, local);
      } catch (const std::exception &) {
      }
      reading.store(false);
    });
    for (int i = 0; i < 2000000 && reading.load(); ++i) {
      if (fdsForPath(path).size() >= 1) {
        // The parse of six lines is over before this wait ends, so the
        // append lands in the reader's short stability wait.
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
        std::ofstream out(path, std::ios::app);
        out << "9 1 0 0 0\n";
        break;
      }
    }
    reader.join();
    sinp::dropLammpsDumpIndex(path);
  }
  fs::remove_all(dir);
}

TEST_CASE("a second open fails when the process is out of descriptors",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "nofile.lammpstrj";
  writeDumpLines(dump, 220000, 220000);
  const auto path = dump.string();
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path, 1, cloud);
  REQUIRE(cloud.nop == 220000);
  // The session keeps the descriptor from the first read. A later read
  // then serves the frame without opening, so the descriptor limit never
  // applies. Drop the session before the limit is set.
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
    maxfd = std::max(maxfd, fd);
    held.push_back(fd);
  }
  // No free descriptor remains, so the next open fails. The limit is
  // restored before the directory is removed: remove_all needs a
  // descriptor, and the tight limit is still in force until this block
  // ends.
  { struct Restore {
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

#ifdef SEAMS_HAS_READCON
#ifdef SEAMS_HAS_READCON_DB
TEST_CASE("a corpus that cannot be opened or holds a long frame",
          "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto locked = dir / "locked";
  fs::create_directories(locked);
  fs::permissions(locked, fs::perms::owner_read | fs::perms::owner_exec);
  REQUIRE(sinp::ingestConCorpus(locked.string(), 1, "missing.con") < 0);
  fs::permissions(locked, fs::perms::owner_all);

  const int nAtoms = 50000;
  std::string atoms;
  atoms.reserve(static_cast<std::size_t>(nAtoms) * 32);
  for (int i = 0; i < nAtoms; ++i) {
    atoms += "   1.0    2.0    3.0 1    ";
    atoms += std::to_string(i);
    atoms += "\n";
  }
  const auto huge = dir / "huge.con";
  std::string text = "Random Number Seed\n{\"con_spec_version\":2}\n"
                     "10.0\t10.0\t10.0\n90.0\t90.0\t90.0\n0 0\n218 0 1\n1\n";
  text += std::to_string(nAtoms);
  text += "\n63.546000\nCu\nCoordinates of Component 1\n";
  text += atoms;
  writeText(huge, text);
  const auto corpus = dir / "corpus";
  fs::create_directories(corpus);
  REQUIRE(sinp::ingestConCorpus(corpus.string(), 7, huge.string()) > 0);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readConCorpus(corpus.string(), "", "", 1, cloud);
  REQUIRE(cloud.nop == nAtoms);
  const auto bare = dir / "bare.con";
  writeText(bare,
            "Random Number Seed\n{\"con_spec_version\":2}\n"
            "10.0\t10.0\t10.0\n90.0\t90.0\t90.0\n0 0\n218 0 1\n1\n0\n"
            "63.546000\nCu\nCoordinates of Component 1\n");
  const auto bareCorpus = dir / "bare-corpus";
  fs::create_directories(bareCorpus);
  if (sinp::ingestConCorpus(bareCorpus.string(), 1, bare.string()) > 0) {
    cloud = sinp::readConCorpus(bareCorpus.string(), "", "", 1, cloud);
  }
  cloud = sinp::readCon(bare.string(), 1, cloud);
  const auto cut = dir / "cut.con";
  writeText(cut,
            "Random Number Seed\n{\"con_spec_version\":2}\n"
            "10.0\t10.0\t10.0\n90.0\t90.0\t90.0\n0 0\n218 0 1\n1\n1\n"
            "63.546000\nCu\nCoordinates of Component 1\n");
  cloud = sinp::readCon(cut.string(), 1, cloud);
  const auto nanCon = dir / "nan.con";
  writeText(nanCon,
            "Random Number Seed\n{\"con_spec_version\":2}\n"
            "10.0\t10.0\t10.0\n90.0\t90.0\t90.0\n0 0\n218 0 1\n1\n2\n"
            "63.546000\nCu\nCoordinates of Component 1\n"
            "   nan    2.0    3.0 1    0\n"
            "   nan    5.0    6.0 1    1\n");
  cloud = sinp::readCon(nanCon.string(), 1, cloud);
  fs::remove_all(dir);
}
#endif
#endif

TEST_CASE("a directory and a truncated timestep header", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(dir.string(), 1, cloud);
  const auto trunc = dir / "trunc.lammpstrj";
  writeText(trunc, "ITEM: TIMESTEP");
  cloud = sinp::readLammpsTrj(trunc.string(), 2, cloud);
  const auto junk = dir / "junk.con";
  writeText(junk, "not a con\n");
  const auto junkCorpus = dir / "junk-corpus";
  fs::create_directories(junkCorpus);
#ifdef SEAMS_HAS_READCON
  if (sinp::ingestConCorpus(junkCorpus.string(), 1, junk.string()) > 0) {
    cloud = sinp::readConCorpus(junkCorpus.string(), "", "", 1, cloud);
  }
#endif
  fs::remove_all(dir);
}

TEST_CASE("a short xyz names the missing atoms", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_short.xyz");
  {
    std::ofstream out(path);
    out << "3\ncomment\nO 0 0 0\n";
  }
  try {
    sinp::readXYZ(path.string());
    FAIL("a short xyz must not come back as a frame");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) ==
            "XYZ frame is truncated: expected 3 atoms, read 1.");
  }
  fs::remove(path);
}

TEST_CASE("an empty xyz names the header", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_empty.xyz");
  {
    std::ofstream out(path);
    out << "";
  }
  try {
    sinp::readXYZ(path.string());
    FAIL("an empty xyz must not come back as a frame");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) == "XYZ header is not an atom count.");
  }
  fs::remove(path);
}

TEST_CASE("a junk xyz header names the failure", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_junk.xyz");
  {
    std::ofstream out(path);
    out << "this is not an xyz file\n";
  }
  try {
    sinp::readXYZ(path.string());
    FAIL("a junk xyz must not come back as a frame");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) == "XYZ header is not an atom count.");
  }
  fs::remove(path);
}

TEST_CASE("a multi-frame xyz returns the first frame", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_multi.xyz");
  {
    std::ofstream out(path);
    out << "1\nframe\nO 1 2 3\n1\nframe\nO 4 5 6\n";
  }
  auto cloud = sinp::readXYZ(path.string());
  REQUIRE(cloud.nop == 1);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(1.0, 1e-9));
  fs::remove(path);
}

TEST_CASE("a short lammps dump names the missing atoms", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_short.lammpstrj");
  {
    std::ofstream out(path);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n3\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  try {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    sinp::readLammpsTrj(path.string(), 1, cloud);
    FAIL("a short dump must not come back as a frame");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) ==
            "LAMMPS dump frame is truncated: expected 3 atoms, read 1.");
  }
  fs::remove(path);
  sinp::dropLammpsDumpIndex(path.string());
}

TEST_CASE("an orthogonal lattice keeps three lengths", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_ortho.xyz");
  {
    std::ofstream out(path);
    out << "1\nLattice=\"20 0 0 0 18 0 0 0 16\"\nO 1 1 1\n";
  }
  auto cloud = sinp::readXYZ(path.string());
  REQUIRE(cloud.box.size() == 3);
  REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(20.0, 1e-9));
  REQUIRE_THAT(cloud.box[1], Catch::Matchers::WithinAbs(18.0, 1e-9));
  REQUIRE_THAT(cloud.box[2], Catch::Matchers::WithinAbs(16.0, 1e-9));
  REQUIRE_THAT(cloud.boxLow[0], Catch::Matchers::WithinAbs(0.0, 1e-9));
  fs::remove(path);
}

TEST_CASE("a broken lattice is named", "[seams_input]") {
  const auto bare = dseams_test::freshFile("dseams_bare.xyz");
  {
    std::ofstream out(bare);
    out << "1\nLattice=10 0 0 0 10 0 0 0 10\nO 0 0 0\n";
  }
  try {
    sinp::readXYZ(bare.string());
    FAIL("a lattice without quotes must fail");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) == "XYZ Lattice is not a cell.");
  }
  fs::remove(bare);
  const auto shear = dseams_test::freshFile("dseams_shear.xyz");
  {
    std::ofstream out(shear);
    out << "1\nLattice=\"10 1 0 0 10 0 0 0 10\"\nO 0 0 0\n";
  }
  try {
    sinp::readXYZ(shear.string());
    FAIL("a non-reduced lattice must fail");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) == "XYZ Lattice is not a reduced cell.");
  }
  fs::remove(shear);
}

TEST_CASE("scaled dump columns and a parallel short frame", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_scaled_cols.lammpstrj");
  {
    std::ofstream out(path);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z xu xsu yu ysu zu zsu\n"
        << "1 1 0.1 0.2 0.3 1 2 3 4 5 6\n";
  }
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud = sinp::readLammpsTrj(path.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(0.1, 1e-9));
  fs::remove(path);
  sinp::dropLammpsDumpIndex(path.string());

  const auto only = dseams_test::freshFile("dseams_xsu.lammpstrj");
  {
    std::ofstream out(only);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type xsu ysu zsu\n"
        << "1 1 4 5 6\n";
  }
  cloud = sinp::readLammpsTrj(only.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(4.0, 1e-9));
  fs::remove(only);
  sinp::dropLammpsDumpIndex(only.string());
}

TEST_CASE("reader edges on real files", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);

  auto expectThrow = [](const fs::path &path, const char *message) {
    try {
      sinp::readXYZ(path.string());
      FAIL("expected a named xyz failure");
    } catch (const std::runtime_error &err) {
      REQUIRE(std::string(err.what()) == message);
    }
  };

  const auto unclosed = dir / "unclosed.xyz";
  {
    std::ofstream out(unclosed);
    out << "1\nLattice=\"10 0 0 0 10 0 0 0 10\nO 0 0 0\n";
  }
  expectThrow(unclosed, "XYZ Lattice is not a cell.");

  const auto shortCell = dir / "shortcell.xyz";
  {
    std::ofstream out(shortCell);
    out << "1\nLattice=\"1 2 3\"\nO 0 0 0\n";
  }
  expectThrow(shortCell, "XYZ Lattice is not a cell.");

  const auto negative = dir / "neg.xyz";
  {
    std::ofstream out(negative);
    out << "1\nLattice=\"-1 0 0 0 10 0 0 0 10\"\nO 0 0 0\n";
  }
  expectThrow(negative, "XYZ Lattice is not a reduced cell.");

  const auto az = dir / "az.xyz";
  {
    std::ofstream out(az);
    out << "1\nLattice=\"10 0 1 0 10 0 0 0 10\"\nO 0 0 0\n";
  }
  expectThrow(az, "XYZ Lattice is not a reduced cell.");

  const auto bz = dir / "bz.xyz";
  {
    std::ofstream out(bz);
    out << "1\nLattice=\"10 0 0 0 10 1 0 0 10\"\nO 0 0 0\n";
  }
  expectThrow(bz, "XYZ Lattice is not a reduced cell.");

  const auto tiny = dir / "tiny.xyz";
  {
    std::ofstream out(tiny);
    out << "1\nLattice=\"0.5 0 0 0 0.5 0 0 0 0.5\"\nO 0 0 0\n";
  }
  auto cloud = sinp::readXYZ(tiny.string());
  REQUIRE(cloud.box.size() == 3);
  REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(0.5, 1e-9));

  const auto asDir = dir / "adir.xyz";
  fs::create_directories(asDir);
  expectThrow(asDir, "XYZ header is not an atom count.");

  const auto cr = dir / "cr.lammpstrj";
  {
    std::ofstream out(cr, std::ios::binary);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n\t1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\r\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z xu xs xsu yu ys ysu zu zs zsu\n"
        << "1 1 0.1 0.2 0.3 1 2 3 4 5 6 7 8 9\n";
  }
  molSys::PointCloud<molSys::Point<double>, double> dump;
  dump = sinp::readLammpsTrj(cr.string(), 1, dump);
  REQUIRE(dump.nop == 1);
  REQUIRE_THAT(dump.pts[0].x, Catch::Matchers::WithinAbs(0.1, 1e-9));
  sinp::dropLammpsDumpIndex(cr.string());

  const auto junkTok = dir / "junk.lammpstrj";
  {
    std::ofstream out(junkTok);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10 extra\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  dump = sinp::readLammpsTrj(junkTok.string(), 1, dump);
  REQUIRE(dump.nop == 1);
  sinp::dropLammpsDumpIndex(junkTok.string());

  const auto empty = dir / "empty.lammpstrj";
  {
    std::ofstream out(empty);
    out << "";
  }
  REQUIRE(sinp::nLammpsFrames(empty.string()) == 0);
  int seen = 0;
  sinp::forEachLammpsFrame(empty.string(), 1, 2, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; },
                           1, -1, 1);
  REQUIRE(seen == 0);
  sinp::dropLammpsDumpIndex(empty.string());

  const auto two = dir / "two.lammpstrj";
  {
    std::ofstream out(two);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        << "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 1 0 0\n";
  }
  seen = 0;
  sinp::forEachLammpsFrame(two.string(), 1, 2, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; },
                           2, 0, 1);
  REQUIRE(seen == 2);
  seen = 0;
  sinp::forEachLammpsFrame(two.string(), 1, 2, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; },
                           1, 5, 2);
  REQUIRE(seen == 0);
  sinp::dropLammpsDumpIndex(two.string());

  {
    std::ofstream out(two, std::ios::app);
    out << " \n";
  }
  dump = sinp::readLammpsTrj(two.string(), 1, dump);
  REQUIRE(dump.nop == 1);
  sinp::dropLammpsDumpIndex(two.string());

  const auto bonds = dir / "rings.txt";
  {
    std::ofstream out(bonds);
    out << "2 Bonds\n1 2 3 0\nabc\n1 2\n4 5 6 7\n";
  }
  auto rings = sinp::readBonds(bonds.string());
  REQUIRE(rings.size() >= 1);

  const auto missing = dir / "no-such.lammpstrj";
  dump = sinp::readLammpsTrj(missing.string(), 1, dump);
  REQUIRE(dump.nop == 0);

  fs::remove_all(dir);
}

TEST_CASE("one tilt axis is still a reduced cell", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const char *lattices[] = {"10 0 0 1 8 0 0 0 9", "10 0 0 0 8 0 1 0 9",
                            "10 0 0 0 8 0 0 1 9"};
  for (const char *lattice : lattices) {
    const auto path = dir / "one.xyz";
    {
      std::ofstream out(path);
      out << "1\nLattice=\"" << lattice << "\"\nO 0 0 0\n";
    }
    auto cloud = sinp::readXYZ(path.string());
    REQUIRE(cloud.box.size() == 6);
  }
  fs::remove_all(dir);
}

TEST_CASE("a lattice axis that is not positive is rejected", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  auto expectThrow = [](const fs::path &path) {
    try {
      sinp::readXYZ(path.string());
      FAIL("expected a reduced-cell failure");
    } catch (const std::runtime_error &err) {
      REQUIRE(std::string(err.what()) == "XYZ Lattice is not a reduced cell.");
    }
  };
  const char *bad[] = {"10 0 0 0 0 0 0 0 10", "10 1 0 0 10 0 0 0 10",
                       "10 0 0 0 10 0 0 0 0"};
  for (const char *lattice : bad) {
    const auto path = dir / "bad.xyz";
    {
      std::ofstream out(path);
      out << "1\nLattice=\"" << lattice << "\"\nO 0 0 0\n";
    }
    expectThrow(path);
  }
  fs::remove_all(dir);
}

TEST_CASE("an unreadable xyz and a partial symbol", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);

  const auto locked = dir / "locked.xyz";
  {
    std::ofstream out(locked);
    out << "1\ncomment\nO 0 0 0\n";
  }
  fs::permissions(locked, fs::perms::none);
  auto cloud = sinp::readXYZ(locked.string());
  REQUIRE(cloud.nop == 0);
  fs::permissions(locked, fs::perms::owner_all);

  const auto partial = dir / "partial.xyz";
  {
    std::ofstream out(partial);
    out << "1\ncomment\n1O 0 0 0\n";
  }
  cloud = sinp::readXYZ(partial.string());
  REQUIRE(cloud.nop == 1);

  const auto words = dir / "words.xyz";
  {
    std::ofstream out(words);
    out << "1\ncomment\nO x 0 0\n";
  }
  try {
    sinp::readXYZ(words.string());
    FAIL("expected a coordinate failure");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) == "XYZ coordinate is not a number.");
  }
  fs::remove_all(dir);
}

TEST_CASE("dump lines cover the field scanner", "[seams_input]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  molSys::PointCloud<molSys::Point<double>, double> cloud;

  const auto spaced = dir / "spaced.lammpstrj";
  {
    std::ofstream out(spaced);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n 1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0\t10\t0\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0.2 0 0\n";
  }
  cloud = sinp::readLammpsTrj(spaced.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(0.2, 1e-9));
  sinp::dropLammpsDumpIndex(spaced.string());

  const auto blankCount = dir / "blank.lammpstrj";
  {
    std::ofstream out(blankCount);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  cloud = sinp::readLammpsTrj(blankCount.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  sinp::dropLammpsDumpIndex(blankCount.string());

  const auto shortBox = dir / "shortbox.lammpstrj";
  {
    std::ofstream out(shortBox);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n5\nzz\t10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  cloud = sinp::readLammpsTrj(shortBox.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  sinp::dropLammpsDumpIndex(shortBox.string());

  const auto noCount = dir / "nocount.lammpstrj";
  {
    std::ofstream out(noCount);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM:\nITEM: FOO bar\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  cloud = sinp::readLammpsTrj(noCount.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  sinp::dropLammpsDumpIndex(noCount.string());

  const auto noZ = dir / "noz.lammpstrj";
  {
    std::ofstream out(noZ);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y\n1 1 0 0\n";
  }
  cloud = sinp::readLammpsTrj(noZ.string(), 1, cloud);
  REQUIRE(cloud.pts.empty());
  sinp::dropLammpsDumpIndex(noZ.string());

  const auto crlf = dir / "crlf.lammpstrj";
  {
    std::ofstream out(crlf, std::ios::binary);
    out << "ITEM: TIMESTEP\r\n0\r\nITEM: NUMBER OF ATOMS\r\n1\r\n"
        << "note ITEM: TIMESTEP embedded\r\n"
        << "ITEM: BOX BOUNDS pp pp pp\r\n0 10\r\n0 10\r\n0 10\r\n"
        << "ITEM: ATOMS id type x y z\r\n1 1 0 0 0\r\n"
        << "ITEM: TIMESTEP\r\n1\r\nITEM: NUMBER OF ATOMS\r\n1\r\n"
        << "ITEM: BOX BOUNDS pp pp pp\r\n0 10\r\n0 10\r\n0 10\r\n"
        << "ITEM: ATOMS id type x y z\r\n1 1 1 0 0\r\n"
        << "ITEM: TIME";
  }
  REQUIRE(sinp::nLammpsFrames(crlf.string()) >= 1);
  cloud = sinp::readLammpsTrj(crlf.string(), 2, cloud);
  REQUIRE(cloud.nop == 1);
  sinp::dropLammpsDumpIndex(crlf.string());

  const auto gap = dir / "gap.lammpstrj";
  {
    std::ofstream out(gap);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        << "\n";
  }
  cloud = sinp::readLammpsTrj(gap.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  cloud = sinp::readLammpsTrj(gap.string(), 2, cloud);
  REQUIRE(cloud.nop == 0);
  sinp::dropLammpsDumpIndex(gap.string());

  const auto flipped = dir / "flip.lammpstrj";
  {
    std::ofstream out(flipped);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  cloud = sinp::readLammpsTrj(flipped.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  const auto stamp = fs::last_write_time(flipped);
  {
    std::fstream io(flipped, std::ios::in | std::ios::out | std::ios::binary);
    io.seekp(0);
    io << "XXXX: TIMESTEP";
  }
  fs::last_write_time(flipped, stamp);
  cloud = sinp::readLammpsTrj(flipped.string(), 1, cloud);
  REQUIRE(cloud.nop == 0);
  sinp::dropLammpsDumpIndex(flipped.string());

  int seen = 0;
  sinp::forEachLammpsFrame(spaced.string(), 1, 1, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; },
                           1, -1, 1);
  REQUIRE(seen == 0);
  seen = 0;
  sinp::forEachLammpsFrame(spaced.string(), 1, 1, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; },
                           2, 0, 1);
  REQUIRE(seen == 1);
  sinp::dropLammpsDumpIndex(spaced.string());

  const auto bonds = dir / "rings.txt";
  {
    std::ofstream out(bonds);
    out << "2 Bonds\n\n1 2 3\n";
  }
  auto rings = sinp::readBonds(bonds.string());
  REQUIRE(rings.size() == 1);
  fs::permissions(bonds, fs::perms::none);
  rings = sinp::readBonds(bonds.string());
  REQUIRE(rings.empty());
  fs::permissions(bonds, fs::perms::owner_all);

  seen = 0;
  sinp::forEachLammpsFrame(spaced.string(), 1, 9, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; },
                           1, 0, 1);
  REQUIRE(seen == 1);
  for (int n = 0; n < 64; ++n) {
    cloud = sinp::readLammpsTrj(noZ.string(), 1, cloud);
  }
  sinp::dropLammpsDumpIndex(noZ.string());

  const auto noType = dir / "notype.lammpstrj";
  {
    std::ofstream out(noType);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id y z\n1 0 0\n"
        << "ITEM: ATOMS id type y z\n1 1 0 0\n"
        << "ITEM: ATOMS id type x z\n1 1 0 0\n"
        << "ITEM: ATOMS id type x y z\n1\n";
  }
  cloud = sinp::readLammpsTrj(noType.string(), 1, cloud);
  REQUIRE(cloud.pts.empty());
  sinp::dropLammpsDumpIndex(noType.string());

  const auto tabs = dir / "tabs.lammpstrj";
  {
    std::ofstream out(tabs);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n";
    for (int n = 0; n < 200; ++n) {
      out << "0\t10\t0\n";
    }
    out << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  cloud = sinp::readLammpsTrj(tabs.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  sinp::dropLammpsDumpIndex(tabs.string());

  const auto proseCount = dir / "prosecount.lammpstrj";
  {
    std::ofstream out(proseCount);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n"
        << "ITEM: BOX BOUNDS pp pp pp\nnotanumber\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 3 0 0\n";
  }
  const auto blanks = dir / "blanks.lammpstrj";
  {
    std::ofstream out(blanks);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n";
    for (int n = 0; n < 80; ++n) {
      out << " \n";
    }
  }
  cloud = sinp::readLammpsTrj(blanks.string(), 1, cloud);
  REQUIRE(cloud.nop == 0);
  sinp::dropLammpsDumpIndex(blanks.string());

  cloud = sinp::readLammpsTrj(proseCount.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(3.0, 1e-9));
  sinp::dropLammpsDumpIndex(proseCount.string());

  const auto follow = dir / "follow.lammpstrj";
  {
    std::ofstream out(follow);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        << "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 1 0 0\n";
  }
  cloud = sinp::readLammpsTrj(follow.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  {
    const auto stamp = fs::last_write_time(follow);
    std::fstream io(follow, std::ios::in | std::ios::out | std::ios::binary);
    std::string bytes((std::istreambuf_iterator<char>(io)), {});
    const auto at = bytes.find("\nITEM: TIMESTEP\n1\n");
    REQUIRE(at != std::string::npos);
    const auto pos = at + 1;
    io.clear();
    io.seekp(static_cast<std::streamoff>(pos));
    io << "XXXX: TIMESTEP";
    io.close();
    fs::last_write_time(follow, stamp);
  }
  cloud = sinp::readLammpsTrj(follow.string(), 2, cloud);
  REQUIRE(cloud.nop == 0);
  sinp::dropLammpsDumpIndex(follow.string());

  const auto bare = dir / "bare.xyz";
  {
    std::ofstream out(bare);
    out << "1\ncomment\n0 0 0\n";
  }
  cloud = sinp::readXYZ(bare.string());
  REQUIRE(cloud.nop == 1);
  {
    std::ofstream out(bare);
    out << "2\ncomment\nO 0 0 0\nO 1\n";
  }
  try {
    sinp::readXYZ(bare.string());
    FAIL("expected a short xyz");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()).find("XYZ frame is truncated") == 0);
  }

  const auto negativeCount = dir / "negcount.lammpstrj";
  {
    std::ofstream out(negativeCount);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n-1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS type x y z\n1 0 0 0\n1\n";
  }
  cloud = sinp::readLammpsTrj(negativeCount.string(), 1, cloud);
  REQUIRE(cloud.nop == 1);
  sinp::dropLammpsDumpIndex(negativeCount.string());

  seen = 0;
  sinp::forEachLammpsFrame(spaced.string(), 1, 1, -1,
                           [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { seen++; },
                           1, 0, 1);
  REQUIRE(seen == 1);

  const auto cut = dir / "cut.lammpstrj";
  {
    std::ofstream out(cut);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n"
        << "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 1 0 0\n";
  }
  try {
    sinp::forEachLammpsFrame(cut.string(), 1, 2, -1,
                             [&](int, molSys::PointCloud<molSys::Point<double>, double> &) {},
                             2, 0, 1);
    FAIL("expected a truncated dump");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()).find("LAMMPS dump frame is truncated") == 0);
  }
  sinp::dropLammpsDumpIndex(cut.string());

  fs::remove_all(dir);
}

TEST_CASE("a triclinic lattice is bounds plus tilt", "[seams_input]") {
  const auto path = dseams_test::freshFile("dseams_tilt.xyz");
  {
    std::ofstream out(path);
    // a=(10,0,0) b=(-2,8,0) c=(-1,-3,9)
    // xy=-2 xz=-1 yz=-3
    // xmin=-3 xmax=0 ymin=-3 ymax=0
    // span = (13, 11, 9), origin = (-3, -3, 0)
    out << "1\nLattice=\"10 0 0 -2 8 0 -1 -3 9\"\nO 1 1 1\n";
  }
  auto cloud = sinp::readXYZ(path.string());
  REQUIRE(cloud.box.size() == 6);
  REQUIRE(cloud.boxLow.size() == 3);
  REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(13.0, 1e-9));
  REQUIRE_THAT(cloud.box[1], Catch::Matchers::WithinAbs(11.0, 1e-9));
  REQUIRE_THAT(cloud.box[2], Catch::Matchers::WithinAbs(9.0, 1e-9));
  REQUIRE_THAT(cloud.box[3], Catch::Matchers::WithinAbs(-2.0, 1e-9));
  REQUIRE_THAT(cloud.box[4], Catch::Matchers::WithinAbs(-1.0, 1e-9));
  REQUIRE_THAT(cloud.box[5], Catch::Matchers::WithinAbs(-3.0, 1e-9));
  REQUIRE_THAT(cloud.boxLow[0], Catch::Matchers::WithinAbs(-3.0, 1e-9));
  REQUIRE_THAT(cloud.boxLow[1], Catch::Matchers::WithinAbs(-3.0, 1e-9));
  REQUIRE_THAT(cloud.boxLow[2], Catch::Matchers::WithinAbs(0.0, 1e-9));
  fs::remove(path);
}

TEST_CASE("a dump stream does not throw ios_base::failure", "[seams_input]") {
  // loadLammpsFrame calls parseLammpsFrameBody on an ifstream or an
  // OwnedFdStream. Neither enables exceptions(), and the body throws
  // std::runtime_error for a finished short frame. ios_base::failure
  // is therefore not a reachable catch around that call.
  std::ifstream bare;
  REQUIRE((bare.exceptions() & std::ios_base::badbit) == std::ios_base::goodbit);
  REQUIRE_NOTHROW(bare.exceptions(std::ios_base::goodbit));
  const auto dir = dseams_test::freshDir();
  const auto dump = dir / "short.lammpstrj";
  writeDumpLines(dump, 4, 1);
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  REQUIRE_THROWS_AS(sinp::readLammpsTrj(dump.string(), 1, cloud), std::runtime_error);
  fs::remove_all(dir);
}

TEST_CASE("a second xyz frame is past the one native frame", "[seams_input]") {
  const auto path = dseams_test::freshFile("one-frame.xyz");
  {
    std::ofstream out(path);
    out << "1\ncomment\nO 0 0 0\n";
  }
  REQUIRE_THROWS_AS(sinp::openFrame(path.string(), 2, -1), std::runtime_error);
  fs::remove(path);
}
