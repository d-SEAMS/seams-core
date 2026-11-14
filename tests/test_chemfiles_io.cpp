// Copyright (c) d-SEAMS developers
#include "fresh_temp.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bop.hpp>
#include <chemfiles.hpp>
#include <franzblau.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>

#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;
using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

std::map<int, int> ringSizes(Cloud cloud, int type) {
  std::map<int, int> hist;
  if (cloud.nop <= 0) {
    return hist;
  }
  auto rows = nneigh::neighListO(3.5, cloud, type);
  auto index = nneigh::neighbourListByIndex(cloud, rows);
  for (const auto &ring : primitive::ringNetwork(index, 6)) {
    hist[static_cast<int>(ring.size())] += 1;
  }
  return hist;
}

std::map<int, int> iceCounts(Cloud cloud, int type) {
  std::map<int, int> hist;
  auto rows = nneigh::neighListO(3.5, cloud, type);
  chill::getCorrelPlus(cloud, rows, false);
  chill::getIceTypePlusNoPrint(cloud, rows, false);
  for (const auto &point : cloud.pts) {
    if (point.type == type) {
      hist[static_cast<int>(point.iceType)] += 1;
    }
  }
  return hist;
}

void writeAtoms(const fs::path &path, const std::string &format,
                const std::vector<std::array<double, 3>> &xyz,
                const std::array<double, 3> &lengths, const std::string &name) {
  chemfiles::Frame frame;
  frame.set_cell(chemfiles::UnitCell(
      chemfiles::Vector3D(lengths[0], lengths[1], lengths[2])));
  for (const auto &p : xyz) {
    frame.add_atom(chemfiles::Atom(name),
                   chemfiles::Vector3D(p[0], p[1], p[2]));
  }
  if (format == std::string("XTC")) {
    frame.set("xtc_precision", 1000000.0);
  }
  chemfiles::Trajectory out(path.string(), 'w', format);
  out.write(frame);
}

std::vector<std::array<double, 3>> oxygenCoords(const std::string &dump,
                                                int type) {
  Cloud native;
  native = sinp::readLammpsTrjO(dump, 1, native, type);
  std::vector<std::array<double, 3>> xyz;
  xyz.reserve(static_cast<std::size_t>(native.nop));
  for (const auto &point : native.pts) {
    xyz.push_back({point.x, point.y, point.z});
  }
  return xyz;
}

} // namespace

TEST_CASE("a caught exception keeps its text", "[chemfiles]") {
  try {
    throw std::runtime_error("XYZ header is not an atom count.");
  } catch (...) {
    REQUIRE(sinp::currentExceptionText() == "XYZ header is not an atom count.");
  }
  try {
    throw 1;
  } catch (...) {
    REQUIRE(sinp::currentExceptionText() == "unhandled error");
  }
}

TEST_CASE("chemfiles matches the native LAMMPS and XYZ readers",
          "[chemfiles]") {
  Cloud native;
  native = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, native, 2);
  Cloud via = sinp::openFrame("traj/exampleTraj.lammpstrj", 1, 2);
  REQUIRE(native.nop > 10);
  REQUIRE(via.nop == native.nop);
  REQUIRE(via.box.size() >= 3);
  for (int k = 0; k < 3; ++k) {
    REQUIRE_THAT(via.box[static_cast<std::size_t>(k)],
                 Catch::Matchers::WithinAbs(native.box[static_cast<std::size_t>(k)],
                                            1e-6));
  }
  for (int i = 0; i < native.nop; ++i) {
    REQUIRE(via.pts[static_cast<std::size_t>(i)].type ==
            native.pts[static_cast<std::size_t>(i)].type);
    REQUIRE(via.pts[static_cast<std::size_t>(i)].atomID ==
            native.pts[static_cast<std::size_t>(i)].atomID);
    REQUIRE_THAT(via.pts[static_cast<std::size_t>(i)].x,
                 Catch::Matchers::WithinAbs(native.pts[static_cast<std::size_t>(i)].x,
                                            1e-4));
    REQUIRE_THAT(via.pts[static_cast<std::size_t>(i)].y,
                 Catch::Matchers::WithinAbs(native.pts[static_cast<std::size_t>(i)].y,
                                            1e-4));
    REQUIRE_THAT(via.pts[static_cast<std::size_t>(i)].z,
                 Catch::Matchers::WithinAbs(native.pts[static_cast<std::size_t>(i)].z,
                                            1e-4));
  }
  REQUIRE(ringSizes(via, 2) == ringSizes(native, 2));
  REQUIRE(iceCounts(via, 2) == iceCounts(native, 2));

  const auto path = dseams_test::freshFile("dseams_chem_native.xyz");
  {
    std::ofstream out(path);
    out << native.nop << "\nbox\n";
    for (const auto &point : native.pts) {
      out << "O " << point.x << " " << point.y << " " << point.z << "\n";
    }
  }
  auto xyzNative = sinp::readXYZ(path.string());
  Cloud xyzChem;
  xyzChem = sinp::readChemfiles(path.string(), 1, xyzChem, -1);
  REQUIRE(xyzChem.nop == xyzNative.nop);
  REQUIRE(xyzChem.nop == native.nop);
  for (int i = 0; i < xyzNative.nop; ++i) {
    REQUIRE_THAT(xyzChem.pts[static_cast<std::size_t>(i)].x,
                 Catch::Matchers::WithinAbs(
                     xyzNative.pts[static_cast<std::size_t>(i)].x, 1e-4));
    REQUIRE_THAT(xyzChem.pts[static_cast<std::size_t>(i)].y,
                 Catch::Matchers::WithinAbs(
                     xyzNative.pts[static_cast<std::size_t>(i)].y, 1e-4));
    REQUIRE_THAT(xyzChem.pts[static_cast<std::size_t>(i)].z,
                 Catch::Matchers::WithinAbs(
                     xyzNative.pts[static_cast<std::size_t>(i)].z, 1e-4));
  }
  fs::remove(path);
}

TEST_CASE("XTC DCD and NetCDF carry the same census as the dump",
          "[chemfiles]") {
  const auto xyz = oxygenCoords("traj/exampleTraj.lammpstrj", 2);
  REQUIRE(xyz.size() > 10);
  Cloud native;
  native = sinp::readLammpsTrjO("traj/exampleTraj.lammpstrj", 1, native, 2);
  const std::array<double, 3> lengths = {
      native.box[0], native.box[1], native.box[2]};
  const auto nativeRings = ringSizes(native, 2);
  const auto dir = dseams_test::freshDir();
  const struct {
    const char *file;
    const char *format;
    double tol;
    bool census;
  } formats[] = {
      {"dseams_chem.dcd", "DCD", 1e-4, true},
      {"dseams_chem.xtc", "XTC", 1e-3, true},
      {"dseams_chem.nc", "Amber NetCDF", 1e-4, true},
      {"dseams_chem.gro", "GRO", 1e-2, false},
      {"dseams_chem.pdb", "PDB", 1e-2, false},
  };
  for (const auto &fmt : formats) {
    const auto path = dir / fmt.file;
    writeAtoms(path, fmt.format, xyz, lengths, "O");
    Cloud got = sinp::openFrame(path.string(), 1, -1);
    INFO(fmt.file);
    REQUIRE(got.nop == native.nop);
    for (int i = 0; i < got.nop; ++i) {
      REQUIRE_THAT(got.pts[static_cast<std::size_t>(i)].x,
                   Catch::Matchers::WithinAbs(
                       native.pts[static_cast<std::size_t>(i)].x, fmt.tol));
    }
    if (fmt.census) {
      REQUIRE(ringSizes(got, got.pts[0].type) == nativeRings);
      REQUIRE(iceCounts(got, got.pts[0].type) == iceCounts(native, 2));
    }
    fs::remove(path);
  }
}

TEST_CASE("water names and a selection split oxygen from hydrogen",
          "[chemfiles]") {
  const auto path = dseams_test::freshFile("dseams_chem_water.xyz");
  {
    std::ofstream out(path);
    out << "3\n\nO 0 0 0\nH 1 0 0\nHW 0 1 0\n";
  }
  Cloud all;
  sinp::ChemfilesRead request;
  request.frame = 1;
  request.typeFilter = -1;
  all = sinp::readChemfiles(path.string(), all, request);
  REQUIRE(all.nop == 3);
  REQUIRE(all.pts[0].type == 8);
  REQUIRE(all.pts[1].type == 1);
  REQUIRE(all.pts[2].type == 1);

  Cloud guessed = sinp::openFrame(path.string(), 1, 0);
  REQUIRE(guessed.nop == 1);
  REQUIRE(guessed.pts[0].type == 8);
  Cloud selected = sinp::openFrame(path.string(), 1, 8, "name O");
  REQUIRE(selected.nop == 1);
  REQUIRE(selected.pts[0].type == 8);

  request.selection = "name H";
  request.guessOxygen = false;
  Cloud hydrogens;
  hydrogens = sinp::readChemfiles(path.string(), hydrogens, request);
  REQUIRE(hydrogens.nop == 1);
  REQUIRE(hydrogens.pts[0].type == 1);

  request.selection = "name O";
  request.typeFilter = 8;
  Cloud both;
  both = sinp::readChemfiles(path.string(), both, request);
  REQUIRE(both.nop == 1);

  request.selection = "name Zz";
  Cloud none;
  none = sinp::readChemfiles(path.string(), none, request);
  REQUIRE(none.nop == 0);

  request.selection = "this is not a selection";
  request.typeFilter = -1;
  Cloud bad;
  bad = sinp::readChemfiles(path.string(), bad, request);
  REQUIRE(bad.nop == 0);

  request.selection.clear();
  request.frame = 0;
  Cloud early;
  early = sinp::readChemfiles(path.string(), early, request);
  REQUIRE(early.nop == 0);

  request.frame = 4;
  Cloud late;
  late = sinp::readChemfiles(path.string(), late, request);
  REQUIRE(late.nop == 0);
  fs::remove(path);
}

TEST_CASE("a numeric dump guesses type 2 and then type 1", "[chemfiles]") {
  const auto path = dseams_test::freshFile("dseams_chem_types.xyz");
  {
    std::ofstream out(path);
    out << "3\n\n1 0 0 0\n1 1 0 0\n3 2 0 0\n";
  }
  Cloud guessed = sinp::openFrame(path.string(), 1, 0);
  REQUIRE(guessed.nop == 2);
  REQUIRE(guessed.pts[0].type == 1);

  {
    std::ofstream out(path);
    out << "2\n\n2 0 0 0\n1 3 0 0\n";
  }
  guessed = sinp::openFrame(path.string(), 1, 0);
  REQUIRE(guessed.nop == 1);
  REQUIRE(guessed.pts[0].type == 2);

  Cloud missing = sinp::openFrame("no/such/trajectory.lammpstrj", 1, -1);
  REQUIRE(missing.nop == 0);
  const auto junkDump = dseams_test::freshFile("dseams_chem_junk.lammpstrj");
  const auto junkXyz = dseams_test::freshFile("dseams_chem_junk.xyz");
  {
    std::ofstream out(junkDump);
    out << "this is not a dump\n";
    std::ofstream xyz(junkXyz);
    xyz << "not an xyz either\n";
  }
  REQUIRE(sinp::openFrame(junkDump.string(), 1, -1).nop == 0);
  REQUIRE(sinp::openFrame(junkDump.string(), 1, 2).nop == 0);
  REQUIRE(sinp::openFrame(junkDump.string(), 1, 0).nop == 0);
  // Ids are 1..N, so the early native path does not apply. The coordinate
  // is not a number, so chemfiles refuses the frame and the type filter
  // falls through to the native oxygen reader.
  const auto oneToN = dseams_test::freshFile("dseams_chem_oneton.lammpstrj");
  {
    std::ofstream out(oneToN);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id type x y z\n1 2 not-a-coordinate 0 0\n";
  }
  REQUIRE(sinp::openFrame(oneToN.string(), 1, 2).nop == 0);
  fs::remove(oneToN);
  REQUIRE_THROWS_AS(sinp::openFrame(junkXyz.string(), 1, -1), std::runtime_error);
  fs::remove(junkDump);
  fs::remove(junkXyz);
  REQUIRE(sinp::chemfilesSteps("no/such/trajectory.lammpstrj") == -1);
  REQUIRE_FALSE(sinp::forEachChemfilesFrame(
      "no/such/trajectory.lammpstrj", 1, 1, -1, "", 0, 1,
      [&](int, Cloud &) {}));
  REQUIRE(sinp::openFrame("noextension", 1, -1).nop == 0);
  REQUIRE(sinp::chemfilesSteps(path.string()) == 1);
  fs::remove(path);
}

TEST_CASE("a range walks every chemfiles step and skips CON", "[chemfiles]") {
  const auto path = dseams_test::freshFile("dseams_chem_range.xyz");
  {
    std::ofstream out(path);
    out << "1\n\nO 0 0 0\n1\n\nO 1 0 0\n1\n\nO 2 0 0\n";
  }
  std::vector<int> seen;
  const bool walked = sinp::forEachChemfilesFrame(
      path.string(), 1, 3, -1, "", 0, 1,
      [&](int frame, Cloud &cloud) {
        seen.push_back(frame);
        REQUIRE(cloud.nop == 1);
      });
  REQUIRE(walked);
  REQUIRE(seen.size() == 3);
  seen.clear();
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 0, -1, "", 0, 1,
      [&](int frame, Cloud &) { seen.push_back(frame); }));
  REQUIRE(seen.size() == 3);

  seen.clear();
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 3, -1, "", 1, 2,
      [&](int frame, Cloud &) { seen.push_back(frame); }));
  REQUIRE(seen.size() == 1);
  REQUIRE(seen[0] == 2);

  int calls = 0;
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 1, -1, "", 3, 2,
      [&](int, Cloud &) { calls += 1; }));
  REQUIRE(calls == 0);
  calls = 0;
  REQUIRE(sinp::forEachChemfilesFrame(
      path.string(), 1, 1, -1, "", 0, 0,
      [&](int, Cloud &) { calls += 1; }));
  REQUIRE(calls == 1);

  REQUIRE_FALSE(sinp::forEachChemfilesFrame(
      "con/tiny_multi_cuh2.con", 1, 1, -1, "", 0, 1,
      [&](int, Cloud &) {}));
#ifdef SEAMS_HAS_READCON
  Cloud via = sinp::openFrame("con/tiny_multi_cuh2.con", 1, -1);
  Cloud direct;
  direct = sinp::readCon("con/tiny_multi_cuh2.con", 1, direct);
  REQUIRE(via.nop == direct.nop);
  REQUIRE(via.nop > 0);
#endif
  fs::remove(path);
}

TEST_CASE("a triclinic cell is stored with three tilt factors", "[chemfiles]") {
  const auto path = dseams_test::freshFile("dseams_chem_tilt.pdb");
  chemfiles::Frame frame;
  frame.set_cell(chemfiles::UnitCell(chemfiles::Vector3D(10, 10, 10),
                                     chemfiles::Vector3D(90, 90, 60)));
  frame.add_atom(chemfiles::Atom("O"), chemfiles::Vector3D(1, 1, 1));
  chemfiles::Trajectory out(path.string(), 'w', "PDB");
  out.write(frame);
  out.close();
  Cloud cloud = sinp::openFrame(path.string(), 1, -1);
  REQUIRE(cloud.nop == 1);
  REQUIRE(cloud.box.size() == 6);
  fs::remove(path);

  const auto pdb = dseams_test::freshFile("dseams_chem_res.pdb");
  {
    std::ofstream out(pdb);
    out << "ATOM      1  OW  SOL     7       1.000   2.000   3.000  1.00  0.00\n"
        << "END\n";
  }
  Cloud residue;
  residue = sinp::openFrame(pdb.string(), 1, -1);
  REQUIRE(residue.nop == 1);
  REQUIRE(residue.pts[0].molID == 7);
  fs::remove(pdb);
}

TEST_CASE("a dcd uses the sibling pdb names", "[chemfiles]") {
  const auto root = dseams_test::freshDir();
  const auto dcd = root / "dseams_chem_names.dcd";
  const auto pdb = root / "dseams_chem_names.pdb";
  {
    chemfiles::Frame frame;
    frame.set_cell(chemfiles::UnitCell(chemfiles::Vector3D(20, 20, 20)));
    frame.add_atom(chemfiles::Atom("X"), chemfiles::Vector3D(1, 1, 1));
    frame.add_atom(chemfiles::Atom("X"), chemfiles::Vector3D(2, 1, 1));
    chemfiles::Trajectory out(dcd.string(), 'w', "DCD");
    out.write(frame);
  }
  {
    std::ofstream out(pdb);
    out << "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1           1\n"
        << "ATOM      1  O   HOH A   1       1.000   1.000   1.000  1.00  0.00           O\n"
        << "ATOM      2  H1  HOH A   1       2.000   1.000   1.000  1.00  0.00           H\n"
        << "END\n";
  }
  Cloud oxygen = sinp::openFrame(dcd.string(), 1, 8);
  REQUIRE(oxygen.nop == 1);
  REQUIRE(oxygen.pts[0].type == 8);
  fs::remove(dcd);
  fs::remove(pdb);
}

TEST_CASE("an extended xyz lattice is the cell", "[chemfiles]") {
  const auto path = dseams_test::freshFile("dseams_chem_lattice.xyz");
  {
    std::ofstream out(path);
    out << "2\n"
        << "Properties=species:S:1:pos:R:3 Lattice=\"20 0 0 0 18 0 0 0 16\"\n"
        << "O 1 1 1\n"
        << "O 6 4 3\n";
  }
  Cloud via = sinp::openFrame(path.string(), 1, -1);
  REQUIRE(via.nop == 2);
  REQUIRE(via.box.size() >= 3);
  REQUIRE_THAT(via.box[0], Catch::Matchers::WithinAbs(20.0, 1e-6));
  REQUIRE_THAT(via.box[1], Catch::Matchers::WithinAbs(18.0, 1e-6));
  REQUIRE_THAT(via.box[2], Catch::Matchers::WithinAbs(16.0, 1e-6));
  REQUIRE_THAT(via.pts[1].x, Catch::Matchers::WithinAbs(6.0, 1e-6));
  fs::remove(path);
}

TEST_CASE("cif and mmcif round trip with an explicit format", "[chemfiles]") {
  const auto root = dseams_test::freshDir();
  const auto cif = root / "dseams_chem_round.cif";
  const auto mmcif = root / "dseams_chem_round.mmcif";
  {
    chemfiles::Frame frame;
    frame.set_cell(chemfiles::UnitCell(chemfiles::Vector3D(10, 11, 12)));
    frame.add_atom(chemfiles::Atom("O"), chemfiles::Vector3D(1, 2, 3));
    chemfiles::Trajectory cifOut(cif.string(), 'w', "CIF");
    cifOut.write(frame);
    chemfiles::Trajectory mmOut(mmcif.string(), 'w', "mmCIF");
    mmOut.write(frame);
  }
  for (const auto &path : {cif, mmcif}) {
    Cloud cloud = sinp::openFrame(path.string(), 1, -1);
    REQUIRE(cloud.nop == 1);
    REQUIRE_THAT(cloud.pts[0].x, Catch::Matchers::WithinAbs(1.0, 1e-4));
    REQUIRE_THAT(cloud.pts[0].y, Catch::Matchers::WithinAbs(2.0, 1e-4));
    REQUIRE_THAT(cloud.pts[0].z, Catch::Matchers::WithinAbs(3.0, 1e-4));
    REQUIRE(cloud.box.size() >= 3);
    REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(10.0, 1e-4));
    REQUIRE_THAT(cloud.box[1], Catch::Matchers::WithinAbs(11.0, 1e-4));
    REQUIRE_THAT(cloud.box[2], Catch::Matchers::WithinAbs(12.0, 1e-4));
    fs::remove(path);
  }
}

TEST_CASE("an empty xyz range and a con path name the failure", "[chemfiles]") {
  const auto empty = dseams_test::freshFile("dseams_empty_range.xyz");
  {
    std::ofstream out(empty);
    out << "";
  }
  int seen = 0;
  try {
    sinp::forEachChemfilesFrame(empty.string(), 1, 3, -1, "", 0, 1,
                                [&](int, molSys::PointCloud<molSys::Point<double>, double> &) {
                                  seen += 1;
                                });
    FAIL("an empty xyz range must fail");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) == "XYZ read failed: the file has no frames");
  }
  REQUIRE(seen == 0);
  fs::remove(empty);
#ifndef SEAMS_HAS_READCON
  try {
    sinp::openFrame((dseams_test::freshFile("missing.con")).string(), 1, -1);
    FAIL("a con path must fail in this build");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()).find("readcon-core") != std::string::npos);
  }
#endif
}

TEST_CASE("the ice XXI cif is one cell of 152 oxygens", "[chemfiles]") {
  std::string path;
  for (const char *candidate :
       {"traj/iceXXI_lee2026.cif", "traj/iceXXI_lee2026.lammpstrj",
        "input/traj/iceXXI_lee2026.cif",
        "input/traj/iceXXI_lee2026.lammpstrj"}) {
    if (fs::exists(candidate)) {
      path = candidate;
      break;
    }
  }
  REQUIRE(!path.empty());
  Cloud cloud = sinp::openFrame(path, 1, -1);
  REQUIRE(cloud.nop == 152);
  REQUIRE(cloud.box.size() >= 3);
  REQUIRE_THAT(cloud.box[0], Catch::Matchers::WithinAbs(20.1966, 1e-3));
  REQUIRE_THAT(cloud.box[1], Catch::Matchers::WithinAbs(20.1966, 1e-3));
  REQUIRE_THAT(cloud.box[2], Catch::Matchers::WithinAbs(7.8912, 1e-3));
}

TEST_CASE("water names and a rejected dump use the matching reader", "[chemfiles]") {
  const auto pdb = dseams_test::freshFile("dseams_names.pdb");
  {
    std::ofstream out(pdb);
    out << "ATOM      1  OW  SOL A   1       1.000   2.000   3.000  1.00  0.00           O\n"
        << "ATOM      2  HW  SOL A   1       1.100   2.000   3.000  1.00  0.00           H\n"
        << "ATOM      3  OH  MOL A   2       4.000   5.000   6.000  1.00  0.00           O\n"
        << "ATOM      4  HO  MOL A   2       4.100   5.000   6.000  1.00  0.00           H\n"
        << "END\n";
  }
  Cloud named = sinp::openFrame(pdb.string(), 1, -1);
  REQUIRE(named.nop == 4);
  REQUIRE(named.pts[0].type == 8);
  REQUIRE(named.pts[1].type == 1);
  fs::remove(pdb);

  const auto junk = dseams_test::freshFile("dseams_not_a_dump.lammpstrj");
  {
    std::ofstream out(junk);
    out << "this is not a dump\n";
  }
  Cloud recovered = sinp::openFrame(junk.string(), 1, -1);
  REQUIRE(recovered.nop == 0);
  fs::remove(junk);

  const auto xyz = dseams_test::freshFile("dseams_one.xyz");
  {
    std::ofstream out(xyz);
    out << "1\ncomment\nO 1 1 1\n";
  }
  int seen = 0;
  sinp::forEachChemfilesFrame(xyz.string(), 1, 0, -1, "", -1, 1,
                              [&](int, molSys::PointCloud<molSys::Point<double>, double> &) {
                                seen += 1;
                              });
  REQUIRE(seen == 0);
  fs::remove(xyz);

  const auto gro = dseams_test::freshFile("dseams_empty.gro");
  {
    std::ofstream out(gro);
    out << "";
  }
  int groFrames = 0;
  const bool walked = sinp::forEachChemfilesFrame(
      gro.string(), 1, 2, -1, "", 0, 1,
      [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { groFrames += 1; });
  REQUIRE((walked == false || groFrames == 0));
  fs::remove(gro);
}

TEST_CASE("name checks and the dump fallback take both sides", "[chemfiles]") {
  const auto pdb = dseams_test::freshFile("dseams_oh_ho.pdb");
  {
    std::ofstream out(pdb);
    out << "ATOM      1  OH  MOL A   1       1.000   1.000   1.000  1.00  0.00\n"
        << "ATOM      2  HO  MOL A   1       1.200   1.000   1.000  1.00  0.00\n"
        << "ATOM      3  XX  MOL A   2       3.000   1.000   1.000  1.00  0.00\n"
        << "END\n";
  }
  Cloud named = sinp::openFrame(pdb.string(), 1, -1);
  REQUIRE(named.nop >= 1);
  REQUIRE(named.pts[0].type == 8);
  fs::remove(pdb);

  const auto junk = dseams_test::freshFile("dseams_bad_steps_dir.pdb");
  fs::remove_all(junk);
  fs::create_directories(junk);
  REQUIRE(sinp::chemfilesSteps(junk.string()) < 0);
  fs::remove_all(junk);

  const auto emptyGro = dseams_test::freshFile("dseams_zero.gro");
  {
    std::ofstream out(emptyGro);
    out << "";
  }
  Cloud none = sinp::openFrame(emptyGro.string(), 1, -1);
  REQUIRE(none.nop == 0);
  int walked = 0;
  sinp::forEachChemfilesFrame(
      emptyGro.string(), 1, 1, -1, "", 0, 1,
      [&](int, molSys::PointCloud<molSys::Point<double>, double> &) { walked += 1; });
  REQUIRE(walked == 0);
  fs::remove(emptyGro);

  const auto dump = dseams_test::freshFile("dseams_only_type1.dump");
  {
    std::ofstream out(dump);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 1 0 0 0\n";
  }
  // type 0 is the guess path: chemfiles may read this, and if it cannot
  // the native reader keeps the type 1 atom after an empty type 2 pass.
  Cloud guessed = sinp::openFrame(dump.string(), 1, 0);
  REQUIRE(guessed.nop == 1);
  fs::remove(dump);
  sinp::dropLammpsDumpIndex(dump.string());

  const auto names = dseams_test::freshFile("dseams_oh_ho.xyz");
  {
    std::ofstream out(names);
    out << "3\ncomment\nOH 0 0 0\nHO 1 0 0\nAr 2 0 0\n";
  }
  Cloud tokens = sinp::openFrame(names.string(), 1, 0);
  REQUIRE(tokens.nop >= 1);
  fs::remove(names);

  const auto emptyXyz = dseams_test::freshFile("dseams_empty_walk.xyz");
  {
    std::ofstream out(emptyXyz);
    out << "";
  }
  auto expectEmpty = [&](const std::function<void()> &fn) {
    try {
      fn();
      FAIL("expected an empty xyz failure");
    } catch (const std::runtime_error &err) {
      REQUIRE(std::string(err.what()) == "XYZ read failed: the file has no frames");
    }
  };
  expectEmpty([&] {
    sinp::forEachChemfilesFrame(
        emptyXyz.string(), 1, 1, -1, "", 0, 1,
        [&](int, molSys::PointCloud<molSys::Point<double>, double> &) {});
  });
  expectEmpty([&] { sinp::readChemfiles(emptyXyz.string(), 1, guessed, -1); });
  fs::remove(emptyXyz);

  const auto badDump = dseams_test::freshFile("dseams_not_a_dump.dump");
  fs::remove_all(badDump);
  fs::create_directories(badDump);
  Cloud fallback = sinp::openFrame(badDump.string(), 1, 0);
  REQUIRE(fallback.nop == 0);
  fs::remove_all(badDump);

  const auto typed = dseams_test::freshFile("dseams_type2.dump");
  {
    std::ofstream out(typed);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
        << "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        << "ITEM: ATOMS id type x y z\n1 2 4 5 6\n";
  }
  Cloud typeTwo = sinp::openFrame(typed.string(), 1, 0);
  REQUIRE(typeTwo.nop == 1);
  REQUIRE(typeTwo.pts[0].type == 2);
  fs::remove(typed);
  sinp::dropLammpsDumpIndex(typed.string());

  const auto words = dseams_test::freshFile("dseams_tokens.xyz");
  {
    std::ofstream out(words);
    out << "6\ncomment\nO 0 0 0\nOW 1 0 0\nOH 2 0 0\n"
        << "H 3 0 0\nHW 4 0 0\nHO 5 0 0\n";
  }
  const auto mixed = dseams_test::freshFile("dseams_mixed.pdb");
  {
    std::ofstream out(mixed);
    out << "ATOM      1   OH MOL A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        << "ATOM      2   HO MOL A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        << "END\n";
  }
  Cloud mixedCloud = sinp::openFrame(mixed.string(), 1, -1);
  REQUIRE(mixedCloud.nop == 2);
  REQUIRE(mixedCloud.pts[0].type == 6);
  fs::remove(mixed);

  const auto smiles = dseams_test::freshFile("dseams_group.smi");
  {
    std::ofstream out(smiles);
    out << "O\n";
  }
  Cloud grouped = sinp::openFrame(smiles.string(), 1, -1);
  REQUIRE(grouped.nop == 1);
  REQUIRE(grouped.pts[0].molID == grouped.pts[0].atomID);
  fs::remove(smiles);

  Cloud spelled = sinp::openFrame(words.string(), 1, -1);
  REQUIRE(spelled.nop == 6);
  REQUIRE(spelled.pts[2].type == 8);
  REQUIRE(spelled.pts[3].type == 1);
  fs::remove(words);

  const auto picked = dseams_test::freshFile("dseams_select.xyz");
  {
    std::ofstream out(picked);
    out << "2\ncomment\nO 0 0 0\nH 1 0 0\n";
  }
  Cloud selected = sinp::openFrame(picked.string(), 1, 8, "name O");
  REQUIRE(selected.nop == 1);
  Cloud allNamed = sinp::openFrame(picked.string(), 1, -1, "name O");
  REQUIRE(allNamed.nop == 1);
  fs::remove(picked);

  const auto locked = dseams_test::freshFile("dseams_locked_fallback.xyz");
  {
    std::ofstream out(locked);
    out << "1\ncomment\nO 0 0 0\n";
  }
  fs::permissions(locked, fs::perms::none);
  Cloud hidden = sinp::openFrame(locked.string(), 1, -1);
  REQUIRE(hidden.nop == 0);
  fs::permissions(locked, fs::perms::owner_all);
  fs::remove(locked);

  const auto asDir = dseams_test::freshFile("dseams_dir_fallback.xyz");
  fs::remove_all(asDir);
  fs::create_directories(asDir);
  try {
    sinp::openFrame(asDir.string(), 1, -1);
    FAIL("a directory is not an xyz frame");
  } catch (const std::runtime_error &err) {
    REQUIRE(std::string(err.what()) == "XYZ header is not an atom count.");
  }
  fs::remove_all(asDir);
}

TEST_CASE("dump molecule ids follow the LAMMPS atom id", "[chemfiles]") {
  const auto dir = dseams_test::freshDir();
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto shuf = dir / "shuf.lammpstrj";
  {
    std::ofstream out(shuf);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id mol type x y z\n"
           "2 9 1 1 0 0\n"
           "1 4 2 0 0 0\n";
  }
  Cloud via = sinp::openFrame(shuf.string(), 1, -1);
  Cloud again = sinp::openFrame(shuf.string(), 1, -1);
  REQUIRE(again.nop == via.nop);
  REQUIRE(via.nop == 2);
  REQUIRE(via.pts[0].atomID == 1);
  REQUIRE(via.pts[0].molID == 4);
  REQUIRE(via.pts[1].atomID == 2);
  REQUIRE(via.pts[1].molID == 9);

  const auto gap = dir / "gap.dump";
  {
    std::ofstream out(gap);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id mol type x y z\n"
           "4 7 2 1 1 1\n"
           "1 9 1 0 0 0\n";
  }
  Cloud gapped = sinp::openFrame(gap.string(), 1, -1);
  REQUIRE(gapped.nop == 2);
  bool saw4 = false;
  for (const auto &pt : gapped.pts) {
    if (pt.atomID == 4) {
      saw4 = true;
      // mol 7 belongs to atom 4. mol 9 belongs to atom 1.
      // A swap of those mol ids fails this check.
      REQUIRE(pt.molID == 7);
      REQUIRE(pt.molID != 9);
      REQUIRE(pt.type == 2);
    }
    if (pt.atomID == 1) {
      REQUIRE(pt.molID == 9);
      REQUIRE(pt.molID != 7);
    }
  }
  REQUIRE(saw4);

  Cloud missing = sinp::openFrame((dir / "missing.lammpstrj").string(), 1, -1);
  REQUIRE(missing.nop == 0);

  const auto dup = dir / "dup.lammpstrj";
  {
    std::ofstream out(dup);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id mol type x y z\n"
           "1 3 2 0 0 0\n"
           "1 3 1 1 0 0\n";
  }
  setenv("SEAMS_MOL_INDEX_LOG", "1", 1);
  Cloud doubled = sinp::openFrame(dup.string(), 1, -1);
  const auto zero = dir / "zero-log.lammpstrj";
  {
    std::ofstream out(zero);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id type x y z\n0 1 0 0 0\n";
  }
  setenv("SEAMS_MOL_INDEX_LOG", "0", 1);
  Cloud loggedZero = sinp::openFrame(zero.string(), 1, -1);
  const auto emptyLog = dir / "empty-log.lammpstrj";
  fs::copy_file(zero, emptyLog);
  setenv("SEAMS_MOL_INDEX_LOG", "", 1);
  Cloud loggedEmpty = sinp::openFrame(emptyLog.string(), 1, -1);
  Cloud frameZero = sinp::openFrame(shuf.string(), 0, -1);
  unsetenv("SEAMS_MOL_INDEX_LOG");
  REQUIRE(frameZero.nop == 0);
  {
    std::ofstream out(shuf, std::ios::app);
    out << "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n1\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id mol type x y z\n1 4 2 0 0 0\n";
  }
  Cloud grown = sinp::openFrame(shuf.string(), 1, -1);
  REQUIRE(grown.pts[0].molID == 4);
  const auto stamped = fs::last_write_time(shuf);
  fs::last_write_time(shuf, stamped + std::chrono::seconds(5));
  Cloud later = sinp::openFrame(shuf.string(), 1, -1);
  REQUIRE(later.pts[0].atomID == 1);
  fs::last_write_time(shuf, stamped + std::chrono::seconds(5) + std::chrono::nanoseconds(2000));
  Cloud laterNs = sinp::openFrame(shuf.string(), 1, -1);
  REQUIRE(laterNs.pts[0].molID == 4);
  fs::remove_all(dir);
}

TEST_CASE("openFrame keeps the reduced sH edge, not the dump span", "[chemfiles]") {
  // BOX BOUNDS is the bounding box. xy = -6.2106057 pads x, so the
  // edge is 24.8424228 rather than the 31.053 span of the bounds.
  Cloud cloud = sinp::openFrame("traj/genice_sH.lammpstrj", 1, 1);
  REQUIRE(cloud.nop == 136);
  double cell[3][3];
  double origin[3];
  nneigh::dumpBoundsToH(cloud.box, cloud.boxLow, cell, origin);
  REQUIRE_THAT(cell[0][0], Catch::Matchers::WithinAbs(24.8424228, 1e-6));
  REQUIRE_THAT(cell[1][1], Catch::Matchers::WithinAbs(10.7570846, 1e-6));
  REQUIRE_THAT(cell[2][2], Catch::Matchers::WithinAbs(20.1312910, 1e-6));
  REQUIRE_THAT(cell[1][0], Catch::Matchers::WithinAbs(-6.2106057, 1e-6));
  REQUIRE_THAT(origin[0], Catch::Matchers::WithinAbs(0.0, 1e-6));
}
