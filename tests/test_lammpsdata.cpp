// Copyright (c) d-SEAMS developers
// A generated LAMMPS data file of three water molecules. Chemfiles names
// the format "LAMMPS Data" and registers no extension, so .lammpsdata is
// readable only when the reader is told that name. The same file is
// counted, read through the C API, and read by the seams CLI.

#include "fresh_temp.hpp"
#include "seams.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace {

namespace fs = std::filesystem;

struct Atom {
  int type;
  double x;
  double y;
  double z;
};

// Three molecules, atom_style full. Type 1 is oxygen, type 2 is hydrogen.
constexpr Atom kAtoms[9] = {
    {1, 1.0, 1.0, 1.0},  {2, 2.0, 1.0, 1.0},  {2, 1.0, 2.0, 1.0},
    {1, 6.0, 1.0, 1.0},  {2, 7.0, 1.0, 1.0},  {2, 6.0, 2.0, 1.0},
    {1, 11.0, 1.0, 1.0}, {2, 12.0, 1.0, 1.0}, {2, 11.0, 2.0, 1.0},
};

constexpr double kBox = 16.0;

const char *kFile =
    "three water molecules\n"
    "\n"
    "9 atoms\n"
    "2 atom types\n"
    "6 bonds\n"
    "1 bond types\n"
    "\n"
    "0.0 16.0 xlo xhi\n"
    "0.0 16.0 ylo yhi\n"
    "0.0 16.0 zlo zhi\n"
    "\n"
    "Masses\n"
    "\n"
    "1 15.9994\n"
    "2 1.00794\n"
    "\n"
    "Atoms # full\n"
    "\n"
    "1 1 1 -0.834 1.0 1.0 1.0\n"
    "2 1 2 0.417 2.0 1.0 1.0\n"
    "3 1 2 0.417 1.0 2.0 1.0\n"
    "4 2 1 -0.834 6.0 1.0 1.0\n"
    "5 2 2 0.417 7.0 1.0 1.0\n"
    "6 2 2 0.417 6.0 2.0 1.0\n"
    "7 3 1 -0.834 11.0 1.0 1.0\n"
    "8 3 2 0.417 12.0 1.0 1.0\n"
    "9 3 2 0.417 11.0 2.0 1.0\n"
    "\n"
    "Bonds\n"
    "\n"
    "1 1 1 2\n"
    "2 1 1 3\n"
    "3 1 4 5\n"
    "4 1 4 6\n"
    "5 1 7 8\n"
    "6 1 7 9\n";

struct CliResult {
  int status = -1;
  std::string out;
};

CliResult cliRead(const fs::path &path) {
  CliResult cli;
  const char *bin = std::getenv("SEAMS_BIN");
  if (bin == nullptr || bin[0] == '\0') {
    cli.out = "SEAMS_BIN is unset";
    return cli;
  }
  const std::string command = "\"" + std::string(bin) +
                              "\" --type -1 read \"" + path.string() +
                              "\" 2>&1";
#if defined(_WIN32)
  FILE *pipe = _popen(command.c_str(), "r");
#else
  FILE *pipe = popen(command.c_str(), "r");
#endif
  if (pipe == nullptr) {
    cli.out = "CLI failed to start";
    return cli;
  }
  char buffer[1024];
  while (std::fgets(buffer, sizeof buffer, pipe) != nullptr) {
    cli.out += buffer;
  }
#if defined(_WIN32)
  cli.status = _pclose(pipe);
#else
  cli.status = pclose(pipe);
#endif
  return cli;
}

void expectWater(const fs::path &path) {
  size_t nframes = 99;
  const seams_status counted = seams_path_nframes(path.c_str(), &nframes);
  const char *countError = seams_last_error();
  const std::string countText =
      countError == nullptr ? std::string() : std::string(countError);
  INFO(countText);
  REQUIRE(counted == SEAMS_OK);
  REQUIRE(nframes == 1);

  seams_system *sys = seams_system_create();
  REQUIRE(sys != nullptr);
  const seams_status read =
      seams_system_read(sys, path.c_str(), 1, -1);
  const char *readError = seams_last_error();
  const std::string readText =
      readError == nullptr ? std::string() : std::string(readError);
  INFO(readText);
  REQUIRE(read == SEAMS_OK);
  REQUIRE(seams_system_count(sys) == 9);

  double xyz[27] = {};
  int types[9] = {};
  int ids[9] = {};
  double box[9] = {};
  REQUIRE(seams_system_positions(sys, xyz, 27) == SEAMS_OK);
  REQUIRE(seams_system_types(sys, types, 9) == SEAMS_OK);
  REQUIRE(seams_system_ids(sys, ids, 9) == SEAMS_OK);
  REQUIRE(seams_system_box(sys, box, 9) == SEAMS_OK);
  int oxygens = 0;
  int hydrogens = 0;
  for (int i = 0; i < 9; ++i) {
    REQUIRE(ids[i] == i + 1);
    REQUIRE(types[i] == kAtoms[i].type);
    if (types[i] == 1) {
      oxygens += 1;
    }
    if (types[i] == 2) {
      hydrogens += 1;
    }
    REQUIRE_THAT(xyz[3 * i],
                 Catch::Matchers::WithinAbs(kAtoms[i].x, 1e-6));
    REQUIRE_THAT(xyz[3 * i + 1],
                 Catch::Matchers::WithinAbs(kAtoms[i].y, 1e-6));
    REQUIRE_THAT(xyz[3 * i + 2],
                 Catch::Matchers::WithinAbs(kAtoms[i].z, 1e-6));
  }
  REQUIRE(oxygens == 3);
  REQUIRE(hydrogens == 6);
  for (int i = 0; i < 3; ++i) {
    REQUIRE_THAT(box[i], Catch::Matchers::WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(box[3 + i], Catch::Matchers::WithinAbs(kBox, 1e-6));
    REQUIRE_THAT(box[6 + i], Catch::Matchers::WithinAbs(0.0, 1e-6));
  }
  seams_system_destroy(sys);

  const CliResult cli = cliRead(path);
  INFO(cli.out);
  REQUIRE(cli.status == 0);
  REQUIRE(cli.out.find("nop 9") != std::string::npos);
  REQUIRE(cli.out.find("frame 1") != std::string::npos);
  REQUIRE(cli.out.find("box 16 16 16") != std::string::npos);
}

} // namespace

TEST_CASE("three water molecules in a lammps data file are one frame",
          "[lammpsdata]") {
  const fs::path dir = dseams_test::freshDir();
  const fs::path lower = dir / "water.lammpsdata";
  const fs::path upper = dir / "water.LAMMPSDATA";
  {
    std::ofstream out(lower);
    REQUIRE(out.good());
    out << kFile;
  }
  {
    std::ofstream out(upper);
    REQUIRE(out.good());
    out << kFile;
  }
  expectWater(lower);
  expectWater(upper);
  std::error_code ec;
  fs::remove_all(dir, ec);
}
