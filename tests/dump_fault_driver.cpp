// Copyright (c) d-SEAMS developers
// Reads one dump under LD_PRELOAD. The mode is argv[1]: eintr, eio, fstat.
#include "fresh_temp.hpp"

#include <mol_sys.hpp>
#include <seams_input.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

static void writeDump(const fs::path &path, int declared, int lines) {
  std::ofstream out(path);
  out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n" << declared
      << "\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
      << "ITEM: ATOMS id type x y z\n";
  for (int i = 0; i < lines; ++i) {
    out << i << " 1 0 0 0\n";
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    return 2;
  }
  const std::string mode = argv[1];
  // The three modes run as parallel tests of one binary.
  const auto dir = dseams_test::freshDir();
  const auto dump = dir / "fault.lammpstrj";
  if (mode == "eintr" || mode == "eio") {
    writeDump(dump, 2, 2);
    const int n = sinp::nLammpsFrames(dump.string());
    sinp::dropLammpsDumpIndex(dump.string());
    fs::remove_all(dir);
    return n == 1 ? 0 : 1;
  }
  if (mode == "fstat") {
    // The fault library fails every fstat after the first while
    // DUMP_FSTAT_EIO is set. libc++'s remove_all calls fstat, so the
    // fault is cleared before the directory is removed.
    const auto cleanup = [&dir] {
      ::unsetenv("DUMP_FSTAT_EIO");
      fs::remove_all(dir);
    };
    writeDump(dump, 4, 1);
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    try {
      cloud = sinp::readLammpsTrj(dump.string(), 1, cloud);
    } catch (const std::exception &) {
      cleanup();
      return 1;
    }
    sinp::dropLammpsDumpIndex(dump.string());
    const bool ok = cloud.pts.size() == 1;
    cleanup();
    return ok ? 0 : 1;
  }
  return 2;
}
