// Copyright (c) d-SEAMS developers
// Reads one dump under LD_PRELOAD. The mode is argv[1]: eintr, eio, fstat.
#include <mol_sys.hpp>
#include <seams_input.hpp>

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
  const auto dir = fs::temp_directory_path() / "dseams-dump-fault";
  fs::remove_all(dir);
  fs::create_directories(dir);
  const auto dump = dir / "fault.lammpstrj";
  if (mode == "eintr" || mode == "eio") {
    writeDump(dump, 2, 2);
    const int n = sinp::nLammpsFrames(dump.string());
    sinp::dropLammpsDumpIndex(dump.string());
    fs::remove_all(dir);
    return n == 1 ? 0 : 1;
  }
  if (mode == "fstat") {
    writeDump(dump, 4, 1);
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    try {
      cloud = sinp::readLammpsTrj(dump.string(), 1, cloud);
    } catch (const std::exception &) {
      fs::remove_all(dir);
      return 1;
    }
    sinp::dropLammpsDumpIndex(dump.string());
    const bool ok = cloud.pts.size() == 1;
    fs::remove_all(dir);
    return ok ? 0 : 1;
  }
  return 2;
}
