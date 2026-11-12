// Copyright (c) d-SEAMS developers
//
// Reads cmake/mmcif_no_entity.mmcif. That file is what chemfiles 0.10.4
// writes for one oxygen with no residue: label_seq_id is "." and the
// atom_site loop has no label_entity_id. The 0.10.4 reader treats "." as
// "parse the entity id" and dereferences the missing column. A fresh read
// throws, and a write followed by a read in the same process segfaults
// inside mmCIFFormat::read. Exit 0 only when the oxygen comes back at
// (1, 2, 3) inside a 10 x 11 x 12 cell.
#include <chemfiles.hpp>

#include <cmath>
#include <iostream>

#ifndef SEAMS_MMCIF_PROBE
#define SEAMS_MMCIF_PROBE ""
#endif

int main(int argc, char **argv) {
  const char *path =
      (argc > 1 && argv[1][0] != '\0') ? argv[1] : SEAMS_MMCIF_PROBE;
  if (path[0] == '\0') {
    std::cerr << "mmCIF probe path is empty\n";
    return 2;
  }
  try {
    chemfiles::Trajectory in(path, 'r', "mmCIF");
    if (in.nsteps() < 1) {
      return 3;
    }
    const auto frame = in.read();
    if (frame.size() != 1) {
      return 4;
    }
    const auto position = frame.positions()[0];
    const auto lengths = frame.cell().lengths();
    if (std::fabs(position[0] - 1.0) > 1e-3 ||
        std::fabs(position[1] - 2.0) > 1e-3 ||
        std::fabs(position[2] - 3.0) > 1e-3) {
      return 5;
    }
    if (std::fabs(lengths[0] - 10.0) > 1e-3 ||
        std::fabs(lengths[1] - 11.0) > 1e-3 ||
        std::fabs(lengths[2] - 12.0) > 1e-3) {
      return 6;
    }
    return 0;
  } catch (const std::exception &err) {
    std::cerr << err.what() << '\n';
    return 1;
  }
}
