// Copyright (c) d-SEAMS developers
//
// Decides whether an installed chemfiles can replace the patched 0.10.4.
//
// It first lists the registered formats. cmake/patches/
// chemfiles-0.10.4-molfiles.patch builds every VMD molfile plugin that
// reads atoms and registers 37 molfile formats. An unpatched 0.10.x
// registers only TRJ, PSF and Molden. Any missing name exits 7 and prints
// how many are present and which are missing, so the build never keeps a
// library with fewer readers. It then counts the frames of
// cmake/molfile_water.rst7. Unpatched Molfile<F>::nsteps() hands rst7 a
// null timestep and the probe dies; the patched count is one frame.
//
// Then it reads cmake/mmcif_no_entity.mmcif. That file is what chemfiles
// 0.10.4 writes for one oxygen with no residue: label_seq_id is "." and the
// atom_site loop has no label_entity_id. The 0.10.4 reader treats "." as
// "parse the entity id" and dereferences the missing column. A fresh read
// throws, and a write followed by a read in the same process segfaults
// inside mmCIFFormat::read. Exit 0 only when the oxygen comes back at
// (1, 2, 3) inside a 10 x 11 x 12 cell.
//
// It then counts the steps of a MOL2 molecule with an empty counts line.
// 0.10.4 indexes the empty split and segfaults, so the probe dies; a
// patched reader throws a format error, which passes.
//
// Last, it reads a LAMMPS Data file whose header holds a line that is
// neither a keyword nor a section name. 0.10.4 guards that with an assert,
// which a release build compiles out, and then spins in unreachable().
// A patched reader throws a format error. An alarm turns the spin into a
// signal death, so the probe fails instead of hanging the configure step.
// Windows has no alarm and skips this step.
#include <chemfiles.hpp>

#include <cmath>
#include <iostream>
#include <set>
#include <string>

#if !defined(_WIN32)
#include <unistd.h>
#endif

#ifndef SEAMS_MMCIF_PROBE
#define SEAMS_MMCIF_PROBE ""
#endif
#ifndef SEAMS_RST7_PROBE
#define SEAMS_RST7_PROBE ""
#endif

// The molfile format names the patched library registers. Keep in step
// with cmake/patches/chemfiles-0.10.4-molfiles.patch and the reader table
// in docs/orgmode/reference/data-formats.org.
static const char *const molfileFormats[] = {
    "TRJ",
    "PSF",
    "Molden",
    "ABINIT",
    "MSI Biograf Format",
    "Scripps Binpos",
    "InsightII car",
    "CHARMM Coordinates",
    "CPMD",
    "Gaussian Cube",
    "DL_POLY_C HISTORY",
    "DL_POLY_4 HISTORY",
    "DL_POLY CONFIG",
    "DESRES Trajectory",
    "GAMESS",
    "js",
    "Maestro File",
    "InsightII MDF",
    "NAMD Binary Coordinates",
    "AMBER7 Parm",
    "AMBER Parm",
    "PQR",
    "AMBER7 Restart",
    "VMD Tinker",
    "VASP_XDATCAR5",
    "VASP_OUTCAR",
    "VASP_POSCAR",
    "VASP_XDATCAR",
    "VASP_xml",
    "VTF structure format",
    "VTF trajectory format",
    "Internal Paratool Format",
    "(Animated) XCrySDen Structure File",
    "DESRES Molecular Structure",
    "Gromacs g96",
    "CPMD 3.x output",
    "Gaussian Logfile (g94,g98,g03)",
};

int main(int argc, char **argv) {
  const char *path =
      (argc > 1 && argv[1][0] != '\0') ? argv[1] : SEAMS_MMCIF_PROBE;
  const char *rst7 =
      (argc > 2 && argv[2][0] != '\0') ? argv[2] : SEAMS_RST7_PROBE;
  if (path[0] == '\0' || rst7[0] == '\0') {
    std::cerr << "mmCIF or rst7 probe path is empty\n";
    return 2;
  }

  std::set<std::string> registered;
  for (const auto &format : chemfiles::formats_list()) {
    registered.insert(format.get().name);
  }
  const auto expected = sizeof molfileFormats / sizeof molfileFormats[0];
  std::size_t present = 0;
  std::string missing;
  for (const char *name : molfileFormats) {
    if (registered.count(name) != 0) {
      present += 1;
    } else {
      missing += missing.empty() ? name : std::string(", ") + name;
    }
  }
  if (present != expected) {
    std::cout << present << " of " << expected
              << " VMD molfile readers; missing " << missing << '\n';
    return 7;
  }
  try {
    chemfiles::Trajectory water(rst7, 'r', "AMBER7 Restart");
    if (water.nsteps() != 1 || water.read().size() != 3) {
      return 8;
    }
  } catch (const std::exception &err) {
    std::cerr << err.what() << '\n';
    return 8;
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
  } catch (const std::exception &err) {
    std::cerr << err.what() << '\n';
    return 1;
  }
  static const char mol2[] = "\n@<TRIPOS>MOLECULE\n";
  try {
    auto truncated = chemfiles::Trajectory::memory_reader(mol2, sizeof mol2 - 1, "MOL2");
    (void)truncated.nsteps();
  } catch (const std::exception &) {
  }
#if !defined(_WIN32)
  static const char lammps[] =
      "bad header\n\n1 atoms\n0.0 0.0 0.0 abc origin\n\nAtoms\n\n";
  alarm(20);
  try {
    auto bad = chemfiles::Trajectory::memory_reader(lammps, sizeof lammps - 1,
                                                    "LAMMPS Data");
    (void)bad.read();
  } catch (const std::exception &) {
  }
  alarm(0);
#endif
  return 0;
}
