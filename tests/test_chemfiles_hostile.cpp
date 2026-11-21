// Copyright (c) d-SEAMS developers
// Every file under input/hostile/ crashed, hung or overflowed a buffer in
// the r2 chemfiles patches: the f1-f8 files are the gate reproducers, the
// fuzz_* files are minimised fuzzer findings named after the r2 source
// line or function that failed, or after what the file declares. Each one
// must now read or fail with a chemfiles error.
#include "fresh_temp.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <chemfiles.hpp>
#include <seams_input.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

namespace fs = std::filesystem;

struct Outcome {
  std::string error;
  std::size_t steps = 0;
  std::size_t natoms = 0;
};

// Opens PATH as FORMAT and reads every step.
Outcome readAll(const std::string &path, const std::string &format) {
  Outcome out;
  try {
    chemfiles::Trajectory trajectory(path, 'r', format);
    out.steps = trajectory.nsteps();
    for (std::size_t i = 0; i < out.steps; ++i) {
      out.natoms = trajectory.read().size();
    }
  } catch (const chemfiles::Error &err) {
    out.error = err.what();
  }
  return out;
}

struct Readable {
  const char *file;
  const char *format;
  std::size_t natoms;
};

struct Refused {
  const char *file;
  const char *format;
  const char *error;
};

void requireRead(const Readable &item) {
  const auto path = std::string("hostile/") + item.file;
  INFO(path << " as " << item.format);
  const auto out = readAll(path, item.format);
  REQUIRE(out.error.empty());
  REQUIRE(out.steps == 1);
  REQUIRE(out.natoms == item.natoms);
}

void requireRefused(const Refused &item) {
  const auto path = std::string("hostile/") + item.file;
  INFO(path << " as " << item.format);
  const auto out = readAll(path, item.format);
  REQUIRE(out.error.find(item.error) != std::string::npos);
}

constexpr const char *kCpmd = "CPMD 3.x output";
constexpr const char *kParm7 = "AMBER7 Parm";
constexpr const char *kTrjErr = "could not open the file";
constexpr const char *kTngErr = "critical error while calling";

} // namespace

TEST_CASE("F1: QM logs read their frame without a null QM timestep",
          "[chemfiles][molfile][hostile]") {
  // Molfile<F> passed a null QM timestep and never asked for its
  // metadata, so the GAMESS, CPMD and Gaussian readers wrote through
  // null when the log holds SCF or wavefunction data.
  for (const auto &item : {Readable{"f1_gamess_surface.log", "GAMESS", 1},
                           Readable{"f1_water_cpmd.out", kCpmd, 3},
                           Readable{"f1_water_cpmd3133.out", kCpmd, 3},
                           Readable{"f1_water_gaussian_mo.log",
                                    "Gaussian Logfile (g94,g98,g03)", 3}}) {
    requireRead(item);
  }
  // The extensions pick the same readers through seams.
  REQUIRE(sinp::openFrame("hostile/f1_water_cpmd.out", 1, -1).nop == 3);
  REQUIRE(sinp::openFrame("hostile/f1_gamess_surface.log", 1, -1).nop == 1);
}

TEST_CASE("F2-F8: the gate reproducers are errors", "[chemfiles][molfile][hostile]") {
  for (const auto &item : {
           // F2 a 4 kB version word overflowed have_cpmd's word buffers
           Refused{"f2_cpmd_long_version.out", kCpmd, "could not open"},
           // F3 %FORMAT search looped at EOF; the header strstr read past fgets
           Refused{"f3_truncated_flag.parm7", kParm7, "could not read the molecule"},
           Refused{"f3_truncated_header.parm7", kParm7, "could not open"},
           // F4 natoms = 2147483647 skipped lines past EOF forever
           Refused{"f4_natoms_max.cube", "Gaussian Cube", "could not open"},
           // F5 a row with fewer columns than the model column was indexed
           Refused{"f5_empty_type_symbol.mmcif", "mmCIF", "has no model number in column 14"},
           // F6 predict() printed a 16-byte field with no terminator
           Refused{"f6_one_byte.mae", "Maestro File", "could not open"},
           // F7 inbuf[strlen(inbuf) - 1] on a line starting with NUL
           Refused{"f7_nul.pqr", "PQR", "could not open"},
           // F8 3 * natoms overflowed before the size check
           Refused{"f8_natoms_wrap.binpos", "Scripps Binpos", "could not open"},
       }) {
    requireRefused(item);
  }
  for (const char *path : {"hostile/f3_truncated_flag.parm7", "hostile/f4_natoms_max.cube",
                           "hostile/f7_nul.pqr", "hostile/f8_natoms_wrap.binpos"}) {
    INFO(path);
    REQUIRE(sinp::chemfilesSteps(path) == -1);
    REQUIRE(sinp::openFrame(path, 1, -1).nop == 0);
  }
}

TEST_CASE("F5: an mmCIF atom with no name or type keeps its columns",
          "[chemfiles][hostile]") {
  // The writer printed empty type_symbol and label values, which left the
  // row one column short; that is how f5_empty_type_symbol.mmcif was made.
  const auto dir = dseams_test::freshDir();
  const auto path = (dir / "empty_type.mmcif").string();
  {
    chemfiles::Frame frame;
    frame.add_atom(chemfiles::Atom("", ""), chemfiles::Vector3D(1, 2, 3));
    frame.add_atom(chemfiles::Atom("O", "O"), chemfiles::Vector3D(4, 5, 6));
    chemfiles::Trajectory out(path, 'w', "mmCIF");
    out.write(frame);
  }
  chemfiles::Trajectory in(path, 'r', "mmCIF");
  REQUIRE(in.nsteps() == 1);
  const auto frame = in.read();
  REQUIRE(frame.size() == 2);
  REQUIRE(frame[0].name().empty());
  REQUIRE(frame[0].type().empty());
  REQUIRE(frame[1].name() == "O");
  REQUIRE_THAT(frame.positions()[1][2], Catch::Matchers::WithinAbs(6.0, 1e-4));
  fs::remove_all(dir);
}

TEST_CASE("fuzzed native and molfile inputs are errors", "[chemfiles][hostile]") {
  for (const auto &item : {
           Refused{"fuzz_cssr_cssr_130.cssr", "CSSR", "invalid atom id 0"},
           Refused{"fuzz_dcd_abort.dcd", "DCD", "invalid fortran record marker"},
           Refused{"fuzz_g96_gromacs_461.g96", "Gromacs g96", "could not open"},
           Refused{"fuzz_gamess_gamessplugin_3453.gms", "GAMESS", "could not open"},
           Refused{"fuzz_gro_tilted_v1.gro", "GRO", "v1(y), v1(z) and v2(z) must be zero"},
           Refused{"fuzz_mmtf_rle_odd.mmtf", "MMTF", "has an odd number of values"},
           Refused{"fuzz_mmtf_rle_overflow.mmtf", "MMTF", "exceed its length 4"},
           Refused{"fuzz_ncrst_netcdf3file_651.ncrst", "Amber Restart",
                   "uses dimension 1610612736, but there are only 5"},
           Refused{"fuzz_ncrst_unreachable_33.ncrst", "Amber Restart",
                   "variable data type 1769239296 is unknown"},
           Refused{"fuzz_tpr_abort.tpr", "TPR", "invalid GROMACS string"},
           Refused{"fuzz_tpr_binaryfile_254.tpr", "TPR", "invalid GROMACS string"},
           Refused{"fuzz_tpr_tpr_657.tpr", "TPR", "TPR symbol table"},
           Refused{"fuzz_tpr_tpr_1260.tpr", "TPR", "TPR atoms (293 entries of 12 bytes)"},
           Refused{"fuzz_trj_gromacs_1259_negative.trj", "TRJ", kTrjErr},
           Refused{"fuzz_trj_gromacs_1378.trj", "TRJ", kTrjErr},
           Refused{"fuzz_trr_abort.trr", "TRR", "has inconsistent sizes"},
           Refused{"fuzz_trr_binaryfile_254.trr", "TRR", "has inconsistent sizes"},
           Refused{"fuzz_trr_trr_108.trr", "TRR", "has inconsistent sizes"},
           Refused{"fuzz_xsf_xsfplugin_338.xsf", "(Animated) XCrySDen Structure File",
                   "could not open"},
           Refused{"fuzz_xtc_xdrfile_268.xtc", "XTC", "invalid small index 0"},
           Refused{"fuzz_xtc_xdrfile_316.xtc", "XTC", "invalid small index 0"},
           Refused{"fuzz_xtc_xdrfile_320.xtc", "XTC", "invalid small index 0"},
           Refused{"fuzz_xtc_xtc_108.xtc", "XTC", "buffer overrun during decompression"},
           Refused{"fuzz_xtc_xtc_123.xtc", "XTC", "buffer overrun during decompression"},
           Refused{"fuzz_dtr_dtrplugin_670.dtr", "DESRES Trajectory", "found no atoms"},
           Refused{"fuzz_dtr_dtrplugin_766.dtr", "DESRES Trajectory", "could not open"},
           Refused{"fuzz_dtr_dtrplugin_818_heap.dtr", "DESRES Trajectory", "found no atoms"},
           Refused{"fuzz_dtr_dtrplugin_818_misaligned.dtr", "DESRES Trajectory",
                   "found no atoms"},
           Refused{"fuzz_dtr_dtrplugin_818_segv.dtr", "DESRES Trajectory", "found no atoms"},
           Refused{"fuzz_dtr_stlvector_1144.dtr", "DESRES Trajectory", "found no atoms"},
           Refused{"fuzz_tng_tngio_2896.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_2923.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_2989.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_5060.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_5132.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_5258.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_6369.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_7489.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_tngio_7660.tng", "TNG", kTngErr},
       }) {
    requireRefused(item);
  }
}

TEST_CASE("fuzzed inputs that are valid files read", "[chemfiles][hostile]") {
  // A Maestro file with no bond block took &bonds[0] of an empty vector,
  // a CONECT row with 20 neighbours overflowed bonds[16], and a TRJ
  // header with a large natoms overflowed 3 * natoms.
  for (const auto &item : {Readable{"fuzz_maestro_no_bonds.mae", "Maestro File", 3},
                           Readable{"fuzz_xbgf_conect20.xbgf", "Internal Paratool Format", 3},
                           Readable{"fuzz_trj_gromacs_1259.trj", "TRJ", 10}}) {
    requireRead(item);
  }
}

TEST_CASE("a line longer than the text buffer keeps the later step offsets",
          "[chemfiles][hostile]") {
  // After the 8 kB line buffer grew, the file position moved by the
  // refilled size instead of the dropped bytes. Every later step offset
  // was wrong: the second XYZ frame failed to parse and SMI spun at EOF.
  const auto dir = dseams_test::freshDir();
  const auto xyz = (dir / "long_comment.xyz").string();
  {
    std::ofstream out(xyz, std::ios::binary);
    out << "1\n" << std::string(10000, 'c') << "\nO 1 2 3\n1\nsecond\nH 4 5 6\n";
  }
  chemfiles::Trajectory trajectory(xyz, 'r', "XYZ");
  REQUIRE(trajectory.nsteps() == 2);
  const auto second = trajectory.read_step(1);
  REQUIRE(second.size() == 1);
  REQUIRE(second[0].name() == "H");
  REQUIRE_THAT(second.positions()[0][0], Catch::Matchers::WithinAbs(4.0, 1e-6));
  const auto first = trajectory.read_step(0);
  REQUIRE(first[0].name() == "O");

  const auto smi = (dir / "long_line.smi").string();
  {
    std::ofstream out(smi, std::ios::binary);
    out << "C\n" << std::string(20000, 'C') << "\nC#N\nCC\n";
  }
  const auto out = readAll(smi, "SMI");
  REQUIRE(out.error.empty());
  REQUIRE(out.steps == 4);
  REQUIRE(out.natoms == 2);
  fs::remove_all(dir);
}

TEST_CASE("inputs from fuzzing the fixed r3 readers are errors", "[chemfiles][hostile]") {
  for (const auto &item : {
           // gradient rows past the atoms counted wrote outside the gradient
           Refused{"fuzz_gamess_gamessplugin_3793.gms", "GAMESS", "found no atoms"},
           // a coordinate row past the atoms counted overflowed the atom array
           Refused{"fuzz_gamess_strncpy.gms", "GAMESS", "could not open"},
           // the byte swap shifted a byte into the sign bit of an int
           Refused{"fuzz_js_endianswap_116.js", "js", "could not open"},
           // a missing ITEM line where NUMBER OF ATOMS belongs was dereferenced
           Refused{"fuzz_lammps_optional_556.lammpstrj", "LAMMPS",
                   "expected an ITEM entry for 'NUMBER OF ATOMS'"},
           // a CONECT record past the bonds counted wrote outside the bonds
           Refused{"fuzz_bgf_bgfplugin_400.bgf", "MSI Biograf Format",
                   "out of bounds atomic index"},
           // the atom counts of a line added up past INT_MAX
           Refused{"fuzz_xdatcar5_vasp5xdatcarplugin_80.VASP_XDATCAR5", "VASP_XDATCAR5",
                   "could not open"},
           // msgpack converted an empty transform matrix through a null pointer
           Refused{"fuzz_mmtf_carray_34.mmtf", "MMTF",
                   "MMTF transform matrix has 0 values instead of 16"},
           // freeing string data used the counts of a later block
           Refused{"fuzz_tng_tngio_8675.tng", "TNG", kTngErr},
           // string data was filled past the strings allocated
           Refused{"fuzz_tng_tngio_5231.tng", "TNG", kTngErr},
           // 100663311 particles declared, positions stored for 15
           Refused{"fuzz_tng_particles_unset.tng", "TNG", kTngErr},
       }) {
    requireRefused(item);
  }

  // A box shape block declaring 2^44 strings is rejected before anything
  // is allocated, and the positions still read.
  const auto tng = readAll("hostile/fuzz_tng_huge_strings.tng", "TNG");
  REQUIRE(tng.error.empty());
  REQUIRE(tng.steps == 10);
  REQUIRE(tng.natoms == 15);

  // The second open of the file in nsteps() parsed a NUL-holding line
  // differently from the first and passed a null handle to the plugin.
  const auto vtf = readAll("hostile/fuzz_vtf_vtfplugin_1326.vtf", "VTF trajectory format");
  REQUIRE(vtf.error.empty());
  REQUIRE(vtf.steps == 0);
}

TEST_CASE("a VTF line ending in a backslash continues on the next line",
          "[chemfiles][molfile][hostile]") {
  const auto dir = dseams_test::freshDir();
  const auto path = (dir / "continued.vtf").string();
  {
    std::ofstream out(path, std::ios::binary);
    out << "# water\natom 0 name O \\\n  type O\r\natom 1:2 name H type H\n"
        << "bond 0:1\nbond 0:2\nunitcell 10 10 10\n\ntimestep\n0 0 0\n1 0 0\n0 1 0\n"
        << "timestep\n0 0 1\n1 0 1\n0 1 1";
  }
  chemfiles::Trajectory trajectory(path, 'r', "VTF trajectory format");
  REQUIRE(trajectory.nsteps() == 2);
  trajectory.read();
  const auto frame = trajectory.read();
  REQUIRE(frame.size() == 3);
  REQUIRE(frame[0].type() == "O");
  REQUIRE(frame.topology().bonds().size() == 2);
  REQUIRE_THAT(frame.positions()[1][2], Catch::Matchers::WithinAbs(1.0, 1e-6));
  fs::remove_all(dir);
}
