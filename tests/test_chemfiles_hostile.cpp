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
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
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
           Refused{"fuzz_dtr_dtrplugin_670.dtr", "DESRES Trajectory",
                   "key records of 268435456 bytes instead of 24"},
           Refused{"fuzz_dtr_dtrplugin_766.dtr", "DESRES Trajectory",
                   "declares zero frames per file"},
           Refused{"fuzz_dtr_dtrplugin_818_heap.dtr", "DESRES Trajectory",
                   "key records of 268435456 bytes instead of 24"},
           Refused{"fuzz_dtr_dtrplugin_818_misaligned.dtr", "DESRES Trajectory",
                   "key records of 268435456 bytes instead of 24"},
           Refused{"fuzz_dtr_dtrplugin_818_segv.dtr", "DESRES Trajectory",
                   "key records of 268435456 bytes instead of 24"},
           Refused{"fuzz_dtr_stlvector_1144.dtr", "DESRES Trajectory",
                   "key records of 0 bytes instead of 24"},
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
           // an empty square brace was indexed at its first char
           Refused{"fuzz_smi_empty_brace.smi", "SMI", "empty square brace"},
           // renumbering function type 2147483647 overflowed an int
           Refused{"fuzz_tpr_tpr_776.tpr", "TPR", "invalid function type 2147483647"},
           // truncated TNG triplet-compressed block read past its buffer
           Refused{"fuzz_tng_coder_oob.tng", "TNG", kTngErr},
           // MMTF atomIdList had a duplicate atom id before the model's atomSkip_
           Refused{"fuzz_mmtf_atomskip.mmtf", "MMTF", "before current model offset"},
           // short CSSR unit cell line passed to string_view::substr(38)
           Refused{"fuzz_cssr_short_line.cssr", "CSSR", "too short"},
           // short SDF bond line passed to string_view::substr(6, 3)
           Refused{"fuzz_sdf_short_bond.sdf", "SDF", "bond line is too small"},
           // msgpack array32 declaring 536M elements in a 5-byte MMTF file
           Refused{"fuzz_mmtf_huge_array.mmtf", "MMTF", "array size overflow"},
           // zero-length msgpack string called memcpy with a null pointer
           Refused{"fuzz_mmtf_empty_str.mmtf", "MMTF", "Expected msgpack type to be MAP"},
           // compressed XTC delta overflowed int32_t
           Refused{"fuzz_xtc_int_overflow.xtc", "XTC", "integer overflow during decompression"},
           // truncated TNG header failed tng_util_trajectory_open after allocating handle
           Refused{"fuzz_tng_open_leak.tng", "TNG", kTngErr},
           // non-MOLDEN header leaked qmdata_t and FILE in open_molden_read
           Refused{"fuzz_molden_open_leak.molden", "Molden", "could not open"},
           // truncated InsightII car header leaked cardata and FILE
           Refused{"fuzz_car_open_leak.car", "InsightII car", "could not open"},
           // invalid CPMD first column leaked cpmddata and FILE
           Refused{"fuzz_cpmd_open_leak.traj", "CPMD", "could not open"},
           // ABINIT binary header with huge nsym hung reading past EOF
           Refused{"fuzz_abinit_bad_counts.bin", "ABINIT", "could not open"},
           // ABINIT binary header with out-of-range typat[i] indexed znucltypat OOB
           Refused{"fuzz_abinit_bad_typat.bin", "ABINIT", "could not open"},
           // AMBER Parm ReadPARM leaked parmstruct and its arrays
           Refused{"fuzz_parm_read_leak.parm", "AMBER Parm", "could not open"},
           // AMBER Parm perturbation file freed prm without clearing it before ~ReadPARM
           Refused{"fuzz_parm_ifpert.parm", "AMBER Parm", "could not open"},
           // VASP POSCAR leaked strdup(lineptr) and overwritten titleline
           Refused{"fuzz_poscar_strdup_leak.poscar", "VASP_POSCAR", "could not open"},
           // non-SQLite DMS file leaked sqlite3 handle on open failure
           Refused{"fuzz_dms_open_leak.dms", "DESRES Molecular Structure", "could not open"},
           // VTF atom record missing name leaked aid_list
           Refused{"fuzz_vtf_aid_leak.vsf", "VTF structure format", "could not open"},
           // MDF bond parse error leaked hasharray and atomnames
           Refused{"fuzz_mdf_bonds_leak.mdf", "InsightII MDF", "could not read bonds"},
           // repeated CPMD ATOMS block overwrote initatoms without freeing it
           Refused{"fuzz_cpmdlog_double_atoms.out", "CPMD 3.x output", "could not open"},
           // Amber NetCDF file declaring 536M atoms. The count does not fit
           // in the file, which is what stops the allocation.
           Refused{"fuzz_ncrst_huge_atoms.ncrst", "Amber Restart", "does not fit"},
           // CONECT line with no atom id passed empty string_view to decode_hybrid36
           Refused{"fuzz_pdb_short_conect.pdb", "PDB", "could not read atomic number"},
           // GAMESS basis set hitting EOF before TOTAL NUMBER leaked basis_set shells
           Refused{"fuzz_gamess_basis_leak.gms", "GAMESS", "could not open"},
           // TNG frame set with overflowing first_frame + n_frames
           Refused{"fuzz_tng_frameset_overflow.tng", "TNG", kTngErr},
           // XTC compressed coordinate thiscoord + minint overflowed int32_t
           Refused{"fuzz_xtc_minint_overflow.xtc", "XTC",
                   "integer overflow during decompression"},
           // LAMMPS Data Velocities section with out-of-range atom index
           Refused{"fuzz_lammpsdata_vel_oob.lmp", "LAMMPS Data",
                   "invalid atom index 2147483647"},
           // TNG compressed block claiming petabyte uncompressed_len before header check
           Refused{"fuzz_tng_uncompress_size.tng", "TNG", kTngErr},
           // TNG variable-atom frame set leaked molecule_cnt_list in tng_num_frame_sets_get
           Refused{"fuzz_tng_molcnt_leak.tng", "TNG", kTngErr},
           // TNG particle data allocation exceeded file size bound (24 TB calloc)
           Refused{"fuzz_tng_alloc_bound.tng", "TNG", kTngErr},
           // GAMESS read_shell_primitives hit EOF before freeing *prim
           Refused{"fuzz_gamess_prim_leak.gms", "GAMESS", "could not open"},
           // MMTF chainsPerModel/groupsPerChain sum exceeded numChains/numGroups
           Refused{"fuzz_mmtf_chain_oob.mmtf", "MMTF", "issue with data"},
           // DMS illegal bond returned without finalizing sqlite3_stmt, leaking db
           Refused{"fuzz_dms_illegal_bond_leak.dms", "DESRES Molecular Structure",
                   "could not read bonds"},
           // TNG cyclic/non-monotonic frame set navigation offsets caused infinite loop
           Refused{"fuzz_tng_frameset_cycle.tng", "TNG", kTngErr},
           // TNG next-frame-set header error leaked tng_gen_block and looped on bad offset
           Refused{"fuzz_tng_block_leak.tng", "TNG", kTngErr},
           // VTF atom id 4294967295 does not fit in the signed index
           Refused{"fuzz_vtf_huge_aid.vtf", "VTF structure format", "could not open"},
           // A declared id whose atom table would exceed 4 GiB
           Refused{"fuzz_vtf_atom_4gib.vtf", "VTF trajectory format", "could not open"},
           // binpos atom count would allocate about 8 GB; the file is shorter
           Refused{"fuzz_binpos_8gb_alloc.binpos", "Scripps Binpos", "could not open"},
           // XTC2 packed stream writes past the output or indexes magic_bits
           Refused{"fuzz_tng_xtc2_heap_write.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_xtc2_index_oob.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_xtc2_index_neg1123.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_xtc2_index_m1.tng", "TNG", kTngErr},
           // XTC3 auxiliary table larger than the coordinate block
           Refused{"fuzz_tng_xtc3_huge_table.tng", "TNG", kTngErr},
           // Streams written by the broken Huffman length sort are not decodable
           Refused{"fuzz_tng_bwlzh_roundtrip.tng", "TNG", kTngErr},
           Refused{"fuzz_tng_xtc3_roundtrip.tng", "TNG", kTngErr},
       }) {
    requireRefused(item);
  }

  // The second TRR frame declares atoms but no positions, velocities or
  // forces, which allocated its atom count without reading anything; like
  // any invalid trailing frame it now ends the trajectory.
  requireRead({"fuzz_trr_no_vectors.trr", "TRR", 10});
  // An '=' with no name before it looped forever in the comment parser.
  requireRead({"fuzz_xyz_extended_equals.xyz", "XYZ", 1});

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

  // A PDB CONECT record appearing before any ATOM/HETATM records must not
  // dereference an empty atom_offsets_ vector.
  const auto pdb = readAll("hostile/fuzz_pdb_conect_empty.pdb", "PDB");
  REQUIRE(pdb.error.empty());
  REQUIRE(pdb.steps == 1);
  REQUIRE(pdb.natoms == 0);

  // A nucleic residue with P and no O3' must not dereference end() in
  // link_standard_residue_bonds().
  const auto pdbNucleic = readAll("hostile/fuzz_pdb_nucleic_no_o3.pdb", "PDB");
  REQUIRE(pdbNucleic.error.empty());
  REQUIRE(pdbNucleic.steps == 1);
  REQUIRE(pdbNucleic.natoms == 2);

  // A Gromacs g96 file whose first timestep fails must free ts->pos in
  // mdio_timestep().
  const auto g96 = readAll("hostile/fuzz_g96_timestep_leak.g96", "Gromacs g96");
  REQUIRE(g96.error.empty());
  REQUIRE(g96.steps == 0);

  // A LAMMPS Data file with whitespace-only trailing comments must not
  // index an empty split() result.
  const auto lmp = readAll("hostile/fuzz_lammpsdata_blank_comment.lmp", "LAMMPS Data");
  REQUIRE(lmp.error.empty());
  REQUIRE(lmp.steps == 1);
  REQUIRE(lmp.natoms == 1);

  // A PDB HELIX record with helix type 0 must not underflow HELIX_TYPES[type - 1].
  const auto pdbHelix = readAll("hostile/fuzz_pdb_helix_zero.pdb", "PDB");
  REQUIRE(pdbHelix.error.empty());
  REQUIRE(pdbHelix.steps == 1);
  REQUIRE(pdbHelix.natoms == 1);

  // A CIF number with a 10-digit exponent must not overflow int and index
  // gemmi's power-of-10 LUT with a negative offset.
  const auto cifExp = readAll("hostile/fuzz_cif_huge_exponent.cif", "CIF");
  REQUIRE(cifExp.error.empty());
  REQUIRE(cifExp.steps == 1);
  REQUIRE(cifExp.natoms == 1);

  // A TRJ frame with no positions (mdts.pos == NULL) must not pass NULL to memcpy.
  const auto trjNullPos = readAll("hostile/fuzz_trj_null_pos.trj", "TRJ");
  REQUIRE(trjNullPos.error.empty());
  REQUIRE(trjNullPos.steps == 0);

  // Opening a nonexistent XBGF file must free bgf and bgf->meta before returning NULL.
  const auto missingXbgf = readAll("hostile/does_not_exist_xbgf.xbgf",
                                   "Internal Paratool Format");
  REQUIRE_FALSE(missingXbgf.error.empty());

  // An extended XYZ Properties repeat the atom line does not carry is an
  // error. The names are not allocated, and the properties already parsed
  // are not dropped on a later, valid repeat.
  const auto xyzRepeat = readAll("hostile/fuzz_xyz_huge_repeat.xyz", "XYZ");
  REQUIRE(xyzRepeat.error.find("do not match") != std::string::npos);
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

TEST_CASE("a VTF range past 200000 atoms is the structure it declares",
          "[chemfiles][molfile][hostile]") {
  const auto dir = dseams_test::freshDir();
  const auto path = (dir / "wide.vtf").string();
  for (int kAtoms : {201000, 330000, 1000002}) {
    {
      std::ofstream out(path, std::ios::binary);
      out << "atom 0:" << (kAtoms - 1) << " name O type O\n\ntimestep\n";
      for (int i = 0; i < kAtoms; ++i) {
        out << "1 2 3\n";
      }
    }
    chemfiles::Trajectory trajectory(path, 'r', "VTF trajectory format");
    REQUIRE(trajectory.nsteps() == 1);
    const auto frame = trajectory.read();
    REQUIRE(frame.size() == static_cast<std::size_t>(kAtoms));
    REQUIRE(frame[0].name() == "O");
    REQUIRE_THAT(frame.positions()[kAtoms - 1][0], Catch::Matchers::WithinAbs(1.0, 1e-6));
    REQUIRE_THAT(frame.positions()[kAtoms - 1][2], Catch::Matchers::WithinAbs(3.0, 1e-6));
  }
  fs::remove_all(dir);
}

TEST_CASE("TNG files written by the fixed encoder keep their coordinates",
          "[chemfiles][hostile]") {
  for (const char *name : {"fuzz_tng_codec_bwlzh.tng", "fuzz_tng_codec_xtc2.tng",
                           "fuzz_tng_codec_xtc3.tng", "fuzz_tng_codec_gzip.tng",
                           "fuzz_tng_codec_none.tng", "fuzz_tng_codec_triplet.tng"}) {
    chemfiles::Trajectory trajectory(std::string("hostile/") + name, 'r', "TNG");
    REQUIRE(trajectory.nsteps() == 3);
    const auto frame = trajectory.read();
    REQUIRE(frame.size() == 40);
    // Atom 0 is on the origin in y and z. x is the only jitter, and the
    // precision is 0.001 nm.
    REQUIRE_THAT(frame.positions()[0][1], Catch::Matchers::WithinAbs(0.0, 0.02));
    REQUIRE_THAT(frame.positions()[0][2], Catch::Matchers::WithinAbs(0.0, 0.02));
  }
  chemfiles::Trajectory moving("hostile/fuzz_tng_codec_bwlzh_move.tng", 'r', "TNG");
  REQUIRE(moving.nsteps() == 4);
  REQUIRE(moving.read().size() == 200);
}

namespace {

constexpr const char *kDtr = "DESRES Trajectory";

void writeBytes(const fs::path &path, const std::string &bytes) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  out << bytes;
}

std::string readBytes(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

std::string be32(std::uint32_t value) {
  return {static_cast<char>(value >> 24), static_cast<char>(value >> 16),
          static_cast<char>(value >> 8), static_cast<char>(value)};
}

// One timekey: the time, offset and size fields are each a low then a
// high big-endian word.
std::string timekey(std::uint32_t timeHi, std::uint32_t offset, std::uint32_t size) {
  return be32(0) + be32(timeHi) + be32(offset) + be32(0) + be32(size) + be32(0);
}

// hostile/dtr_clickme is a two-frame, eight-atom DESRES trajectory: 256
// frames per file, both 4096-byte frames in frame000000000, a metadata
// frame, not_hashed/.ddparams "0 0" and the empty clickme.dtr.
fs::path copyDtr(const fs::path &into, const std::string &name = "copy.dtr") {
  const auto dtr = into / name;
  fs::create_directories(dtr.parent_path());
  fs::copy("hostile/dtr_clickme", dtr, fs::copy_options::recursive);
  return dtr;
}

void requireTwoFrames(const fs::path &path) {
  INFO(path.string());
  chemfiles::Trajectory trajectory(path.string(), 'r', kDtr);
  REQUIRE(trajectory.nsteps() == 2);
  const auto frame = trajectory.read();
  REQUIRE(frame.size() == 8);
  REQUIRE_THAT(frame.positions()[0][0], Catch::Matchers::WithinAbs(0.0, 1e-4));
  REQUIRE(trajectory.read().size() == 8);
}

void requireDtrError(const fs::path &path, const std::string &expected) {
  INFO(path.string() << " expecting: " << expected);
  const auto out = readAll(path.string(), kDtr);
  REQUIRE(out.error.find(expected) != std::string::npos);
}

} // namespace

TEST_CASE("a DESRES trajectory is a valid directory or the clickme.dtr inside it",
          "[chemfiles][hostile]") {
  requireTwoFrames("hostile/dtr_clickme");
  requireTwoFrames("hostile/dtr_clickme/clickme.dtr");

  const auto dir = dseams_test::freshDir();
  SECTION("the old top-level .ddparams, or none, is a flat directory") {
    const auto dtr = copyDtr(dir);
    fs::remove(dtr / "not_hashed" / ".ddparams");
    writeBytes(dtr / ".ddparams", "0 0\n");
    requireTwoFrames(dtr / "clickme.dtr");
    fs::remove(dtr / ".ddparams");
    fs::remove(dtr / "metadata");
    requireTwoFrames(dtr);
  }
  SECTION("frames in the subdirectory .ddparams hashes them to") {
    // printf frame000000000 | cksum gives 3066110831; with 4 and 3
    // directories dtrplugin looks in 3066110831 % 4 = 3 and
    // (3066110831 / 4) % 3 = 2.
    const auto dtr = copyDtr(dir);
    writeBytes(dtr / "not_hashed" / ".ddparams", "4 3\n");
    fs::create_directories(dtr / "003" / "002");
    fs::rename(dtr / "frame000000000", dtr / "003" / "002" / "frame000000000");
    requireTwoFrames(dtr / "clickme.dtr");
    fs::rename(dtr / "003" / "002" / "frame000000000", dtr / "003" / "frame000000000");
    requireDtrError(dtr, "its first frame file 003/002/frame000000000 is missing");
  }
  SECTION("one frame per file") {
    const auto dtr = copyDtr(dir);
    const auto frames = readBytes(dtr / "frame000000000");
    writeBytes(dtr / "frame000000000", frames.substr(0, 4096));
    writeBytes(dtr / "frame000000001", frames.substr(4096));
    writeBytes(dtr / "timekeys", be32(0x4445534b) + be32(1) + be32(24) +
                                     timekey(0, 0, 4096) + timekey(0x3ff00000, 0, 4096));
    requireTwoFrames(dtr);
  }
  SECTION("an outer directory named like the marker") {
    requireTwoFrames(copyDtr(dir / "clickme.dtr.d", "run.dtr") / "clickme.dtr");
  }
  fs::remove_all(dir);
}

TEST_CASE("a clickme.dtr outside a valid DESRES directory is refused",
          "[chemfiles][hostile]") {
  const auto dir = dseams_test::freshDir();
  const auto lone = dir / "clickme.dtr";
  writeBytes(lone, "");
  requireDtrError(lone, "is a clickme.dtr outside a DESRES trajectory");
  requireDtrError(lone, "it has no timekeys file");

  // the marker next to a broken timekeys file is not a trajectory either
  writeBytes(dir / "timekeys", be32(0x4445534b) + be32(0) + be32(24));
  requireDtrError(lone, "is a clickme.dtr outside a DESRES trajectory");
  requireDtrError(lone, "declares zero frames per file");

  // the marker is empty; a clickme.dtr with content is not the marker
  const auto dtr = copyDtr(dir);
  writeBytes(dtr / "clickme.dtr", "x");
  requireDtrError(dtr / "clickme.dtr", "a DESRES trajectory is a directory");
  fs::remove_all(dir);
}

TEST_CASE("an empty .dtr file that is not the marker is empty",
          "[chemfiles][hostile]") {
  const auto dir = dseams_test::freshDir();
  for (const char *name : {"empty.dtr", "run.dtr", "notclickme.dtr", "clickme.dtr.dtr"}) {
    writeBytes(dir / name, "");
    requireDtrError(dir / name, "is empty");
  }
  // an empty .dtr is empty even inside a valid trajectory directory
  const auto dtr = copyDtr(dir);
  writeBytes(dtr / "other.dtr", "");
  requireDtrError(dtr / "other.dtr", "is empty");
  writeBytes(dir / "data.dtr", "not a directory");
  requireDtrError(dir / "data.dtr", "a DESRES trajectory is a directory");
  fs::remove_all(dir);
}

TEST_CASE("a DESRES directory is checked the way dtrplugin reads it",
          "[chemfiles][hostile]") {
  const auto dir = dseams_test::freshDir();
  const std::string prologue = be32(0x4445534b) + be32(256) + be32(24);
  const auto keys = timekey(0, 0, 4096) + timekey(0x3ff00000, 4096, 4096);
  const auto broken = [&](const char *what, const std::function<void(const fs::path &)> &edit,
                          const std::string &expected) {
    INFO(what);
    const auto dtr = copyDtr(dir, std::string(what) + ".dtr");
    edit(dtr);
    requireDtrError(dtr, "is not a DESRES trajectory: ");
    requireDtrError(dtr, expected);
    requireDtrError(dtr / "clickme.dtr", "is a clickme.dtr outside a DESRES trajectory");
    requireDtrError(dtr / "clickme.dtr", expected);
  };
  broken("no-timekeys", [](const fs::path &d) { fs::remove(d / "timekeys"); },
         "it has no timekeys file");
  broken("timekeys-dir",
         [](const fs::path &d) {
           fs::remove(d / "timekeys");
           fs::create_directory(d / "timekeys");
         },
         "it has no timekeys file");
  broken("short-prologue", [&](const fs::path &d) { writeBytes(d / "timekeys", "DESK"); },
         "shorter than the 12-byte prologue");
  broken("magic", [&](const fs::path &d) { writeBytes(d / "timekeys", "XESK" + prologue.substr(4) + keys); },
         "does not start with the DESRES timekeys magic number");
  broken("zero-fpf",
         [&](const fs::path &d) {
           writeBytes(d / "timekeys", be32(0x4445534b) + be32(0) + be32(24) + keys);
         },
         "declares zero frames per file");
  broken("key-size",
         [&](const fs::path &d) {
           writeBytes(d / "timekeys", be32(0x4445534b) + be32(256) + be32(23) + keys);
         },
         "declares key records of 23 bytes instead of 24");
  broken("partial-key",
         [&](const fs::path &d) { writeBytes(d / "timekeys", prologue + keys.substr(0, 47)); },
         "ends inside a key record");
  broken("no-frames", [&](const fs::path &d) { writeBytes(d / "timekeys", prologue); },
         "its timekeys file has no frames");
  broken("zero-size-key",
         [&](const fs::path &d) { writeBytes(d / "timekeys", prologue + timekey(0, 0, 0)); },
         "its first timekey does not describe a frame");
  broken("ddparams-text",
         [](const fs::path &d) { writeBytes(d / "not_hashed" / ".ddparams", "x\n"); },
         "its not_hashed/.ddparams is not two directory counts from 0 to 4096");
  broken("ddparams-range",
         [](const fs::path &d) { writeBytes(d / "not_hashed" / ".ddparams", "5000 0\n"); },
         "its not_hashed/.ddparams is not two directory counts from 0 to 4096");
  broken("ddparams-top",
         [](const fs::path &d) {
           fs::remove(d / "not_hashed" / ".ddparams");
           writeBytes(d / ".ddparams", "-1 0\n");
         },
         "its .ddparams is not two directory counts from 0 to 4096");
  broken("metadata", [](const fs::path &d) { writeBytes(d / "metadata", "x"); },
         "its metadata is not a DESRES frame");
  broken("no-frame", [](const fs::path &d) { fs::remove(d / "frame000000000"); },
         "its first frame file frame000000000 is missing");
  broken("short-frame",
         [](const fs::path &d) {
           writeBytes(d / "frame000000000", readBytes(d / "frame000000000").substr(0, 100));
         },
         "its first frame file is shorter than its timekey says");
  broken("frame-magic",
         [](const fs::path &d) {
           auto bytes = readBytes(d / "frame000000000");
           bytes[3] = 'X';
           writeBytes(d / "frame000000000", bytes);
         },
         "its first frame is not a DESRES frame");
  fs::remove_all(dir);
}

TEST_CASE("a hostile DESRES frame in a valid directory still fails in the plugin",
          "[chemfiles][hostile]") {
  // The fuzzed 818 frames carry the frame magic and a corrupt header. In a
  // directory that passes the structural checks they reach dtrplugin,
  // whose frame reader refuses them, leaving no atoms.
  const auto dir = dseams_test::freshDir();
  for (const char *name : {"fuzz_dtr_dtrplugin_818_heap.dtr", "fuzz_dtr_dtrplugin_818_misaligned.dtr",
                           "fuzz_dtr_dtrplugin_818_segv.dtr"}) {
    INFO(name);
    const auto dtr = copyDtr(dir, name);
    const auto frame = readBytes(fs::path("hostile") / name / "frame000000000");
    REQUIRE(frame.size() == 96);
    writeBytes(dtr / "frame000000000", frame);
    writeBytes(dtr / "timekeys", be32(0x4445534b) + be32(1) + be32(24) + timekey(0, 0, 96));
    const auto out = readAll((dtr / "clickme.dtr").string(), kDtr);
    REQUIRE(out.error.find("found no atoms") != std::string::npos);
  }
  fs::remove_all(dir);
}

TEST_CASE("extended XYZ loads every column of a wide property",
          "[chemfiles][hostile]") {
  const auto dir = dseams_test::freshDir();
  const auto soap = (dir / "soap.xyz").string();
  const auto pair = (dir / "pair.xyz").string();
  {
    std::ofstream out(soap, std::ios::binary);
    out << "4\nProperties=species:S:1:pos:R:3:soap:R:5000\n";
    for (int atom = 0; atom < 4; ++atom) {
      out << "O 0 0 0";
      for (int i = 0; i < 5000; ++i) {
        out << ' ' << (atom * 5000 + i);
      }
      out << '\n';
    }
  }
  {
    std::ofstream out(pair, std::ios::binary);
    out << "1\nProperties=species:S:1:pos:R:3:a:R:3000:b:R:3000\nH 1 2 3";
    for (int i = 0; i < 6000; ++i) {
      out << ' ' << i;
    }
    out << '\n';
  }
  chemfiles::Trajectory soapTraj(soap, 'r', "XYZ");
  const auto soapFrame = soapTraj.read();
  REQUIRE(soapFrame.size() == 4);
  const auto &soapProps = soapFrame[3].properties();
  REQUIRE(soapProps.size() == 5000);
  REQUIRE(soapProps.get("soap_0"));
  REQUIRE_FALSE(soapProps.get("soap_5000"));
  REQUIRE(soapProps.get("soap_4999").value().as_double() == 3.0 * 5000 + 4999);
  chemfiles::Trajectory pairTraj(pair, 'r', "XYZ");
  const auto pairFrame = pairTraj.read();
  const auto &pairProps = pairFrame[0].properties();
  REQUIRE(pairProps.size() == 6000);
  REQUIRE(pairProps.get("a_0"));
  REQUIRE(pairProps.get("a_2999"));
  REQUIRE(pairProps.get("b_0"));
  REQUIRE(pairProps.get("b_2999").value().as_double() == 5999.0);
  fs::remove_all(dir);
}
