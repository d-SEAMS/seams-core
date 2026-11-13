//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//-----------------------------------------------------------------------------------

// Walk a LAMMPS dump with the threaded CHILL+, q6 and shell kernels.
// `interface` writes layer profiles, front positions and the front velocity.
// `nucleation` writes the largest ice cluster and the critical-size estimate.

#include <interface_order.hpp>
#include <science_front.hpp>
#include <seams_input.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

struct Args {
  std::string cmd;
  std::string dump;
  std::string outDir = "paper/tex/figures/data/science";
  double ps = 10.0;
  double bin = 1.0;
  int threads = 0;
  int type = 1;
};

void usage() {
  std::cerr << "science-front interface DUMP [--ps 10] [--bin 1] [--out-dir DIR] [--threads N]\n"
            << "science-front nucleation DUMP [--out-dir DIR] [--threads N] [--type 1]\n";
}

std::string need(int &i, int argc, char **argv, const char *name) {
  if (i + 1 >= argc) {
    throw std::runtime_error(std::string("missing value for ") + name);
  }
  return argv[++i];
}

bool parseArgs(int argc, char **argv, Args &args) {
  std::vector<std::string> pos;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    if (key == "--ps") {
      args.ps = std::stod(need(i, argc, argv, "--ps"));
    } else if (key == "--bin") {
      args.bin = std::stod(need(i, argc, argv, "--bin"));
    } else if (key == "--out-dir") {
      args.outDir = need(i, argc, argv, "--out-dir");
    } else if (key == "--threads") {
      args.threads = std::stoi(need(i, argc, argv, "--threads"));
    } else if (key == "--type") {
      args.type = std::stoi(need(i, argc, argv, "--type"));
    } else if (key == "--help") {
      usage();
      return false;
    } else if (!key.empty() && key[0] == '-') {
      throw std::runtime_error("unknown option " + key);
    } else {
      pos.push_back(key);
    }
  }
  if (pos.size() != 2) {
    usage();
    return false;
  }
  args.cmd = pos[0];
  args.dump = pos[1];
  return true;
}

std::vector<long long> readTimesteps(const std::string &path) {
  std::ifstream in(path);
  std::vector<long long> steps;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) != 0) {
      continue;
    }
    std::string next;
    if (!std::getline(in, next)) {
      break;
    }
    steps.push_back(std::stoll(next));
  }
  return steps;
}

std::ofstream openOut(const std::string &path) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write " + path);
  }
  out.setf(std::ios::fixed);
  out.precision(6);
  return out;
}

void writeAtoms(std::ostream &os, const std::vector<science::AtomRow> &rows, int frame,
                const std::vector<double> &box) {
  os << "frame id x y z chill chill_ovito shell q6 q6bar in_extended\n";
  if (box.size() >= 3) {
    os << "# box " << box[0] << " " << box[1] << " " << box[2] << "\n";
  }
  for (const auto &row : rows) {
    os << frame << " " << row.id << " " << row.x << " " << row.y << " " << row.z << " "
       << row.chill << " " << row.chillOvito << " " << row.shell << " " << row.q6 << " "
       << row.q6bar << " " << row.inExtended << "\n";
  }
}

int runInterface(const Args &args) {
  const auto steps = readTimesteps(args.dump);
  std::vector<science::Layers> frames;
  std::vector<double> box;
  sinp::forEachLammpsFrame(
      args.dump, 1, 0, args.type,
      [&](int, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
        if (box.empty()) {
          box = cloud.box;
        }
        frames.push_back(science::layersOf(cloud, -1, args.bin));
        const auto &row = frames.back();
        std::cerr << "interface frame " << frames.size() << " n " << row.n << " thickness "
                  << row.thickness << "\n";
      },
      1);
  auto layers = openOut(args.outDir + "/interface-layers.csv");
  layers << "frame time_ps axis mid n cubic hexagonal interfacial liquid "
            "shell_cubic shell_hex shell_stack shell_other q6 q6bar\n";
  auto fronts = openOut(args.outDir + "/interface-fronts.csv");
  fronts << "frame time_ps n axis length defined thickness pos0 width0 pos1 width1\n";
  std::vector<double> time;
  std::vector<double> pos0;
  std::vector<double> pos1;
  std::vector<double> thick;
  for (std::size_t f = 0; f < frames.size(); ++f) {
    const auto &row = frames[f];
    const double t = static_cast<double>(f) * args.ps;
    for (const auto &bin : row.rows) {
      layers << (f + 1) << " " << t << " " << row.axis << " " << bin.mid << " " << bin.n << " "
             << bin.cubic << " " << bin.hexagonal << " " << bin.interfacial << " " << bin.liquid
             << " " << bin.shellCubic << " " << bin.shellHex << " " << bin.shellStack << " "
             << bin.shellOther << " " << bin.q6 << " " << bin.q6bar << "\n";
    }
    const double p0 = row.fronts.size() > 0 && row.fronts[0].defined ? row.fronts[0].position : 0;
    const double w0 = row.fronts.size() > 0 && row.fronts[0].defined ? row.fronts[0].width : 0;
    const double p1 = row.fronts.size() > 1 && row.fronts[1].defined ? row.fronts[1].position : 0;
    const double w1 = row.fronts.size() > 1 && row.fronts[1].defined ? row.fronts[1].width : 0;
    fronts << (f + 1) << " " << t << " " << row.n << " " << row.axis << " " << row.length << " "
           << (row.defined ? 1 : 0) << " " << row.thickness << " " << p0 << " " << w0 << " " << p1
           << " " << w1 << "\n";
    if (row.defined) {
      time.push_back(t);
      pos0.push_back(p0);
      pos1.push_back(p1);
      thick.push_back(row.thickness);
    }
  }
  const auto v0 = science::bootstrapSlope(time, pos0, 2000, 27011);
  const auto v1 = science::bootstrapSlope(time, pos1, 2000, 27012);
  const auto vt = science::bootstrapSlope(time, thick, 2000, 27013);
  auto velocity = openOut(args.outDir + "/interface-velocity.csv");
  velocity << "quantity slope_A_per_ps lo hi n_frames\n";
  velocity << "front0 " << v0.slope << " " << v0.lo << " " << v0.hi << " " << v0.n << "\n";
  velocity << "front1 " << v1.slope << " " << v1.lo << " " << v1.hi << " " << v1.n << "\n";
  velocity << "thickness " << vt.slope << " " << vt.lo << " " << vt.hi << " " << vt.n << "\n";
  // Half the thickness rate is the mean advance of the two fronts when
  // the ice slab thickens from both sides.
  velocity << "mean_front " << (0.5 * vt.slope) << " " << (0.5 * vt.lo) << " " << (0.5 * vt.hi)
           << " " << vt.n << "\n";

  if (!frames.empty()) {
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(args.dump, 1, cloud, args.type);
    auto atoms = openOut(args.outDir + "/interface-atoms.csv");
    writeAtoms(atoms, science::labelAtoms(cloud, science::kQ6BarIceCut), 1, cloud.box);
    const auto study = iface::studyDump(args.dump, 1, science::kCutoffA, args.bin, false);
    auto causes = openOut(args.outDir + "/interface-causes.csv");
    causes << "n disagree neighbour_test q3_shell open_edge mixed\n";
    causes << study.n << " " << study.nDisagree << " "
           << study.causeCount[static_cast<int>(iface::Cause::neighbour)] << " "
           << study.causeCount[static_cast<int>(iface::Cause::q3Shell)] << " "
           << study.causeCount[static_cast<int>(iface::Cause::openEdge)] << " "
           << study.causeCount[static_cast<int>(iface::Cause::mixed)] << "\n";
  }
  (void)steps;
  return frames.empty() ? 1 : 0;
}

void writeCritical(std::ostream &os, const char *name, const science::CriticalSize &size,
                   const science::NucleusFrame &at, const science::ClusterPart &part) {
  const double denom = static_cast<double>(part.cubic + part.hexagonal);
  const double cubicity = denom > 0.0 ? static_cast<double>(part.cubic) / denom : 0.0;
  const double stack =
      part.size > 0 ? static_cast<double>(part.shellStack) / static_cast<double>(part.size) : 0.0;
  os << name << " " << size.nStar << " " << size.sigma << " " << size.lo << " " << size.hi << " "
     << size.onset << " " << size.window << " " << (size.upperBound ? 1 : 0) << " "
     << (size.defined ? 1 : 0) << " " << at.frame << " " << part.size << " " << part.cubic << " "
     << part.hexagonal << " " << part.interfacial << " " << part.shellStack << " " << cubicity
     << " " << stack << "\n";
}

int runNucleation(const Args &args) {
  const auto steps = readTimesteps(args.dump);
  std::vector<science::NucleusFrame> rows;
  sinp::forEachLammpsFrame(
      args.dump, 1, 0, args.type,
      [&](int frame, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
        auto rec = science::nucleusOf(cloud, frame, science::kQ6BarIceCut);
        rows.push_back(rec);
        std::cerr << "nucleus frame " << frame << " ext " << rec.extended.size << " core "
                  << rec.core.size << " q6 " << rec.q6.size << " shell " << rec.shell.size
                  << "\n";
      },
      1);
  auto clusters = openOut(args.outDir + "/nucleation-clusters.csv");
  clusters << "frame timestep nop chill_cubic chill_hex chill_inter chill_liquid "
              "shell_cubic shell_hex shell_stack shell_other q6_ice q6bar_mean "
              "stack_as_cubic stack_as_hex stack_as_inter stack_as_liquid "
              "polytype_both polytype_agree "
              "n_core c_core h_core i_core ss_core "
              "n_ext c_ext h_ext i_ext ss_ext "
              "n_shell c_shell h_shell i_shell ss_shell "
              "n_q6 c_q6 h_q6 i_q6 ss_q6\n";
  auto sizeOf = [](const science::ClusterPart &part) { return static_cast<double>(part.size); };
  std::vector<double> ext;
  std::vector<double> core;
  std::vector<double> shell;
  std::vector<double> q6;
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const auto &r = rows[i];
    const long long step =
        i < steps.size() ? steps[i] : 0;
    clusters << r.frame << " " << step << " " << r.nop << " " << r.chillCubic << " " << r.chillHex
             << " " << r.chillInter << " " << r.chillLiquid << " " << r.shellCubic << " "
             << r.shellHex << " " << r.shellStack << " " << r.shellOther << " " << r.q6Ice << " "
             << r.q6barMean << " " << r.stackAsCubic << " " << r.stackAsHex << " "
             << r.stackAsInter << " " << r.stackAsLiquid << " " << r.polytypeBoth << " "
             << r.polytypeAgree << " " << r.core.size << " " << r.core.cubic << " "
             << r.core.hexagonal << " " << r.core.interfacial << " " << r.core.shellStack << " "
             << r.extended.size << " " << r.extended.cubic << " " << r.extended.hexagonal << " "
             << r.extended.interfacial << " " << r.extended.shellStack << " " << r.shell.size
             << " " << r.shell.cubic << " " << r.shell.hexagonal << " " << r.shell.interfacial
             << " " << r.shell.shellStack << " " << r.q6.size << " " << r.q6.cubic << " "
             << r.q6.hexagonal << " " << r.q6.interfacial << " " << r.q6.shellStack << "\n";
    ext.push_back(sizeOf(r.extended));
    core.push_back(sizeOf(r.core));
    shell.push_back(sizeOf(r.shell));
    q6.push_back(sizeOf(r.q6));
  }
  const auto cExt = science::criticalSize(ext, 2000, 27021);
  const auto cCore = science::criticalSize(core, 2000, 27022);
  const auto cShell = science::criticalSize(shell, 2000, 27023);
  const auto cQ6 = science::criticalSize(q6, 2000, 27024);
  auto summary = openOut(args.outDir + "/nucleation-summary.csv");
  summary << "definition n_star sigma lo hi onset window upper_bound defined frame "
             "size cubic hexagonal interfacial shell_stack cubicity stack_fraction\n";
  auto at = [&](const science::CriticalSize &size) -> const science::NucleusFrame & {
    const int index = size.onset >= 0 && size.onset < static_cast<int>(rows.size()) ? size.onset : 0;
    return rows[static_cast<std::size_t>(index)];
  };
  if (!rows.empty()) {
    writeCritical(summary, "chill_extended", cExt, at(cExt), at(cExt).extended);
    writeCritical(summary, "chill_core", cCore, at(cCore), at(cCore).core);
    writeCritical(summary, "shell", cShell, at(cShell), at(cShell).shell);
    writeCritical(summary, "q6bar", cQ6, at(cQ6), at(cQ6).q6);
    const int onset = cExt.onset >= 0 ? cExt.onset + 1 : 1;
    molSys::PointCloud<molSys::Point<double>, double> cloud;
    cloud = sinp::readLammpsTrjO(args.dump, onset, cloud, args.type);
    const std::vector<double> box = cloud.box;
    auto atoms = openOut(args.outDir + "/nucleus-atoms.csv");
    writeAtoms(atoms, science::labelAtoms(std::move(cloud), science::kQ6BarIceCut), onset, box);
    const auto study = iface::studyDump(args.dump, onset, science::kCutoffA, 1.0, false);
    auto causes = openOut(args.outDir + "/nucleation-causes.csv");
    causes << "frame n disagree neighbour_test q3_shell open_edge mixed\n";
    causes << onset << " " << study.n << " " << study.nDisagree << " "
           << study.causeCount[static_cast<int>(iface::Cause::neighbour)] << " "
           << study.causeCount[static_cast<int>(iface::Cause::q3Shell)] << " "
           << study.causeCount[static_cast<int>(iface::Cause::openEdge)] << " "
           << study.causeCount[static_cast<int>(iface::Cause::mixed)] << "\n";
  }
  return rows.empty() ? 1 : 0;
}

} // namespace

int main(int argc, char **argv) {
  Args args;
  try {
    if (!parseArgs(argc, argv, args)) {
      return argc <= 1 ? 2 : 0;
    }
  } catch (const std::exception &err) {
    std::cerr << err.what() << "\n";
    return 2;
  }
#ifdef SEAMS_HAS_OPENMP
  if (args.threads > 0) {
    omp_set_num_threads(args.threads);
  }
#endif
  std::filesystem::create_directories(args.outDir);
  try {
    if (args.cmd == "interface") {
      return runInterface(args);
    }
    if (args.cmd == "nucleation") {
      return runNucleation(args);
    }
  } catch (const std::exception &err) {
    std::cerr << err.what() << "\n";
    return 1;
  }
  usage();
  return 2;
}
