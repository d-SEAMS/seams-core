//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

// Walk a LAMMPS dump and write one row per frame: primitive rings, both
// CHILL+ rules, Lechner-Dellago q6-bar, and diamond symmetry labels.
// A summary names the ring lead, the onset of the final cluster climb,
// and the shift between definitions.

#include <nucleation_pathway.hpp>
#include <seams_input.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Record = nucleation::FrameRecord;

struct Args {
  std::string dump;
  std::string out = "-";
  std::string summary;
  int type = 1;
  int first = 1;
  int last = 0;
  int threads = 0;
  double cutoff = nucleation::kCutoffA;
  int rings = nucleation::kRingDepth;
  double q6 = nucleation::kQ6BarIceCut;
};

void usage() {
  std::cerr
      << "nucleation-pathway DUMP [--out TSV] [--summary FILE]\n"
      << "  [--type 1] [--frame 1] [--last 0] [--threads 0]\n"
      << "  [--cutoff 3.5] [--rings 7] [--q6 0.36]\n";
}

bool parseArgs(int argc, char **argv, Args &args) {
  std::vector<std::string> pos;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    auto need = [&](const char *name) -> std::string {
      if (i + 1 >= argc) {
        throw std::runtime_error(std::string("missing value for ") + name);
      }
      return argv[++i];
    };
    if (key == "--out") {
      args.out = need("--out");
    } else if (key == "--summary") {
      args.summary = need("--summary");
    } else if (key == "--type") {
      args.type = std::stoi(need("--type"));
    } else if (key == "--frame") {
      args.first = std::stoi(need("--frame"));
    } else if (key == "--last") {
      args.last = std::stoi(need("--last"));
    } else if (key == "--threads") {
      args.threads = std::stoi(need("--threads"));
    } else if (key == "--cutoff") {
      args.cutoff = std::stod(need("--cutoff"));
    } else if (key == "--rings") {
      args.rings = std::stoi(need("--rings"));
    } else if (key == "--q6") {
      args.q6 = std::stod(need("--q6"));
    } else if (key == "--help") {
      usage();
      return false;
    } else if (key[0] == '-') {
      throw std::runtime_error("unknown option " + key);
    } else {
      pos.push_back(key);
    }
  }
  if (pos.size() != 1) {
    usage();
    return false;
  }
  args.dump = pos[0];
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

const char *kHeader =
    "frame\ttimestep\tnop\tring4\tring5\tring6\tring7\t"
    "chill_cubic\tchill_hex\tchill_inter\tchill_clath\tchill_interclath\tchill_water\t"
    "ovito_cubic\tovito_hex\tovito_inter\tovito_clath\tovito_interclath\tovito_water\t"
    "sym_cubic\tsym_hex\t"
    "q6bar_mean\tq6bar_p50\tq6bar_ice\t"
    "n_chill_core\tc_chill_core\th_chill_core\t"
    "n_chill_ext\tc_chill_ext\th_chill_ext\t"
    "n_ovito_core\tc_ovito_core\th_ovito_core\t"
    "n_ovito_ext\tc_ovito_ext\th_ovito_ext\t"
    "n_q6\tc_q6\th_q6\t"
    "n_sym\tc_sym\th_sym\t"
    "n_q6solid\tc_q6solid\th_q6solid\n";

void writeRow(std::ostream &os, const Record &r) {
  os << r.frame << '\t' << r.timestep << '\t' << r.nop << '\t' << r.ringCount[4]
     << '\t' << r.ringCount[5] << '\t' << r.ringCount[6] << '\t' << r.ringCount[7]
     << '\t' << r.chillCubic << '\t' << r.chillHex << '\t' << r.chillInter << '\t'
     << r.chillClath << '\t' << r.chillInterClath << '\t' << r.chillWater << '\t'
     << r.ovitoCubic << '\t' << r.ovitoHex << '\t' << r.ovitoInter << '\t'
     << r.ovitoClath << '\t' << r.ovitoInterClath << '\t' << r.ovitoWater << '\t'
     << r.symCubic << '\t' << r.symHex << '\t' << r.q6barMean << '\t' << r.q6barP50
     << '\t' << r.q6barIce << '\t' << r.nChillCore << '\t' << r.cChillCore << '\t'
     << r.hChillCore << '\t' << r.nChillExt << '\t' << r.cChillExt << '\t'
     << r.hChillExt << '\t' << r.nOvitoCore << '\t' << r.cOvitoCore << '\t'
     << r.hOvitoCore << '\t' << r.nOvitoExt << '\t' << r.cOvitoExt << '\t'
     << r.hOvitoExt << '\t' << r.nQ6 << '\t' << r.cQ6 << '\t' << r.hQ6 << '\t'
     << r.nSym << '\t' << r.cSym << '\t' << r.hSym << '\t' << r.nQ6Solid << '\t'
     << r.cQ6Solid << '\t' << r.hQ6Solid << '\n';
}

double fraction(int num, int den) {
  if (den <= 0) {
    return 0;
  }
  return static_cast<double>(num) / static_cast<double>(den);
}

void writeLead(std::ostream &os, const char *name, const nucleation::Lead &lead,
               double psPerFrame) {
  os << "ring_lead\t" << name << '\t';
  if (!lead.ok) {
    os << "ok\t0\n";
    return;
  }
  os << "ok\t1\thalf_frames\t" << lead.halfFrames << "\tascent_frames\t"
     << lead.ascentFrames << "\tcorr_lag\t" << lead.corrLag << "\tcorr\t" << lead.corr
     << "\thalf_ps\t" << lead.halfFrames * psPerFrame << "\tascent_ps\t"
     << lead.ascentFrames * psPerFrame << "\tcorr_ps\t" << lead.corrLag * psPerFrame
     << '\n';
}

void writeNucleus(std::ostream &os, const char *name, const std::vector<Record> &rows,
                  int (*sizeOf)(const Record &), int (*cubicOf)(const Record &),
                  int (*hexOf)(const Record &), double psPerFrame) {
  std::vector<double> series(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    series[i] = static_cast<double>(sizeOf(rows[i]));
  }
  const nucleation::Ascent climb = nucleation::finalAscent(series);
  os << "nucleus\t" << name << '\t';
  if (climb.index < 0) {
    os << "ok\t0\n";
    return;
  }
  const Record &at = rows[static_cast<std::size_t>(climb.index)];
  const Record &top = rows[static_cast<std::size_t>(climb.maxIndex)];
  const nucleation::HalfRise half = nucleation::halfRise(series);
  os << "ok\t1\tclimbed\t" << (climb.climbed ? 1 : 0) << "\tonset_frame\t" << at.frame
     << "\tonset_timestep\t" << at.timestep << "\tonset_size\t" << sizeOf(at)
     << "\tonset_stack\t" << nucleation::stackFraction(cubicOf(at), hexOf(at))
     << "\tmax_frame\t" << top.frame << "\tmax_timestep\t" << top.timestep
     << "\tmax_size\t" << sizeOf(top) << "\tmax_stack\t"
     << nucleation::stackFraction(cubicOf(top), hexOf(top)) << "\tonset_ps\t"
     << (at.frame - rows.front().frame) * psPerFrame;
  if (!half.ok) {
    os << "\thalf_ok\t0\n";
    return;
  }
  const Record &mid = rows[static_cast<std::size_t>(half.index)];
  os << "\thalf_ok\t1\thalf_frame\t" << mid.frame << "\thalf_timestep\t" << mid.timestep
     << "\thalf_size\t" << sizeOf(mid) << "\thalf_stack\t"
     << nucleation::stackFraction(cubicOf(mid), hexOf(mid)) << "\thalf_ps\t"
     << (mid.frame - rows.front().frame) * psPerFrame << '\n';
}

struct NucleusView {
  const char *name;
  int (*sizeOf)(const Record &);
  int (*cubicOf)(const Record &);
  int (*hexOf)(const Record &);
};

void writeSummary(std::ostream &os, const std::vector<Record> &rows, double q6Cut,
                  long long stepDelta, double stepFs) {
  const double psPerFrame = stepDelta * stepFs * 1e-3;
  os << "frames\t" << rows.size() << '\n';
  os << "molecules\t" << (rows.empty() ? 0 : rows.front().nop) << '\n';
  os << "cutoff_A\t" << nucleation::kCutoffA << '\n';
  os << "ring_depth\t" << nucleation::kRingDepth << '\n';
  os << "q6bar_cut\t" << q6Cut << '\n';
  os << "dump_step_delta\t" << stepDelta << '\n';
  os << "assumed_step_fs\t" << stepFs << '\n';
  os << "ps_per_frame\t" << psPerFrame << '\n';
  if (rows.empty()) {
    return;
  }
  std::vector<double> ringFrac(rows.size());
  std::vector<double> ringCount(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const double n5 = static_cast<double>(rows[i].ringCount[5]);
    const double n6 = static_cast<double>(rows[i].ringCount[6]);
    ringFrac[i] = (n5 + n6) > 0 ? n6 / (n5 + n6) : 0;
    ringCount[i] = n6;
  }
  auto leadAgainst = [&](const char *name, auto iceNum) {
    std::vector<double> ice(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
      ice[i] = fraction(iceNum(rows[i]), rows[i].nop);
    }
    writeLead(os, name, nucleation::ringLead(ringFrac, ice), psPerFrame);
  };
  leadAgainst("chill_ext", [](const Record &r) {
    return r.chillCubic + r.chillHex + r.chillInter;
  });
  leadAgainst("chill_core",
              [](const Record &r) { return r.chillCubic + r.chillHex; });
  leadAgainst("ovito_ext", [](const Record &r) {
    return r.ovitoCubic + r.ovitoHex + r.ovitoInter;
  });
  leadAgainst("ovito_core",
              [](const Record &r) { return r.ovitoCubic + r.ovitoHex; });
  leadAgainst("q6bar", [](const Record &r) { return r.q6barIce; });
  leadAgainst("symmetry", [](const Record &r) { return r.symCubic + r.symHex; });
  writeLead(os, "chill_ext_n6",
            nucleation::ringLead(ringCount, [&] {
              std::vector<double> ice(rows.size());
              for (std::size_t i = 0; i < rows.size(); ++i) {
                ice[i] = fraction(rows[i].chillCubic + rows[i].chillHex +
                                      rows[i].chillInter,
                                  rows[i].nop);
              }
              return ice;
            }()),
            psPerFrame);

  const NucleusView views[] = {
      {"chill_core", [](const Record &r) { return r.nChillCore; },
       [](const Record &r) { return r.cChillCore; },
       [](const Record &r) { return r.hChillCore; }},
      {"chill_ext", [](const Record &r) { return r.nChillExt; },
       [](const Record &r) { return r.cChillExt; },
       [](const Record &r) { return r.hChillExt; }},
      {"ovito_core", [](const Record &r) { return r.nOvitoCore; },
       [](const Record &r) { return r.cOvitoCore; },
       [](const Record &r) { return r.hOvitoCore; }},
      {"ovito_ext", [](const Record &r) { return r.nOvitoExt; },
       [](const Record &r) { return r.cOvitoExt; },
       [](const Record &r) { return r.hOvitoExt; }},
      {"q6bar", [](const Record &r) { return r.nQ6; },
       [](const Record &r) { return r.cQ6; },
       [](const Record &r) { return r.hQ6; }},
      {"symmetry", [](const Record &r) { return r.nSym; },
       [](const Record &r) { return r.cSym; },
       [](const Record &r) { return r.hSym; }},
  };
  struct Saved {
    const char *name;
    int frame;
    int size;
    bool climbed;
  };
  std::vector<Saved> saved;
  for (const auto &view : views) {
    writeNucleus(os, view.name, rows, view.sizeOf, view.cubicOf, view.hexOf,
                 psPerFrame);
    std::vector<double> series(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) {
      series[i] = static_cast<double>(view.sizeOf(rows[i]));
    }
    const auto half = nucleation::halfRise(series);
    Saved one;
    one.name = view.name;
    one.climbed = half.ok;
    one.frame = half.ok ? rows[static_cast<std::size_t>(half.index)].frame : -1;
    one.size = half.ok ? view.sizeOf(rows[static_cast<std::size_t>(half.index)]) : 0;
    saved.push_back(one);
  }
  for (std::size_t a = 0; a < saved.size(); ++a) {
    for (std::size_t b = a + 1; b < saved.size(); ++b) {
      if (!saved[a].climbed || !saved[b].climbed) {
        os << "shift\t" << saved[a].name << '\t' << saved[b].name << "\tok\t0\n";
        continue;
      }
      const auto shift = nucleation::definitionShift(
          saved[a].size, saved[a].frame, saved[b].size, saved[b].frame);
      os << "shift\t" << saved[a].name << '\t' << saved[b].name << "\tok\t1\td_size\t"
         << shift.dSize << "\td_frame\t" << shift.dFrame << "\td_ps\t"
         << shift.dFrame * psPerFrame << '\n';
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  Args args;
  try {
    if (!parseArgs(argc, argv, args)) {
      return argc == 1 ? 2 : 0;
    }
  } catch (const std::exception &err) {
    std::cerr << err.what() << '\n';
    return 2;
  }
#ifdef SEAMS_HAS_OPENMP
  if (args.threads > 0) {
    omp_set_num_threads(args.threads);
  }
#endif
  const auto steps = readTimesteps(args.dump);
  std::vector<Record> rows;
  sinp::forEachLammpsFrame(
      args.dump, args.first, args.last, args.type,
      [&](int frame, molSys::PointCloud<molSys::Point<double>, double> &cloud) {
        auto rec = nucleation::analyzeFrame(cloud, args.cutoff, args.rings, args.type,
                                             args.q6, frame);
        if (frame >= 1 && frame <= static_cast<int>(steps.size())) {
          rec.timestep = steps[static_cast<std::size_t>(frame - 1)];
        }
        rows.push_back(rec);
        std::cerr << "frame " << frame << " nop " << rec.nop << " chill "
                  << rec.chillCubic + rec.chillHex << " ovito "
                  << rec.ovitoCubic + rec.ovitoHex << " q6 " << rec.q6barIce
                  << " sym " << rec.symCubic + rec.symHex << " n6 " << rec.ringCount[6]
                  << '\n';
      },
      1);
  std::ofstream file;
  std::ostream *os = &std::cout;
  if (args.out != "-") {
    file.open(args.out);
    if (!file) {
      std::cerr << "cannot write " << args.out << '\n';
      return 1;
    }
    os = &file;
  }
  *os << kHeader;
  for (const auto &row : rows) {
    writeRow(*os, row);
  }
  long long delta = 0;
  if (steps.size() >= 2) {
    delta = steps[1] - steps[0];
    for (std::size_t i = 2; i < steps.size(); ++i) {
      if (steps[i] - steps[i - 1] != delta) {
        delta = 0;
        break;
      }
    }
  }
  if (!args.summary.empty()) {
    std::ofstream summary(args.summary);
    if (!summary) {
      std::cerr << "cannot write " << args.summary << '\n';
      return 1;
    }
    writeSummary(summary, rows, args.q6, delta, 10.0);
  }
  return 0;
}
