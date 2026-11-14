// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// Trajectory analysis tool. analyzeFrame is the per-frame hook.

#include "gromacs/fileio/trxio.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

#include "seams.h"
#include "../frame_census.hpp"
#include "../triclinic_box.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

class SeamsIce : public gmx::TrajectoryAnalysisModule {
 public:
  void initOptions(gmx::IOptionsContainer *options, gmx::TrajectoryAnalysisSettings *settings) override {
    static const char *const help[] = {
        "CHILL+ on a trajectory. Coordinates are converted from nm to Angstrom.",
        "analyzeFrame is the in-run hook used for every frame.",
        "-select picks the sites that are classified. A multi-site water",
        "model needs its oxygens only, for example -s topol.tpr -select 'name OW'.",
        "-dump appends one LAMMPS dump frame per analysed frame with the",
        "coordinates handed to d-SEAMS, in Angstrom, and the CHILL+ label."};
    settings->setHelpText(help);
    // Making molecules whole would move sites away from the coordinates mdrun wrote.
    settings->setRmPBC(false);
    options->addOption(
        gmx::DoubleOption("cutoff").store(&cutoffA_).defaultValue(3.5).description("Cutoff in Angstrom"));
    options->addOption(gmx::SelectionOption("select")
                           .store(&sel_)
                           .defaultSelectionText("all")
                           .onlyAtoms()
                           .description("Sites to classify"));
    options->addOption(
        gmx::StringOption("dump").store(&dump_).description("LAMMPS dump of the classified frames"));
  }

  void initAnalysis(const gmx::TrajectoryAnalysisSettings &, const gmx::TopologyInformation &) override {}

  void analyzeFrame(int, const t_trxframe &fr, t_pbc *, gmx::TrajectoryAnalysisModuleData *pdata) override {
    if (!fr.bX || !fr.bBox) {
      return;
    }
    const gmx::Selection &sel = pdata->parallelSelection(sel_);
    const int n = sel.posCount();
    if (n <= 0) {
      return;
    }
    std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
      const auto &x = sel.position(i).x();
      xyz[static_cast<std::size_t>(3 * i)] = static_cast<double>(x[0]) * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 1)] = static_cast<double>(x[1]) * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 2)] = static_cast<double>(x[2]) * 10.0;
    }
    double raw[3][3];
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) {
        raw[a][b] = static_cast<double>(fr.box[a][b]);
      }
    }
    double box9[9];
    gromacs_box_to_lammps(raw, box9, 10.0);
    // frame_census calls seams_chill_plus.
    const auto t0 = std::chrono::steady_clock::now();
    const FrameCensus census = frame_census(xyz.data(), n, box9, cutoffA_);
    const auto t1 = std::chrono::steady_clock::now();
    seconds_ += std::chrono::duration<double>(t1 - t0).count();
    frames_ += 1;
    if (!census.ok) {
      failed_ = true;
      return;
    }
    if (!dump_.empty()) {
      writeDump(fr.bStep ? static_cast<long long>(fr.step) : frames_ - 1, sel, xyz, box9,
                census.labels);
    }
    ice_ += census.cubic;
    atoms_ += n;
    hexagonal_ = census.hexagonal;
    six_ = census.six;
    q6_ = census.q6;
  }

  void finishAnalysis(int) override {}

  void writeOutput() override {
    const double each = frames_ > 0 ? seconds_ / frames_ : 0.0;
    std::printf("gromacs frames %d ice %d atoms %d overhead_s %.6e\n", frames_, ice_, atoms_, each);
    std::printf("gromacs hexagonal %d six-rings %d q6 %.6f\n", hexagonal_, six_, q6_);
    if (frames_ < 1 || failed_) {
      std::fprintf(stderr, "gromacs analysis did not classify a frame\n");
      std::exit(1);
    }
  }

 private:
  void writeDump(long long step, const gmx::Selection &sel, const std::vector<double> &xyz,
                 const double box9[9], const std::vector<int> &labels) {
    FILE *fp = std::fopen(dump_.c_str(), frames_ == 1 ? "w" : "a");
    if (fp == nullptr) {
      std::fprintf(stderr, "cannot write %s\n", dump_.c_str());
      std::exit(1);
    }
    const int n = sel.posCount();
    const double xy = box9[6], xz = box9[7], yz = box9[8];
    const double xlo = box9[0] + std::min({0.0, xy, xz, xy + xz});
    const double xhi = box9[0] + box9[3] + std::max({0.0, xy, xz, xy + xz});
    const double ylo = box9[1] + std::min(0.0, yz), yhi = box9[1] + box9[4] + std::max(0.0, yz);
    std::fprintf(fp, "ITEM: TIMESTEP\n%lld\nITEM: NUMBER OF ATOMS\n%d\n", step, n);
    std::fprintf(fp, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n%.17g %.17g %.17g\n%.17g %.17g %.17g\n"
                     "%.17g %.17g %.17g\n",
                 xlo, xhi, xy, ylo, yhi, xz, box9[2], box9[2] + box9[5], yz);
    std::fprintf(fp, "ITEM: ATOMS id type x y z label\n");
    for (int i = 0; i < n; ++i) {
      std::fprintf(fp, "%d 1 %.17g %.17g %.17g %d\n", sel.position(i).atomIndices()[0] + 1,
                   xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2], labels[i]);
    }
    std::fclose(fp);
  }

  gmx::Selection sel_;
  std::string dump_;
  double cutoffA_ = 3.5;
  int frames_ = 0;
  int ice_ = 0;
  int atoms_ = 0;
  int hexagonal_ = 0;
  int six_ = 0;
  double q6_ = 0.0;
  bool failed_ = false;
  double seconds_ = 0.0;
};

}  // namespace

int main(int argc, char *argv[]) {
  return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<SeamsIce>(argc, argv);
}
