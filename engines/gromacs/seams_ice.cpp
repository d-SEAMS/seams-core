// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// Trajectory analysis tool. analyzeFrame is the per-frame hook.

#include "gromacs/fileio/trxio.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

#include "seams.h"
#include "../frame_census.hpp"
#include "../triclinic_box.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

class SeamsIce : public gmx::TrajectoryAnalysisModule {
 public:
  void initOptions(gmx::IOptionsContainer *options, gmx::TrajectoryAnalysisSettings *settings) override {
    static const char *const help[] = {
        "CHILL+ on a trajectory. Coordinates are converted from nm to Angstrom.",
        "analyzeFrame is the in-run hook used for every frame."};
    settings->setHelpText(help);
    options->addOption(
        gmx::DoubleOption("cutoff").store(&cutoffA_).defaultValue(3.5).description("Cutoff in Angstrom"));
  }

  void initAnalysis(const gmx::TrajectoryAnalysisSettings &, const gmx::TopologyInformation &) override {}

  void analyzeFrame(int, const t_trxframe &fr, t_pbc *, gmx::TrajectoryAnalysisModuleData *) override {
    if (!fr.bX || !fr.bBox || fr.natoms <= 0) {
      return;
    }
    const int n = fr.natoms;
    std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
      xyz[static_cast<std::size_t>(3 * i)] = static_cast<double>(fr.x[i][0]) * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 1)] = static_cast<double>(fr.x[i][1]) * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 2)] = static_cast<double>(fr.x[i][2]) * 10.0;
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
