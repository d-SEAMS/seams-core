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
    std::vector<int> labels(static_cast<std::size_t>(n));
    int counts[9] = {};
    const auto t0 = std::chrono::steady_clock::now();
    seams_system *sys = seams_system_create();
    const auto set = seams_system_set_frame(sys, xyz.data(), n, box9, nullptr, nullptr);
    const auto code =
        seams_chill_plus(sys, cutoffA_, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels.data(), counts);
    seams_system_destroy(sys);
    const auto t1 = std::chrono::steady_clock::now();
    seconds_ += std::chrono::duration<double>(t1 - t0).count();
    frames_ += 1;
    if (set == SEAMS_OK && code == SEAMS_OK) {
      ice_ += counts[SEAMS_LABEL_CUBIC];
      atoms_ += n;
    }
  }

  void finishAnalysis(int) override {}

  void writeOutput() override {
    const double each = frames_ > 0 ? seconds_ / frames_ : 0.0;
    std::printf("gromacs frames %d ice %d atoms %d overhead_s %.6e\n", frames_, ice_, atoms_, each);
    if (frames_ < 1 || ice_ != atoms_) {
      std::fprintf(stderr, "gromacs analysis did not label every atom cubic\n");
      std::exit(1);
    }
  }

 private:
  double cutoffA_ = 3.5;
  int frames_ = 0;
  int ice_ = 0;
  int atoms_ = 0;
  double seconds_ = 0.0;
};

}  // namespace

int main(int argc, char *argv[]) {
  return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<SeamsIce>(argc, argv);
}
