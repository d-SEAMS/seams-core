// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
// C++ driver for the same step the OpenMM reporter takes.
// It links libOpenMM when the headers are on the include path.
// The Python reporter is the path exercised in the e2e script.
// This file is the C++ plugin-shaped entry: one context, two steps.

#include "OpenMM.h"
#include "seams.h"

#include <cstdio>
#include <vector>

int main() {
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  std::vector<OpenMM::Vec3> positions;
  for (const auto &shift : fcc) {
    for (const auto &site : basis) {
      positions.emplace_back((shift[0] + site[0]) * a * 0.1, (shift[1] + site[1]) * a * 0.1,
                             (shift[2] + site[2]) * a * 0.1);
    }
  }
  const int n = static_cast<int>(positions.size());
  OpenMM::System system;
  for (int i = 0; i < n; ++i) system.addParticle(16.0);
  const double box = a * 0.1;
  system.setDefaultPeriodicBoxVectors(OpenMM::Vec3(box, 0, 0), OpenMM::Vec3(0, box, 0),
                                      OpenMM::Vec3(0, 0, box));
  OpenMM::CustomExternalForce *force = new OpenMM::CustomExternalForce("0");
  for (int i = 0; i < n; ++i) force->addParticle(i);
  system.addForce(force);
  OpenMM::VerletIntegrator integrator(0.001);
  OpenMM::Platform &platform = OpenMM::Platform::getPlatformByName("Reference");
  OpenMM::Context context(system, integrator, platform);
  context.setPositions(positions);
  int ice = 0;
  int calls = 0;
  for (int step = 0; step < 2; ++step) {
    integrator.step(1);
    OpenMM::State state = context.getState(OpenMM::State::Positions);
    const auto &pos = state.getPositions();
    std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
      xyz[static_cast<std::size_t>(3 * i)] = pos[i][0] * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 1)] = pos[i][1] * 10.0;
      xyz[static_cast<std::size_t>(3 * i + 2)] = pos[i][2] * 10.0;
    }
    const double box9[9] = {0, 0, 0, a, a, a, 0, 0, 0};
    std::vector<int> labels(static_cast<std::size_t>(n));
    int counts[9] = {};
    seams_system *sys = seams_system_create();
    seams_system_set_frame(sys, xyz.data(), n, box9, nullptr, nullptr);
    seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels.data(), counts);
    seams_system_destroy(sys);
    ice += counts[0];
    calls += 1;
  }
  if (calls != 2 || ice != calls * n) {
    std::fprintf(stderr, "openmm c++ ice %d calls %d\n", ice, calls);
    return 1;
  }
  std::printf("openmm-cpp calls %d ice %d\n", calls, ice);
  return 0;
}
