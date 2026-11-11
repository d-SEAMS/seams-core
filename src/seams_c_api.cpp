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

#include <seams_c_api.h>

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>

#include <exception>
#include <string>
#include <vector>

int seams_chill_plus(const double *xyz, int n, const double *box, int *labels) {
  return seams_chill_plus_graph(xyz, n, box, 3.5, nullptr, labels);
}

namespace {

int chillLabels(molSys::PointCloud<molSys::Point<double>, double> &cloud,
                double cutoff, const std::string &name, const char *rule,
                int *labels) {
  const auto nList =
      name == "cutoff"
          ? nneigh::neighListO(cutoff, cloud, 1)
          : nneigh::kNearestNeighbourList(cloud, 4, cutoff + 1.5, 1, true);
  chill::applyChillPlus(cloud, nList,
                        chill::chillPlusRule(rule == nullptr ? "" : rule),
                        false);
  for (int i = 0; i < cloud.nop; i++) {
    labels[i] = static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
  }
  return 0;
}

molSys::PointCloud<molSys::Point<double>, double>
cloudFromXyz(const double *xyz, int n) {
  molSys::PointCloud<molSys::Point<double>, double> cloud;
  cloud.nop = n;
  cloud.currentFrame = 1;
  cloud.pts.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    auto &p = cloud.pts[static_cast<std::size_t>(i)];
    p.x = xyz[3 * i];
    p.y = xyz[3 * i + 1];
    p.z = xyz[3 * i + 2];
    p.type = 1;
    p.atomID = i + 1;
    p.molID = i + 1;
    cloud.idIndexMap[i + 1] = i;
  }
  return cloud;
}

} // namespace

int seams_chill_plus_graph(const double *xyz, int n, const double *box,
                           double cutoff, const char *graph, int *labels) {
  if (xyz == nullptr || box == nullptr || labels == nullptr || n <= 0 ||
      !(cutoff > 0.0)) {
    return 1;
  }
  const std::string name =
      (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
  if (name != "cutoff" && name != "mutual") {
    return 1;
  }
  if (!(box[3] > 0.0 && box[4] > 0.0 && box[5] > 0.0)) {
    return 1;
  }
  auto cloud = cloudFromXyz(xyz, n);
  cloud.box = {box[3], box[4], box[5]};
  cloud.boxLow = {box[0], box[1], box[2]};
  return chillLabels(cloud, cutoff, name, nullptr, labels);
}

int seams_chill_plus_rule(const double *xyz, int n, const double *box,
                          double cutoff, const char *graph, const char *rule,
                          int *labels) {
  if (xyz == nullptr || box == nullptr || labels == nullptr || n <= 0 ||
      !(cutoff > 0.0)) {
    return 1;
  }
  const std::string name =
      (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
  if (name != "cutoff" && name != "mutual") {
    return 1;
  }
  if (!(box[3] > 0.0 && box[4] > 0.0 && box[5] > 0.0)) {
    return 1;
  }
  try {
    auto cloud = cloudFromXyz(xyz, n);
    cloud.box = {box[3], box[4], box[5]};
    cloud.boxLow = {box[0], box[1], box[2]};
    return chillLabels(cloud, cutoff, name, rule, labels);
  } catch (const std::exception &) {
    return 1;
  }
}

int seams_chill_plus_graph_lammps(const double *xyz, int n,
                                  const double *box9, double cutoff,
                                  const char *graph, int *labels) {
  if (xyz == nullptr || box9 == nullptr || labels == nullptr || n <= 0 ||
      !(cutoff > 0.0)) {
    return 1;
  }
  const std::string name =
      (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
  if (name != "cutoff" && name != "mutual") {
    return 1;
  }
  auto cloud = cloudFromXyz(xyz, n);
  if (!nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, box9[0], box9[1],
                                    box9[2], box9[3], box9[4], box9[5],
                                    box9[6], box9[7], box9[8])) {
    return 1;
  }
  return chillLabels(cloud, cutoff, name, nullptr, labels);
}
