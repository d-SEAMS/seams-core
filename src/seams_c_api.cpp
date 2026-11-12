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
#include <water_select.hpp>

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

namespace {

thread_local std::string g_c_error;

void note_c_error(const std::string &message) { g_c_error = message; }

void clear_c_error() { g_c_error.clear(); }

bool sameType(const int *types, int n, int &type0) {
  if (types == nullptr || n <= 0) {
    type0 = 1;
    return true;
  }
  type0 = types[0];
  for (int i = 1; i < n; ++i) {
    if (types[i] != type0) {
      return false;
    }
  }
  return true;
}

molSys::PointCloud<molSys::Point<double>, double>
cloudFromSpecies(const double *xyz, int n, const int *types, const int *mols) {
  auto cloud = cloudFromXyz(xyz, n);
  for (int i = 0; i < n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.type = types != nullptr ? types[i] : 1;
    point.molID = mols != nullptr ? mols[i] : i + 1;
  }
  return cloud;
}

bool applyLammpsBox(molSys::PointCloud<molSys::Point<double>, double> &cloud,
                    const double *box9) {
  return nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, box9[0], box9[1],
                                      box9[2], box9[3], box9[4], box9[5],
                                      box9[6], box9[7], box9[8]);
}

water::Selection
selectSpecies(const molSys::PointCloud<molSys::Point<double>, double> &cloud,
              int oxygen_type, int mol_mode, const int *mol_ids, int n_mol_ids) {
  water::Request request;
  if (oxygen_type > 0) {
    request.type = oxygen_type;
  }
  if (mol_mode == 1) {
    request.mol_mode = water::MolMode::Water;
  } else if (mol_mode == 2) {
    if (mol_ids == nullptr || n_mol_ids <= 0) {
      throw std::invalid_argument("molecule selection is empty");
    }
    request.mol_mode = water::MolMode::Ids;
    request.molecules.assign(mol_ids, mol_ids + n_mol_ids);
  } else if (mol_mode != 0) {
    throw std::invalid_argument("molecule mode is 0, 1 (water) or 2 (ids)");
  }
  return water::selectOxygens(cloud, request);
}

} // namespace

const char *seams_c_error(void) { return g_c_error.c_str(); }

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

namespace {

int prepareSpecies(const double *xyz, int n, const double *box9,
                   const int *types, const int *mols, int oxygen_type,
                   int mol_mode, const int *mol_ids, int n_mol_ids,
                   double cutoff, const char *graph, std::string &name,
                   molSys::PointCloud<molSys::Point<double>, double> &full,
                   water::Selection &selection, bool &uniformPath) {
  clear_c_error();
  uniformPath = false;
  if (xyz == nullptr || box9 == nullptr || n <= 0 || !(cutoff > 0.0)) {
    note_c_error("xyz, box, and a positive atom count are required");
    return 1;
  }
  name = (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
  if (name != "cutoff" && name != "mutual" && !name.empty()) {
    note_c_error("graph is cutoff or mutual");
    return 1;
  }
  if (name.empty()) {
    name = "cutoff";
  }
  int type0 = 1;
  const bool uniform = sameType(types, n, type0);
  if (uniform && mol_mode == 0 && (oxygen_type == 0 || oxygen_type == type0)) {
    uniformPath = true;
    full = cloudFromXyz(xyz, n);
    if (!applyLammpsBox(full, box9)) {
      note_c_error("box edges do not make a cell");
      return 1;
    }
    return 0;
  }
  if (!uniform && oxygen_type <= 0 && mol_mode == 0) {
    note_c_error("mixed species need a water type or molecule");
    return 1;
  }
  try {
    full = cloudFromSpecies(xyz, n, types, mols);
    if (!applyLammpsBox(full, box9)) {
      note_c_error("box edges do not make a cell");
      return 1;
    }
    selection = selectSpecies(full, oxygen_type, mol_mode, mol_ids, n_mol_ids);
  } catch (const std::exception &ex) {
    note_c_error(ex.what());
    return 1;
  }
  return 0;
}

molSys::PointCloud<molSys::Point<double>, double>
typeOneCopy(const molSys::PointCloud<molSys::Point<double>, double> &cloud) {
  auto out = cloudFromXyz(nullptr, 0);
  out.nop = cloud.nop;
  out.currentFrame = 1;
  out.box = cloud.box;
  out.boxLow = cloud.boxLow;
  out.pts.resize(static_cast<std::size_t>(std::max(cloud.nop, 0)));
  for (int i = 0; i < cloud.nop; ++i) {
    auto &point = out.pts[static_cast<std::size_t>(i)];
    const auto &src = cloud.pts[static_cast<std::size_t>(i)];
    point.x = src.x;
    point.y = src.y;
    point.z = src.z;
    point.type = 1;
    point.atomID = i + 1;
    point.molID = i + 1;
    out.idIndexMap[point.atomID] = i;
  }
  return out;
}

} // namespace

int seams_chill_plus_species(const double *xyz, int n, const double *box9,
                             const int *types, const int *mols,
                             int oxygen_type, int mol_mode, const int *mol_ids,
                             int n_mol_ids, double cutoff, const char *graph,
                             int *labels, int *ion_ice, int *ion_liquid,
                             double *ion_front) {
  if (labels == nullptr) {
    note_c_error("labels are required");
    return 1;
  }
  std::string name;
  molSys::PointCloud<molSys::Point<double>, double> full;
  water::Selection selection;
  bool uniformPath = false;
  if (prepareSpecies(xyz, n, box9, types, mols, oxygen_type, mol_mode, mol_ids,
                     n_mol_ids, cutoff, graph, name, full, selection,
                     uniformPath) != 0) {
    return 1;
  }
  if (uniformPath) {
    const int rc = chillLabels(full, cutoff, name, nullptr, labels);
    if (ion_ice != nullptr) {
      *ion_ice = 0;
    }
    if (ion_liquid != nullptr) {
      *ion_liquid = 0;
    }
    if (ion_front != nullptr) {
      *ion_front = -1.0;
    }
    return rc;
  }
  auto oxygen = water::subsetCloud(full, selection);
  auto scored = typeOneCopy(oxygen);
  std::vector<int> subLabels(static_cast<std::size_t>(std::max(scored.nop, 0)),
                             6);
  if (scored.nop > 0 &&
      chillLabels(scored, cutoff, name, nullptr, subLabels.data()) != 0) {
    note_c_error("CHILL+ failed on the water oxygens");
    return 1;
  }
  for (int i = 0; i < n; ++i) {
    labels[i] = static_cast<int>(molSys::atom_state_type::unclassified);
  }
  for (int k = 0; k < selection.n_oxygen; ++k) {
    const int index = selection.index[static_cast<std::size_t>(k)];
    labels[index] = subLabels[static_cast<std::size_t>(k)];
    oxygen.pts[static_cast<std::size_t>(k)].iceType =
        static_cast<molSys::atom_state_type>(subLabels[static_cast<std::size_t>(k)]);
  }
  const auto ions = water::collectIons(full, selection);
  const auto report = water::ionFront(oxygen, ions);
  if (ion_ice != nullptr) {
    *ion_ice = report.ice;
  }
  if (ion_liquid != nullptr) {
    *ion_liquid = report.liquid;
  }
  if (ion_front != nullptr) {
    *ion_front = report.front;
  }
  return 0;
}

int seams_oxygen_cutoff_degree(const double *xyz, int n, const double *box9,
                               const int *types, const int *mols,
                               int oxygen_type, int mol_mode,
                               const int *mol_ids, int n_mol_ids, double cutoff,
                               int *degree) {
  if (degree == nullptr) {
    note_c_error("degree is required");
    return 1;
  }
  std::string name;
  molSys::PointCloud<molSys::Point<double>, double> full;
  water::Selection selection;
  bool uniformPath = false;
  if (prepareSpecies(xyz, n, box9, types, mols, oxygen_type, mol_mode, mol_ids,
                     n_mol_ids, cutoff, "cutoff", name, full, selection,
                     uniformPath) != 0) {
    return 1;
  }
  for (int i = 0; i < n; ++i) {
    degree[i] = -1;
  }
  molSys::PointCloud<molSys::Point<double>, double> scored;
  std::vector<int> back;
  if (uniformPath) {
    scored = full;
    back.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      back[static_cast<std::size_t>(i)] = i;
    }
  } else {
    auto oxygen = water::subsetCloud(full, selection);
    scored = typeOneCopy(oxygen);
    back = selection.index;
  }
  const auto rows = nneigh::neighListO(cutoff, scored, 1);
  if (static_cast<int>(rows.size()) < scored.nop) {
    note_c_error("neighbour list is empty");
    return 1;
  }
  for (int k = 0; k < scored.nop; ++k) {
    const auto &row = rows[static_cast<std::size_t>(k)];
    const int count = row.empty() ? 0 : static_cast<int>(row.size()) - 1;
    degree[back[static_cast<std::size_t>(k)]] = count;
  }
  return 0;
}
