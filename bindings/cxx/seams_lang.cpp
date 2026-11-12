// Copyright (c) d-SEAMS developers
#include "seams_lang.h"

#include <bop.hpp>
#include <chill_plus_scheme.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <order_parameter.hpp>
#include <rdf.hpp>
#include <rdf2d.hpp>
#include <ring.hpp>
#include <seams_c_api.h>
#include <seams_input.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>
#include <topo_two_dim.hpp>
#include <tum_offload.hpp>
#include <water_select.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstring>
#include <exception>
#include <fstream>
#include <generic.hpp>
#include <memory>
#include <set>
#include <string>
#include <vector>

#ifndef SEAMS_LANG_VERSION
#define SEAMS_LANG_VERSION "0"
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

thread_local std::string g_error;

void set_error(const std::string &message) { g_error = message; }

// Called from a handler. A standard exception keeps its text.
// Any other type is reported as unknown and does not leave the handler.
void remember_exception() {
  try {
    throw;
  } catch (const std::exception &ex) {
    set_error(ex.what());
  } catch (...) {
    set_error("unknown error");
  }
}

void clear_error() { g_error.clear(); }

// A gro frame is a title, an atom count, that many atom lines, and a box.
// A missing file or a short frame is -1.
int countGroFrames(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    return -1;
  }
  int frames = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (!std::getline(in, line)) {
      return -1;
    }
    const char *first = line.data();
    const char *last = first + line.size();
    while (first < last && (*first == ' ' || *first == '\t')) {
      ++first;
    }
    int n = 0;
    const auto parsed = std::from_chars(first, last, n);
    if (parsed.ec != std::errc{} || n < 0) {
      return -1;
    }
    for (int i = 0; i < n; ++i) {
      if (!std::getline(in, line)) {
        return -1;
      }
    }
    if (!std::getline(in, line)) {
      return -1;
    }
    frames += 1;
  }
  return frames > 0 ? frames : -1;
}

// MODEL records are one frame each. A file with atoms and no MODEL is
// one frame. A missing file is -1.
int countPdbFrames(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    return -1;
  }
  int models = 0;
  int atoms = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("MODEL", 0) == 0) {
      models += 1;
    } else if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) {
      atoms += 1;
    }
  }
  if (models > 0) {
    return models;
  }
  return atoms > 0 ? 1 : -1;
}

// One xyz frame is a count, a comment, and that many atom lines.
// A missing file, a non-numeric count, or a short frame is -1.
int countXyzFrames(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    return -1;
  }
  int frames = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    int n = 0;
    const char *first = line.data();
    const char *last = first + line.size();
    const auto parsed = std::from_chars(first, last, n);
    if (parsed.ec != std::errc{} || parsed.ptr != last || n < 0) {
      return -1;
    }
    if (!std::getline(in, line)) {
      return -1;
    }
    for (int i = 0; i < n; ++i) {
      if (!std::getline(in, line)) {
        return -1;
      }
    }
    frames += 1;
  }
  return frames > 0 ? frames : -1;
}

std::string extension_of(const std::string &path) {
  const auto dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return {};
  }
  return path.substr(dot + 1);
}

Cloud load_path(const std::string &path, int frame, int typeI) {
  // A missing CON path is an error. readCon clears the cloud and returns
  // it, which would look like an empty frame.
  if (extension_of(path) == "con" && !gen::file_exists(path)) {
    throw std::runtime_error("cannot read " + path);
  }
  return sinp::openFrame(path, frame, typeI);
}

int resolve_type(const Cloud &cloud, int requested) {
  if (requested > 0) {
    return requested;
  }
  if (cloud.nop == 0) {
    return 1;
  }
  return cloud.pts[0].type;
}

bool keep_atom(int atom_type, int requested) {
  return requested <= 0 || atom_type == requested;
}

void write_ice_labels(seams_frame *frame, int type, int *labels, int counts9[9]);

std::vector<std::vector<int>> bond_list(Cloud &cloud, double cutoff, int typ,
                                        const std::string &graph, int k) {
  const double cand = cutoff + 1.5;
  if (graph == "cutoff") {
    return nneigh::neighListO(cutoff, cloud, typ);
  }
  if (graph == "knn-union" || graph == "union") {
    return nneigh::kNearestNeighbourList(cloud, k, cand, typ, false);
  }
  if (graph == "mutual" || graph == "knn" || graph == "knn-mutual" ||
      graph.empty()) {
    return nneigh::kNearestNeighbourList(cloud, k, cand, typ, true);
  }
  throw std::invalid_argument(
      "unknown bond graph '" + graph +
      "' (use mutual, cutoff, knn-union, or seeded)");
}

std::vector<std::vector<int>>
six_of(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }
  return six;
}

void fill_box(const Cloud &cloud, double box9[9]) {
  for (int i = 0; i < 9; ++i) {
    box9[i] = 0.0;
  }
  box9[0] = cloud.boxLow[0];
  box9[1] = cloud.boxLow[1];
  box9[2] = cloud.boxLow[2];
  box9[3] = cloud.box[0];
  box9[4] = cloud.box[1];
  box9[5] = cloud.box[2];
  if (cloud.box.size() >= 6) {
    box9[6] = cloud.box[3];
    box9[7] = cloud.box[4];
    box9[8] = cloud.box[5];
  }
}

void apply_box(Cloud &cloud, const double *box9) {
  cloud.boxLow = {box9[0], box9[1], box9[2]};
  const bool sheared = box9[6] != 0.0 || box9[7] != 0.0 || box9[8] != 0.0;
  if (sheared) {
    cloud.box = {box9[3], box9[4], box9[5], box9[6], box9[7], box9[8]};
  } else {
    cloud.box = {box9[3], box9[4], box9[5]};
  }
}

} // namespace

struct seams_frame {
  Cloud cloud;
  std::vector<double> xyz;
  std::vector<int> types;
  std::vector<int> ids;
  std::vector<int> labels;
  std::vector<double> ql;
  std::vector<double> qlbar;
  std::vector<int> per_atom;
  std::vector<double> rdf_r;
  std::vector<double> rdf_g;
  std::vector<int> ring_counts;
  std::vector<std::vector<int>> wide_list;
  std::vector<double> wide_xyz;
  double wide_box[9] = {};
  double wide_cutoff = -1.0;
  double wide_skin = -1.0;
  int wide_type = 0;
  std::string wide_graph;
  bool wide_ready = false;
  // Empty, or one entry per kept water oxygen. A list shorter than nop
  // means later analyses run on that subset and leave ions unclassified.
  std::vector<int> oxygen_index;

  void sync_arrays() {
    const int n = cloud.nop;
    xyz.resize(static_cast<std::size_t>(std::max(n, 0)) * 3);
    types.resize(static_cast<std::size_t>(std::max(n, 0)));
    ids.resize(static_cast<std::size_t>(std::max(n, 0)));
    for (int i = 0; i < n; ++i) {
      const auto &point = cloud.pts[static_cast<std::size_t>(i)];
      xyz[static_cast<std::size_t>(3 * i)] = point.x;
      xyz[static_cast<std::size_t>(3 * i + 1)] = point.y;
      xyz[static_cast<std::size_t>(3 * i + 2)] = point.z;
      types[static_cast<std::size_t>(i)] = point.type;
      ids[static_cast<std::size_t>(i)] = point.atomID;
    }
  }
};

const char *seams_lang_version(void) { return SEAMS_LANG_VERSION; }

const char *seams_last_error(void) { return g_error.c_str(); }

int seams_path_nframes(const char *path) {
  clear_error();
  if (path == nullptr) {
    set_error("path is null");
    return -1;
  }
  try {
    if (!gen::file_exists(path)) {
      set_error("cannot read " + std::string(path));
      return -1;
    }
    const auto ext = extension_of(path);
    if (ext == "lammpstrj" || ext == "dump") {
      const int frames = sinp::nLammpsFrames(path);
      if (frames <= 0) {
        set_error("unreadable dump");
        return -1;
      }
      return frames;
    }
    // Chemfiles first. A miss falls through to the native counters.
    const int steps = sinp::chemfilesSteps(path);
    if (steps >= 0) {
      return steps;
    }
    if (ext == "gro") {
      const int frames = countGroFrames(path);
      if (frames < 0) {
        set_error("unreadable trajectory");
        return -1;
      }
      return frames;
    }
    if (ext == "pdb") {
      const int frames = countPdbFrames(path);
      if (frames < 0) {
        set_error("unreadable trajectory");
        return -1;
      }
      return frames;
    }
    if (ext == "dcd") {
      set_error("unreadable trajectory");
      return -1;
    }
    if (ext == "xyz") {
      std::ifstream in(path);
      if (!in) {
        set_error("cannot read " + std::string(path));
        return -1;
      }
      const int frames = countXyzFrames(path);
      if (frames > 0) {
        return frames;
      }
      // A readable file that is not a counted xyz block is still one
      // native frame. configuration.xyz is a count, a box, and coordinates.
      return 1;
    }
    return sinp::nLammpsFrames(path);
  } catch (...) {
    remember_exception();
    return -1;
  }
}

seams_frame *seams_frame_read(const char *path, int frame, int type_filter) {
  clear_error();
  if (path == nullptr) {
    set_error("path is null");
    return nullptr;
  }
  try {
    auto owned = std::make_unique<seams_frame>();
    seams_frame *out = owned.get();
    out->cloud = load_path(path, frame, type_filter);
    if (out->cloud.currentFrame <= 0) {
      out->cloud.currentFrame = frame > 0 ? frame : 1;
    }
    out->sync_arrays();
    return owned.release();
  } catch (...) {
    remember_exception();
    return nullptr;
  }
}

seams_frame *seams_frame_from_arrays(const double *xyz, int n,
                                    const double *box9, const int *types,
                                    const int *ids) {
  clear_error();
  if (xyz == nullptr || box9 == nullptr || n <= 0) {
    set_error("xyz, box, and a positive atom count are required");
    return nullptr;
  }
  try {
    auto owned = std::make_unique<seams_frame>();
    seams_frame *out = owned.get();
    Cloud &cloud = out->cloud;
    cloud.nop = n;
    cloud.currentFrame = 1;
    apply_box(cloud, box9);
    cloud.pts.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      auto &point = cloud.pts[static_cast<std::size_t>(i)];
      point.x = xyz[3 * i];
      point.y = xyz[3 * i + 1];
      point.z = xyz[3 * i + 2];
      point.type = types != nullptr ? types[i] : 1;
      point.atomID = ids != nullptr ? ids[i] : i + 1;
      point.molID = point.atomID;
      cloud.idIndexMap[point.atomID] = i;
    }
    out->sync_arrays();
    return owned.release();
  } catch (...) {
    remember_exception();
    return nullptr;
  }
}

seams_frame *seams_frame_from_species(const double *xyz, int n,
                                     const double *box9, const int *types,
                                     const int *ids, const int *mols) {
  seams_frame *out = seams_frame_from_arrays(xyz, n, box9, types, ids);
  if (out == nullptr || mols == nullptr) {
    return out;
  }
  for (int i = 0; i < n; ++i) {
    out->cloud.pts[static_cast<std::size_t>(i)].molID = mols[i];
  }
  out->sync_arrays();
  return out;
}

int seams_frame_select(seams_frame *frame, int type, const int *mols, int nmol) {
  clear_error();
  if (frame == nullptr) {
    set_error("frame is null");
    return 1;
  }
  try {
    water::Request request;
    if (type > 0) {
      request.type = type;
    } else if (type < 0) {
      throw std::invalid_argument("atom type must be positive");
    }
    if (mols != nullptr && nmol > 0) {
      request.mol_mode = water::MolMode::Ids;
      request.molecules.assign(mols, mols + nmol);
    } else if (type <= 0) {
      request.mol_mode = water::MolMode::Water;
    }
    const auto selection = water::selectOxygens(frame->cloud, request);
    frame->oxygen_index = selection.index;
    return 0;
  } catch (const std::exception &ex) {
    set_error(ex.what());
    return 1;
  } catch (...) {
    set_error("unknown error");
    return 1;
  }
}

void seams_frame_free(seams_frame *frame) { delete frame; }

int seams_frame_nop(const seams_frame *frame) {
  return frame == nullptr ? 0 : frame->cloud.nop;
}

int seams_frame_index(const seams_frame *frame) {
  return frame == nullptr ? 0 : frame->cloud.currentFrame;
}

void seams_frame_box(const seams_frame *frame, double box9[9]) {
  if (frame == nullptr || box9 == nullptr) {
    return;
  }
  fill_box(frame->cloud, box9);
}

const double *seams_frame_positions(const seams_frame *frame) {
  if (frame == nullptr || frame->xyz.empty()) {
    return nullptr;
  }
  return frame->xyz.data();
}

const int *seams_frame_types(const seams_frame *frame) {
  if (frame == nullptr || frame->types.empty()) {
    return nullptr;
  }
  return frame->types.data();
}

const int *seams_frame_ids(const seams_frame *frame) {
  if (frame == nullptr || frame->ids.empty()) {
    return nullptr;
  }
  return frame->ids.data();
}

namespace {

bool drops_atoms(const seams_frame *frame) {
  return frame != nullptr && !frame->oxygen_index.empty() &&
         static_cast<int>(frame->oxygen_index.size()) < frame->cloud.nop;
}

seams_frame oxygen_subframe(const seams_frame *frame) {
  water::Selection selection;
  selection.index = frame->oxygen_index;
  selection.n_oxygen = static_cast<int>(selection.index.size());
  selection.keep.assign(static_cast<std::size_t>(std::max(frame->cloud.nop, 0)),
                        0);
  for (int index : selection.index) {
    if (index >= 0 && index < frame->cloud.nop) {
      selection.keep[static_cast<std::size_t>(index)] = 1;
    }
  }
  seams_frame sub;
  sub.cloud = water::subsetCloud(frame->cloud, selection);
  return sub;
}

int classify_chill(seams_frame *frame, double cutoff, int type,
                   const char *graph, bool plus, const char *rule, int *labels,
                   int counts9[9]) {
  clear_error();
  if (frame == nullptr) {
    set_error("frame is null");
    return 1;
  }
  if (drops_atoms(frame)) {
    seams_frame sub = oxygen_subframe(frame);
    std::vector<int> subLabels(
        static_cast<std::size_t>(std::max(sub.cloud.nop, 0)));
    int subCounts[9] = {};
    const int rc = classify_chill(&sub, cutoff, type, graph, plus, rule,
                                  sub.cloud.nop > 0 ? subLabels.data() : nullptr,
                                  subCounts);
    if (rc != 0) {
      return rc;
    }
    const int n = frame->cloud.nop;
    const int unclassified =
        static_cast<int>(molSys::atom_state_type::unclassified);
    frame->labels.assign(static_cast<std::size_t>(std::max(n, 0)), unclassified);
    for (std::size_t k = 0; k < frame->oxygen_index.size(); ++k) {
      const int index = frame->oxygen_index[k];
      if (index >= 0 && index < n) {
        frame->labels[static_cast<std::size_t>(index)] =
            subLabels[k];
      }
    }
    if (labels != nullptr && n > 0) {
      std::memcpy(labels, frame->labels.data(),
                  static_cast<std::size_t>(n) * sizeof(int));
    }
    if (counts9 != nullptr) {
      std::memcpy(counts9, subCounts, sizeof(subCounts));
    }
    return 0;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    const std::string which =
        (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
    if (which != "cutoff" && which != "mutual") {
      throw std::invalid_argument("unknown CHILL graph '" + which +
                                  "' (use cutoff or mutual)");
    }
    const std::string ruleText = rule == nullptr ? "" : rule;
    const bool published =
        ruleText.empty() || ruleText == "chill+" || ruleText == "published";
    if (plus && which == "cutoff" && published) {
      std::vector<int> subset;
      const std::vector<int> *sub = nullptr;
      if (typ > 0) {
        subset.reserve(static_cast<std::size_t>(frame->cloud.nop));
        for (int i = 0; i < frame->cloud.nop; ++i) {
          if (frame->cloud.pts[static_cast<std::size_t>(i)].type == typ) {
            subset.push_back(i);
          }
        }
        sub = &subset;
      }
      std::vector<int> offsets;
      std::vector<int> cols;
      if (nneigh::orthoFullCsr(frame->cloud, sub, cutoff, offsets, cols, false)) {
        chill::chillPlusFromCsr(frame->cloud, offsets, cols);
        write_ice_labels(frame, type, labels, counts9);
        return 0;
      }
    }
    auto nList = bond_list(frame->cloud, cutoff, typ, which, 4);
    if (plus) {
      chill::applyChillPlus(frame->cloud, nList, chill::chillPlusRule(ruleText),
                            false);
    } else {
      chill::getCorrel(frame->cloud, nList, false);
      chill::getIceTypeNoPrint(frame->cloud, nList, false);
    }
    const int n = frame->cloud.nop;
    frame->labels.assign(static_cast<std::size_t>(std::max(n, 0)), 6);
    int counts[9] = {};
    for (int i = 0; i < n; ++i) {
      const auto &point = frame->cloud.pts[static_cast<std::size_t>(i)];
      int label = static_cast<int>(point.iceType);
      if (!keep_atom(point.type, type)) {
        label = static_cast<int>(molSys::atom_state_type::unclassified);
      }
      frame->labels[static_cast<std::size_t>(i)] = label;
      if (keep_atom(point.type, type)) {
        counts[label] += 1;
      }
    }
    if (labels != nullptr) {
      if (n > 0) {
        std::memcpy(labels, frame->labels.data(),
                    static_cast<std::size_t>(n) * sizeof(int));
      }
    }
    if (counts9 != nullptr) {
      std::memcpy(counts9, counts, sizeof(counts));
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

std::vector<std::vector<int>> index_graph(seams_frame *frame, double cutoff,
                                          int typ, const std::string &graph,
                                          int k) {
  if (graph == "seeded") {
    throw std::invalid_argument("seeded is a cage graph, not a ring graph");
  }
  auto ids = bond_list(frame->cloud, cutoff, typ, graph, k);
  return nneigh::neighbourListByIndex(frame->cloud, ids);
}

// Pairs of a wide list whose current minimum-image distance is inside cutoff.
// The same strict test as neighListO: a pair on the cutoff is not a neighbour.
std::vector<std::vector<int>>
pairs_within(const Cloud &cloud, const std::vector<std::vector<int>> &wide,
             double cutoff) {
  const double cut2 = cutoff * cutoff;
  const gen::FracBox box = gen::makeFracBox(cloud);
  std::vector<std::vector<int>> out(wide.size());
  for (std::size_t row = 0; row < wide.size(); ++row) {
    if (wide[row].empty()) {
      continue;
    }
    auto &kept = out[row];
    kept.push_back(wide[row][0]);
    const auto self = cloud.idIndexMap.find(wide[row][0]);
    if (self == cloud.idIndexMap.end()) {
      continue;
    }
    const auto &here = cloud.pts[static_cast<std::size_t>(self->second)];
    for (std::size_t j = 1; j < wide[row].size(); ++j) {
      const auto other = cloud.idIndexMap.find(wide[row][j]);
      if (other == cloud.idIndexMap.end()) {
        continue;
      }
      const auto &there = cloud.pts[static_cast<std::size_t>(other->second)];
      const double dist2 = gen::fracDistSq(box, here.x, here.y, here.z, there.x,
                                           there.y, there.z);
      if (dist2 < cut2) {
        kept.push_back(wide[row][j]);
      }
    }
  }
  return out;
}

double max_drift(const Cloud &cloud, const std::vector<double> &previous) {
  const gen::FracBox box = gen::makeFracBox(cloud);
  double furthest = 0.0;
  const int n = cloud.nop;
  for (int i = 0; i < n; ++i) {
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    const std::size_t base = static_cast<std::size_t>(3 * i);
    const double dist2 =
        gen::fracDistSq(box, point.x, point.y, point.z, previous[base],
                        previous[base + 1], previous[base + 2]);
    if (dist2 > furthest) {
      furthest = dist2;
    }
  }
  return std::sqrt(furthest);
}

bool wide_reusable(const seams_frame *frame, double cutoff, double skin, int typ,
                   const std::string &graph, const double *box9) {
  if (frame == nullptr || !frame->wide_ready || !(skin > 0.0) ||
      frame->wide_graph != graph || frame->wide_type != typ ||
      frame->wide_cutoff != cutoff || frame->wide_skin != skin) {
    return false;
  }
  for (int k = 0; k < 9; ++k) {
    if (frame->wide_box[k] != box9[k]) {
      return false;
    }
  }
  if (frame->wide_xyz.size() !=
      static_cast<std::size_t>(std::max(frame->cloud.nop, 0)) * 3) {
    return false;
  }
  return max_drift(frame->cloud, frame->wide_xyz) <= 0.5 * skin;
}

void remember_wide(seams_frame *frame, std::vector<std::vector<int>> wide,
                   double cutoff, double skin, int typ, const std::string &graph,
                   const double *box9) {
  frame->wide_list = std::move(wide);
  const int n = frame->cloud.nop;
  frame->wide_xyz.resize(static_cast<std::size_t>(std::max(n, 0)) * 3);
  for (int i = 0; i < n; ++i) {
    const auto &point = frame->cloud.pts[static_cast<std::size_t>(i)];
    const std::size_t base = static_cast<std::size_t>(3 * i);
    frame->wide_xyz[base] = point.x;
    frame->wide_xyz[base + 1] = point.y;
    frame->wide_xyz[base + 2] = point.z;
  }
  std::memcpy(frame->wide_box, box9, sizeof(frame->wide_box));
  frame->wide_cutoff = cutoff;
  frame->wide_skin = skin;
  frame->wide_type = typ;
  frame->wide_graph = graph;
  frame->wide_ready = true;
}

void write_ice_labels(seams_frame *frame, int type, int *labels, int counts9[9]) {
  const int n = frame->cloud.nop;
  frame->labels.assign(static_cast<std::size_t>(std::max(n, 0)), 6);
  int counts[9] = {};
  for (int i = 0; i < n; ++i) {
    const auto &point = frame->cloud.pts[static_cast<std::size_t>(i)];
    int label = static_cast<int>(point.iceType);
    if (!keep_atom(point.type, type)) {
      label = static_cast<int>(molSys::atom_state_type::unclassified);
    }
    frame->labels[static_cast<std::size_t>(i)] = label;
    if (keep_atom(point.type, type)) {
      counts[label] += 1;
    }
  }
  if (labels != nullptr && n > 0) {
    std::memcpy(labels, frame->labels.data(),
                static_cast<std::size_t>(n) * sizeof(int));
  }
  if (counts9 != nullptr) {
    std::memcpy(counts9, counts, sizeof(counts));
  }
}

} // namespace

int seams_frame_set_positions(seams_frame *frame, const double *xyz, int n,
                             const double *box9) {
  clear_error();
  if (frame == nullptr || xyz == nullptr || box9 == nullptr || n <= 0) {
    set_error("positions need a frame, a box, and a positive atom count");
    return 1;
  }
  if (frame->cloud.nop != n ||
      frame->cloud.pts.size() != static_cast<std::size_t>(n)) {
    set_error("position update needs the same atom count");
    return 1;
  }
  try {
    apply_box(frame->cloud, box9);
    frame->xyz.resize(static_cast<std::size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
      auto &point = frame->cloud.pts[static_cast<std::size_t>(i)];
      const std::size_t base = static_cast<std::size_t>(3 * i);
      point.x = xyz[base];
      point.y = xyz[base + 1];
      point.z = xyz[base + 2];
      frame->xyz[base] = point.x;
      frame->xyz[base + 1] = point.y;
      frame->xyz[base + 2] = point.z;
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_chill_rings(seams_frame *frame, double cutoff, int type,
                           const char *graph, int max_ring, double skin,
                           int *labels, int counts9[9], int *ring_counts,
                           int ring_counts_len, seams_phase_times *times) {
  clear_error();
  if (frame == nullptr || labels == nullptr || counts9 == nullptr ||
      ring_counts == nullptr || max_ring < 3 ||
      ring_counts_len < max_ring + 1 || !(cutoff > 0.0)) {
    set_error("CHILL+ and rings need a frame, a cutoff, and count buffers");
    return 1;
  }
  try {
    using clock = std::chrono::steady_clock;
    const int typ = resolve_type(frame->cloud, type);
    const std::string which =
        (graph == nullptr || graph[0] == '\0') ? std::string("cutoff") : graph;
    double box9[9];
    fill_box(frame->cloud, box9);
    const auto t0 = clock::now();
    int reused = 0;
    std::vector<std::vector<int>> bonds;
    if (which == "cutoff" && skin > 0.0 &&
        wide_reusable(frame, cutoff, skin, typ, which, box9)) {
      bonds = pairs_within(frame->cloud, frame->wide_list, cutoff);
      reused = 1;
    } else if (which == "cutoff" && skin > 0.0) {
      // Leave the OpenMP team as it is. Changing its width sleeps the
      // workers and the ring search pays to create them again.
      nneigh::preferSerialCutoff(true);
      auto wide = nneigh::neighListO(cutoff + skin, frame->cloud, typ);
      nneigh::preferSerialCutoff(false);
      remember_wide(frame, wide, cutoff, skin, typ, which, box9);
      bonds = pairs_within(frame->cloud, frame->wide_list, cutoff);
    } else if (which == "cutoff") {
      frame->wide_ready = false;
      std::vector<int> subset;
      const std::vector<int> *sub = nullptr;
      if (typ > 0) {
        subset.reserve(static_cast<std::size_t>(frame->cloud.nop));
        for (int i = 0; i < frame->cloud.nop; ++i) {
          if (frame->cloud.pts[static_cast<std::size_t>(i)].type == typ) {
            subset.push_back(i);
          }
        }
        sub = &subset;
      }
      std::vector<int> offsets;
      std::vector<int> cols;
      if (!nneigh::orthoFullCsr(frame->cloud, sub, cutoff, offsets, cols, false)) {
        bonds = bond_list(frame->cloud, cutoff, typ, which, 4);
      } else {
        const auto t1 = clock::now();
        chill::chillPlusFromCsr(frame->cloud, offsets, cols);
        write_ice_labels(frame, type, labels, counts9);
        const auto t2 = clock::now();
        const auto rings = primitive::ringNetwork(offsets, cols, max_ring);
        std::fill(ring_counts, ring_counts + ring_counts_len, 0);
        for (const auto &ring : rings) {
          const int size = static_cast<int>(ring.size());
          if (size >= 0 && size < ring_counts_len) {
            ring_counts[size] += 1;
          }
        }
        const auto t3 = clock::now();
        if (times != nullptr) {
          const auto seconds = [](clock::duration span) {
            return std::chrono::duration<double>(span).count();
          };
          times->neighbours = seconds(t1 - t0);
          times->chill = seconds(t2 - t1);
          times->rings = seconds(t3 - t2);
          times->reused = 0;
        }
        return 0;
      }
    } else {
      frame->wide_ready = false;
      nneigh::preferSerialCutoff(true);
      bonds = bond_list(frame->cloud, cutoff, typ, which, 4);
      nneigh::preferSerialCutoff(false);
    }
    const auto t1 = clock::now();
    chill::applyChillPlus(frame->cloud, bonds, chill::chillPlusRule(""), false);
    write_ice_labels(frame, type, labels, counts9);
    const auto t2 = clock::now();
    const auto indexed = nneigh::neighbourListByIndex(frame->cloud, bonds);
    const auto rings = primitive::ringNetwork(indexed, max_ring);
    std::fill(ring_counts, ring_counts + ring_counts_len, 0);
    for (const auto &ring : rings) {
      const int size = static_cast<int>(ring.size());
      if (size >= 0 && size < ring_counts_len) {
        ring_counts[size] += 1;
      }
    }
    const auto t3 = clock::now();
    if (times != nullptr) {
      const auto seconds = [](clock::duration span) {
        return std::chrono::duration<double>(span).count();
      };
      times->neighbours = seconds(t1 - t0);
      times->chill = seconds(t2 - t1);
      times->rings = seconds(t3 - t2);
      times->reused = reused;
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_chill_plus(seams_frame *frame, double cutoff, int type,
                           const char *graph, int *labels, int counts9[9]) {
  return classify_chill(frame, cutoff, type, graph, true, nullptr, labels,
                        counts9);
}

int seams_frame_chill_plus_rule(seams_frame *frame, double cutoff, int type,
                                const char *graph, const char *rule,
                                int *labels, int counts9[9]) {
  return classify_chill(frame, cutoff, type, graph, true, rule, labels,
                        counts9);
}

int seams_frame_chill(seams_frame *frame, double cutoff, int type,
                      const char *graph, int *labels, int counts9[9]) {
  return classify_chill(frame, cutoff, type, graph, false, nullptr, labels,
                        counts9);
}

int seams_frame_rings(seams_frame *frame, double cutoff, int type,
                      const char *graph, int max_depth, int *counts,
                      int counts_len) {
  clear_error();
  if (frame == nullptr || counts == nullptr || max_depth < 3 ||
      counts_len < max_depth + 1) {
    set_error("rings need a frame, max_depth >= 3, and counts_len");
    return 1;
  }
  if (drops_atoms(frame)) {
    seams_frame sub = oxygen_subframe(frame);
    return seams_frame_rings(&sub, cutoff, type, graph, max_depth, counts,
                             counts_len);
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    const std::string which = graph == nullptr ? "mutual" : graph;
    if (which == "cutoff") {
      std::vector<int> subset;
      const std::vector<int> *sub = nullptr;
      if (typ > 0) {
        subset.reserve(static_cast<std::size_t>(frame->cloud.nop));
        for (int i = 0; i < frame->cloud.nop; ++i) {
          if (frame->cloud.pts[static_cast<std::size_t>(i)].type == typ) {
            subset.push_back(i);
          }
        }
        sub = &subset;
      }
      std::vector<int> offsets;
      std::vector<int> cols;
      if (nneigh::orthoFullCsr(frame->cloud, sub, cutoff, offsets, cols, false)) {
        const auto rings = primitive::ringNetwork(offsets, cols, max_depth);
        std::fill(counts, counts + counts_len, 0);
        for (const auto &ring : rings) {
          const int size = static_cast<int>(ring.size());
          if (size >= 0 && size < counts_len) {
            counts[size] += 1;
          }
        }
        return 0;
      }
    }
    auto idx = index_graph(frame, cutoff, typ, which, 4);
    const auto rings = primitive::ringNetwork(idx, max_depth);
    std::fill(counts, counts + counts_len, 0);
    for (const auto &ring : rings) {
      counts[static_cast<int>(ring.size())] += 1;
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_neighbours(seams_frame *frame, double cutoff, int type,
                            const char *graph) {
  clear_error();
  if (frame == nullptr || !(cutoff > 0.0)) {
    set_error("neighbours need a frame and a positive cutoff");
    return 1;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    const std::string which = graph == nullptr ? "cutoff" : graph;
    auto rows = bond_list(frame->cloud, cutoff, typ, which, 4);
    return rows.size() == static_cast<std::size_t>(std::max(frame->cloud.nop, 0))
               ? 0
               : 1;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_steinhardt(seams_frame *frame, double cutoff, int type,
                           int order_l, double *ql, double *qlbar) {
  return seams_frame_steinhardt_graph(frame, cutoff, type, order_l, "mutual", ql,
                                      qlbar);
}

int seams_frame_steinhardt_graph(seams_frame *frame, double cutoff, int type,
                                 int order_l, const char *graph, double *ql,
                                 double *qlbar) {
  clear_error();
  if (frame == nullptr) {
    set_error("frame is null");
    return 1;
  }
  if (drops_atoms(frame)) {
    seams_frame sub = oxygen_subframe(frame);
    std::vector<double> subQl(static_cast<std::size_t>(std::max(sub.cloud.nop, 0)));
    std::vector<double> subBar(subQl.size());
    const int rc = seams_frame_steinhardt_graph(
        &sub, cutoff, type, order_l, graph, sub.cloud.nop > 0 ? subQl.data() : nullptr,
        sub.cloud.nop > 0 ? subBar.data() : nullptr);
    if (rc != 0) {
      return rc;
    }
    const int n = frame->cloud.nop;
    frame->ql.assign(static_cast<std::size_t>(std::max(n, 0)), 0.0);
    frame->qlbar.assign(frame->ql.size(), 0.0);
    for (std::size_t k = 0; k < frame->oxygen_index.size(); ++k) {
      const int index = frame->oxygen_index[k];
      if (index >= 0 && index < n) {
        frame->ql[static_cast<std::size_t>(index)] = subQl[k];
        frame->qlbar[static_cast<std::size_t>(index)] = subBar[k];
      }
    }
    if (ql != nullptr && n > 0) {
      std::memcpy(ql, frame->ql.data(), static_cast<std::size_t>(n) * sizeof(double));
    }
    if (qlbar != nullptr && n > 0) {
      std::memcpy(qlbar, frame->qlbar.data(),
                  static_cast<std::size_t>(n) * sizeof(double));
    }
    return 0;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    const std::string which =
        (graph == nullptr || graph[0] == '\0') ? "mutual" : graph;
    chill::SteinhardtQl result;
    if (which == "cutoff") {
      result = chill::steinhardtQlCutoff(frame->cloud, cutoff, order_l, typ, true);
    } else if (which == "mutual" || which == "knn" || which == "knn-mutual") {
      const double cand = cutoff + 1.5;
      auto nList =
          nneigh::kNearestNeighbourList(frame->cloud, 4, cand, typ, true);
      result = chill::steinhardtQl(frame->cloud, nList, order_l);
    } else {
      throw std::invalid_argument("steinhardt graph is cutoff or mutual");
    }
    const int n = frame->cloud.nop;
    frame->ql = result.ql;
    frame->qlbar = result.qlBar;
    if (ql != nullptr && n > 0) {
      std::memcpy(ql, frame->ql.data(),
                  static_cast<std::size_t>(n) * sizeof(double));
    }
    if (qlbar != nullptr && n > 0) {
      std::memcpy(qlbar, frame->qlbar.data(),
                  static_cast<std::size_t>(n) * sizeof(double));
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_cages(seams_frame *frame, double cutoff, int type, int k,
                      const char *graph, int complete, int *per_atom,
                      int *n_hex, int *n_cubic, int *n_water, int *n_both) {
  clear_error();
  if (frame == nullptr) {
    set_error("frame is null");
    return 1;
  }
  if (k <= 0) {
    k = 4;
  }
  if (drops_atoms(frame)) {
    seams_frame sub = oxygen_subframe(frame);
    std::vector<int> subPer(static_cast<std::size_t>(std::max(sub.cloud.nop, 0)));
    int ih = 0;
    int ic = 0;
    int iw = 0;
    int ib = 0;
    const int rc = seams_frame_cages(
        &sub, cutoff, type, k, graph, complete,
        sub.cloud.nop > 0 ? subPer.data() : nullptr, &ih, &ic, &iw, &ib);
    if (rc != 0) {
      return rc;
    }
    const int n = frame->cloud.nop;
    frame->per_atom.assign(static_cast<std::size_t>(std::max(n, 0)), 0);
    for (std::size_t s = 0; s < frame->oxygen_index.size() && s < subPer.size();
         ++s) {
      const int index = frame->oxygen_index[s];
      if (index >= 0 && index < n) {
        frame->per_atom[static_cast<std::size_t>(index)] = subPer[s];
      }
    }
    if (per_atom != nullptr && n > 0) {
      std::memcpy(per_atom, frame->per_atom.data(),
                  static_cast<std::size_t>(n) * sizeof(int));
    }
    if (n_hex != nullptr) {
      *n_hex = ih;
    }
    if (n_cubic != nullptr) {
      *n_cubic = ic;
    }
    if (n_water != nullptr) {
      *n_water = iw;
    }
    if (n_both != nullptr) {
      *n_both = ib;
    }
    return 0;
  }
  try {
    const int n = frame->cloud.nop;
    const int typ = resolve_type(frame->cloud, type);
    const std::string which = graph == nullptr ? "seeded" : graph;
    std::vector<char> hc(static_cast<std::size_t>(std::max(n, 0)), 0);
    std::vector<char> ddc(static_cast<std::size_t>(std::max(n, 0)), 0);

    auto tally_flags = [&](const std::vector<char> &hc_flags,
                           const std::vector<char> &ddc_flags) {
      int ih = 0;
      int ic = 0;
      int water = 0;
      int both = 0;
      frame->per_atom.assign(static_cast<std::size_t>(n), 0);
      for (int i = 0; i < n; ++i) {
        const auto &point = frame->cloud.pts[static_cast<std::size_t>(i)];
        if (!keep_atom(point.type, type)) {
          continue;
        }
        const bool h = hc_flags[static_cast<std::size_t>(i)] != 0;
        const bool d = ddc_flags[static_cast<std::size_t>(i)] != 0;
        // Per-atom codes: 0 water, 1 hexagonal only, 2 double-diamond only,
        // 3 both. The four counters are disjoint and sum to the kept atoms.
        const int code = h && d ? 3 : h ? 1 : d ? 2 : 0;
        frame->per_atom[static_cast<std::size_t>(i)] = code;
        if (h && d) {
          ++both;
        } else if (h) {
          ++ih;
        } else if (d) {
          ++ic;
        } else {
          ++water;
        }
      }
      if (per_atom != nullptr && n > 0) {
        std::memcpy(per_atom, frame->per_atom.data(),
                    static_cast<std::size_t>(n) * sizeof(int));
      }
      if (n_hex != nullptr) {
        *n_hex = ih;
      }
      if (n_cubic != nullptr) {
        *n_cubic = ic;
      }
      if (n_water != nullptr) {
        *n_water = water;
      }
      if (n_both != nullptr) {
        *n_both = both;
      }
    };

    if (n == 0) {
      tally_flags(hc, ddc);
      return 0;
    }

    const double cand = cutoff + 1.5;
    if (which == "seeded") {
      // tum::cageCounts reads one graph and has no ring completion, so the
      // seeded assignment stays on the host whatever SEAMS_OFFLOAD says.
      // The same rule as the CLI (GEM-6).
      auto graphs = nneigh::kNearestNeighbourPair(frame->cloud, k, cand, typ);
      auto idx_strict = nneigh::neighbourListByIndex(frame->cloud, graphs.first);
      auto idx_union = nneigh::neighbourListByIndex(frame->cloud, graphs.second);
      auto six_strict = six_of(ring::contractibleRings(
          primitive::ringNetwork(idx_strict, 6), frame->cloud));
      auto six_union = six_of(ring::contractibleRings(
          primitive::ringNetwork(idx_union, 6), frame->cloud));
      const auto aff = ring::seededCageAffiliation(
          six_strict, idx_strict, six_union, idx_union, complete != 0);
      const int m = std::min(n, static_cast<int>(aff.hc.size()));
      for (int i = 0; i < m; ++i) {
        hc[static_cast<std::size_t>(i)] = aff.hc[static_cast<std::size_t>(i)];
        ddc[static_cast<std::size_t>(i)] = aff.ddc[static_cast<std::size_t>(i)];
      }
      tally_flags(hc, ddc);
      return 0;
    }

    auto ids = bond_list(frame->cloud, cutoff, typ, which, k);
    auto idx = nneigh::neighbourListByIndex(frame->cloud, ids);
    if (tum::preferOffload()) {
      const auto counts = tum::cageCounts(idx);
      const int m = std::min(n, static_cast<int>(counts.atomHc.size()));
      for (int i = 0; i < m; ++i) {
        hc[static_cast<std::size_t>(i)] =
            counts.atomHc[static_cast<std::size_t>(i)] != 0;
        ddc[static_cast<std::size_t>(i)] =
            counts.atomDdc[static_cast<std::size_t>(i)] != 0;
      }
    } else {
      auto six = six_of(ring::contractibleRings(
          primitive::ringNetwork(idx, 6), frame->cloud));
      const auto aff = ring::cageAffiliation(six, idx);
      for (std::size_t r = 0; r < six.size(); ++r) {
        for (const int atom : six[r]) {
          hc[static_cast<std::size_t>(atom)] =
              hc[static_cast<std::size_t>(atom)] || aff.hc[r];
          ddc[static_cast<std::size_t>(atom)] =
              ddc[static_cast<std::size_t>(atom)] || aff.ddc[r];
        }
      }
    }
    tally_flags(hc, ddc);
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_signature(seams_frame *frame, double cutoff, int type,
                          const char *graph, const char *spec, int *n_cages,
                          int *n_atoms) {
  return seams_frame_signature_edge(frame, cutoff, type, graph, spec, nullptr,
                                    n_cages, n_atoms);
}

int seams_frame_signature_edge(seams_frame *frame, double cutoff, int type,
                               const char *graph, const char *spec,
                               const char *edge, int *n_cages, int *n_atoms) {
  clear_error();
  if (frame == nullptr || spec == nullptr || n_cages == nullptr ||
      n_atoms == nullptr) {
    set_error("signature needs a frame, a census, and count pointers");
    return 1;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    const std::string which =
        (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
    auto idx = index_graph(frame, cutoff, typ, which, 4);
    const auto signature = cage::Signature::parse(spec);
    const int depth = std::max(signature.maxRingSize(), 3);
    const auto rings = ring::cageRings(primitive::ringNetwork(idx, depth),
                                       frame->cloud, cutoff);
    const auto window = cage::parseEdgeWindow(edge == nullptr ? "" : edge);
    std::vector<std::array<double, 3>> xyz;
    const std::vector<std::array<double, 3>> *xyzPtr = nullptr;
    const std::vector<double> *boxPtr = nullptr;
    if (window.enabled) {
      xyz.resize(static_cast<std::size_t>(frame->cloud.nop));
      for (int i = 0; i < frame->cloud.nop; ++i) {
        const auto &p = frame->cloud.pts[static_cast<std::size_t>(i)];
        xyz[static_cast<std::size_t>(i)] = {p.x, p.y, p.z};
      }
      xyzPtr = &xyz;
      boxPtr = &frame->cloud.box;
    }
    const auto foundOf = cage::findBySignatures(
        rings, idx, std::vector<cage::Signature>{signature}, {}, xyzPtr,
        boxPtr, window);
    const auto &found = foundOf.at(0);
    std::set<int> atoms;
    for (const auto &cage_i : found) {
      atoms.insert(cage_i.vertices.begin(), cage_i.vertices.end());
    }
    *n_cages = static_cast<int>(found.size());
    *n_atoms = static_cast<int>(atoms.size());
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_rdf_nbins(double rmax, int nbins) {
  if (nbins > 0) {
    return nbins;
  }
  return std::max(1, static_cast<int>(std::lround(rmax / 0.1)));
}

int seams_frame_rdf(const seams_frame *frame, double rmax, int nbins,
                    int type_i, int type_j, double *r, double *g, int *count,
                    double *volume) {
  clear_error();
  if (frame == nullptr || r == nullptr || g == nullptr || count == nullptr) {
    set_error("rdf needs a frame and r, g, count buffers");
    return 1;
  }
  if (nbins <= 0) {
    set_error("rdf needs a positive bin count");
    return 1;
  }
  try {
    const int bins = nbins;
    const int typ_i = resolve_type(frame->cloud, type_i);
    const int typ_j = type_j > 0 ? type_j : typ_i;
    const auto histogram =
        rdf::partialRdf(frame->cloud, typ_i, typ_j, rmax, bins);
    const int n = static_cast<int>(histogram.r.size());
    for (int i = 0; i < bins; ++i) {
      r[i] = 0.0;
      g[i] = 0.0;
      count[i] = 0;
    }
    for (int i = 0; i < n && i < bins; ++i) {
      r[i] = histogram.r[static_cast<std::size_t>(i)];
      g[i] = histogram.g[static_cast<std::size_t>(i)];
      count[i] = histogram.count[static_cast<std::size_t>(i)];
    }
    if (volume != nullptr) {
      *volume = histogram.volume;
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_rdf2d_nbin(double cutoff, double binwidth) {
  if (!(binwidth > 0.0)) {
    binwidth = 0.05;
  }
  int nbin = static_cast<int>(cutoff / binwidth);
  return nbin < 1 ? 1 : nbin;
}

int seams_frame_rdf2d(seams_frame *frame, double cutoff, double binwidth,
                      const char *output_dir, double *r, double *g) {
  clear_error();
  if (frame == nullptr || r == nullptr || g == nullptr) {
    set_error("rdf2d needs a frame and r, g buffers");
    return 1;
  }
  try {
    const int nbin = seams_frame_rdf2d_nbin(cutoff, binwidth);
    if (output_dir != nullptr && output_dir[0] != '\0') {
      std::string dir(output_dir);
      if (dir.back() != '/') {
        dir.push_back('/');
      }
      std::vector<double> values;
      const int frame_no = frame->cloud.currentFrame > 0
                               ? frame->cloud.currentFrame
                               : 1;
      frame->cloud.currentFrame = frame_no;
      rdf2::rdf2Danalysis_AA(dir, values, frame->cloud, cutoff, binwidth,
                             frame_no, frame_no);
      for (int i = 0; i < nbin; ++i) {
        r[i] = binwidth * (static_cast<double>(i) + 0.5);
        g[i] = i < static_cast<int>(values.size()) ? values[static_cast<std::size_t>(i)]
                                                   : 0.0;
      }
      return 0;
    }
    auto histogram = rdf2::sampleRDF_AA(frame->cloud, cutoff, binwidth, nbin);
    std::vector<double> values(static_cast<std::size_t>(nbin), 0.0);
    auto lengths = rdf2::getSystemLengths(frame->cloud);
    rdf2::normalizeRDF(frame->cloud.nop, values, histogram, binwidth, nbin,
                       lengths, 1);
    for (int i = 0; i < nbin; ++i) {
      r[i] = binwidth * (static_cast<double>(i) + 0.5);
      g[i] = values[static_cast<std::size_t>(i)];
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_prisms(seams_frame *frame, double cutoff, int type,
                       int max_depth, const char *output_dir,
                       seams_prism_row *rows, int *n_rows) {
  clear_error();
  if (frame == nullptr || max_depth < 3) {
    set_error("prisms need a frame and max_depth >= 3");
    return 1;
  }
  const int expect = max_depth - 2;
  if (rows == nullptr) {
    if (n_rows != nullptr) {
      *n_rows = expect;
    }
    return 0;
  }
  if (n_rows == nullptr || *n_rows < expect) {
    set_error("prism row buffer is short");
    return 1;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    auto ids = bond_list(frame->cloud, cutoff, typ, "mutual", 4);
    auto idx = nneigh::neighbourListByIndex(frame->cloud, ids);
    const auto rings = primitive::ringNetwork(idx, max_depth);
    if (output_dir != nullptr && output_dir[0] != '\0') {
      std::string dir(output_dir);
      if (dir.back() != '/') {
        dir.push_back('/');
      }
      int atom_id = 0;
      const int frame_no =
          frame->cloud.currentFrame > 0 ? frame->cloud.currentFrame : 1;
      ring::prismAnalysis(dir, rings, idx, frame->cloud, max_depth, atom_id,
                          frame_no, frame_no, false);
    }
    for (int size = 3; size <= max_depth; ++size) {
      auto one = ring::getSingleRingSize(rings, size);
      seams_prism_row &row = rows[size - 3];
      row.size = size;
      row.perfect = 0;
      row.imperfect = 0;
      if (one.empty()) {
        continue;
      }
      std::vector<ring::strucType> ring_type(one.size());
      std::vector<double> rmsd(static_cast<std::size_t>(frame->cloud.nop), -1.0);
      int perfect = 0;
      int imperfect = 0;
      ring::findPrisms(one, ring_type, perfect, imperfect, idx, frame->cloud,
                       rmsd, false);
      row.perfect = perfect;
      row.imperfect = imperfect;
    }
    *n_rows = expect;
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_monolayer(seams_frame *frame, double cutoff, int type,
                          int max_depth, double sheet_area,
                          const char *output_dir, int *counts, int counts_len) {
  clear_error();
  if (frame == nullptr || counts == nullptr || max_depth < 3 ||
      counts_len < max_depth + 1) {
    set_error("monolayer rings need a frame and a counts buffer");
    return 1;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    auto ids = bond_list(frame->cloud, cutoff, typ, "cutoff", 4);
    auto idx = nneigh::neighbourListByIndex(frame->cloud, ids);
    auto rings = primitive::ringNetwork(idx, max_depth);
    if (output_dir != nullptr && output_dir[0] != '\0') {
      std::string dir(output_dir);
      if (dir.back() != '/') {
        dir.push_back('/');
      }
      const int frame_no =
          frame->cloud.currentFrame > 0 ? frame->cloud.currentFrame : 1;
      ring::polygonRingAnalysis(dir, rings, idx, frame->cloud, max_depth,
                                sheet_area, frame_no);
    }
    std::fill(counts, counts + counts_len, 0);
    for (const auto &ring : rings) {
      counts[static_cast<int>(ring.size())] += 1;
    }
    return 0;
  } catch (...) {
    remember_exception();
    return 1;
  }
}

int seams_frame_symmetry(seams_frame *frame, double shell_rmsd, int *labels) {
  clear_error();
  if (frame == nullptr || labels == nullptr) {
    set_error("symmetry needs a frame and a label buffer");
    return 1;
  }
  const int n = frame->cloud.nop;
  if (n <= 0) {
    return 0;
  }
  std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
  for (int i = 0; i < n; ++i) {
    const auto &point = frame->cloud.pts[static_cast<std::size_t>(i)];
    xyz[static_cast<std::size_t>(3 * i)] = point.x;
    xyz[static_cast<std::size_t>(3 * i + 1)] = point.y;
    xyz[static_cast<std::size_t>(3 * i + 2)] = point.z;
  }
  double lattice[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  const auto &box = frame->cloud.box;
  if (box.size() >= 3) {
    lattice[0] = box[0];
    lattice[4] = box[1];
    lattice[8] = box[2];
  }
  if (box.size() >= 6) {
    lattice[3] = box[3];
    lattice[6] = box[4];
    lattice[7] = box[5];
  }
  const int code =
      seams_symmetry_classify_tol(xyz.data(), n, lattice, shell_rmsd, labels);
  if (code != 0) {
    set_error("symmetry classify failed");
  }
  return code;
}

int seams_frames_jump_tau90(const seams_frame *first, const seams_frame *second,
                            double dt, int oxygen_type, int hydrogen_type,
                            double *tau) {
  clear_error();
  if (first == nullptr || second == nullptr || tau == nullptr) {
    set_error("jump tau needs two frames and an output");
    return 1;
  }
  *tau = topoparam::jumpRotorTau90(first->cloud, second->cloud, dt, oxygen_type,
                                   hydrogen_type);
  return 0;
}
