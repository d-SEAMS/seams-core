//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
//
// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the MIT License as published by
// the Open Source Initiative.
//
// A copy of the MIT License is included in the LICENSE file of this repository.
// You should have received a copy of the MIT License along with this program.
// If not, see <https://opensource.org/licenses/MIT>.
//-----------------------------------------------------------------------------------

#include "seams.h"

#include <bond.hpp>
#include <bop.hpp>
#include <brine_front.hpp>
#include <cluster.hpp>
#include <density.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <chill_plus_scheme.hpp>
#include <edge_set.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <local_symmetry.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <nucleation_pathway.hpp>
#include <rdf.hpp>
#include <topo_fingerprint.hpp>
#include <seams_c_api.h>
#include <seams_input.hpp>
#ifdef SEAMS_HAS_READCON
#include <readcon-core.hpp>
#endif
#include <site.hpp>
#include <thread_team.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>
#include <tum_offload.hpp>

#include <array>
#include <algorithm>
#include <charconv>
#include <climits>
#include <cmath>
#include <cstring>
#include <exception>
#include <fstream>
#include <memory>
#include <new>
#include <set>
#include <string>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

struct seams_object {
  std::string error;
};

thread_local std::string g_error;
thread_local seams_object *g_current = nullptr;
int g_threads = 0;

struct handle_guard {
  seams_object *previous;
  explicit handle_guard(seams_object *object) : previous(g_current) {
    g_current = object;
    if (object != nullptr) {
      object->error.clear();
    }
  }
  ~handle_guard() { g_current = previous; }
  handle_guard(const handle_guard &) = delete;
  handle_guard &operator=(const handle_guard &) = delete;
};

void note_error() {
  if (g_current == nullptr) {
    return;
  }
  try {
    g_current->error = g_error;
  } catch (const std::bad_alloc &) {
  }
}

#ifdef SEAMS_HAS_OPENMP
int saved_omp_threads() {
  static int saved = -1;
  if (saved < 0) {
    saved = omp_get_max_threads();
    if (saved < 1) {
      saved = 1;
    }
  }
  return saved;
}
#endif

void clear_error() { g_error.clear(); }

seams_status fail(seams_status code, const char *message) {
  g_error = message == nullptr ? "" : message;
  note_error();
  return code;
}

void remember() {
  try {
    throw;
  } catch (const std::bad_alloc &) {
    g_error = "out of memory";
  } catch (const std::exception &ex) {
    g_error = ex.what();
  } catch (...) {
    g_error = "unknown error";
  }
  note_error();
}

void use_threads() {
#ifdef SEAMS_HAS_OPENMP
  if (g_threads > 0) {
    omp_set_num_threads(g_threads);
  }
#endif
  seams::adoptThreads();
}

bool finite_box(const double *box9) {
  for (int i = 0; i < 9; ++i) {
    if (!std::isfinite(box9[i])) {
      return false;
    }
  }
  return true;
}

Cloud unified_copy(const Cloud &cloud, int type, int &typ) {
  Cloud work = cloud;
  if (type > 0) {
    typ = type;
    return work;
  }
  if (work.nop <= 0) {
    typ = 1;
    return work;
  }
  const int first = work.pts[0].type;
  bool mixed = first <= 0;
  for (int i = 1; i < work.nop && !mixed; ++i) {
    if (work.pts[static_cast<std::size_t>(i)].type != first) {
      mixed = true;
    }
  }
  if (mixed) {
    for (auto &point : work.pts) {
      point.type = 1;
    }
    typ = 1;
    return work;
  }
  typ = first > 0 ? first : 1;
  if (typ != first) {
    for (auto &point : work.pts) {
      point.type = typ;
    }
  }
  return work;
}

std::vector<std::vector<int>> bond_rows(const Cloud &cloud, double cutoff, int typ,
                                        int graph) {
  if (graph == SEAMS_GRAPH_MUTUAL) {
    return nneigh::kNearestNeighbourList(cloud, 4, cutoff + 1.5, typ, true);
  }
  if (graph != SEAMS_GRAPH_CUTOFF) {
    throw std::invalid_argument("graph is cutoff or mutual");
  }
  return nneigh::neighListO(cutoff, cloud, typ);
}

} // namespace

struct seams_system : seams_object {
  Cloud cloud;
  bool ready = false;
};

struct seams_neighbours : seams_object {
  std::vector<int> offsets;
  std::vector<int> cols;
  int n = 0;
};

struct seams_rings : seams_object {
  std::vector<std::vector<int>> rings;
  std::vector<int> symmetry;
  int n_atoms = 0;
};

struct seams_cage_set : seams_object {
  std::vector<int> per_atom;
  int n_atoms = 0;
  int n_hex = 0;
  int n_cubic = 0;
  int n_water = 0;
  int n_both = 0;
};

struct seams_brine : seams_object {
  site::FrontSnapshot snapshot;
};

struct seams_topo : seams_object {
  std::string frame_key;
  std::string method;
  std::vector<std::string> atom_keys;
  std::vector<int> ring_census;
  int hops = 0;
};

namespace {

std::string extension_of(const std::string &path) {
  const auto dot = path.find_last_of('.');
  if (dot == std::string::npos || dot + 1 >= path.size()) {
    return {};
  }
  return path.substr(dot + 1);
}

int count_gro_frames(const std::string &path) {
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

int count_pdb_frames(const std::string &path) {
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

int count_xyz_frames(const std::string &path) {
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

#ifdef SEAMS_HAS_READCON
int count_con_frames(const std::string &path) {
  int frames = 0;
  readcon::ConFrameIterator cursor(path);
  for (auto &&frame : cursor) {
    (void)frame;
    frames += 1;
  }
  return frames;
}
#endif

std::vector<std::vector<int>> six_of(const std::vector<std::vector<int>> &rings) {
  std::vector<std::vector<int>> six;
  for (const auto &ring : rings) {
    if (ring.size() == 6) {
      six.push_back(ring);
    }
  }
  return six;
}

bool buffer_fits(size_t len, int count) {
  return count >= 0 && len >= static_cast<std::size_t>(count);
}

void write_ice_labels(const Cloud &work, int type, int *labels, int counts9[9]) {
  const int n = work.nop;
  int counts[9] = {};
  for (int i = 0; i < n; ++i) {
    const auto &point = work.pts[static_cast<std::size_t>(i)];
    int label = static_cast<int>(point.iceType);
    if (type > 0 && point.type != type) {
      label = SEAMS_LABEL_UNCLASSIFIED;
    }
    if (label < 0 || label > 8) {
      label = SEAMS_LABEL_UNCLASSIFIED;
    }
    labels[i] = label;
    counts[label] += 1;
  }
  if (counts9 != nullptr) {
    std::memcpy(counts9, counts, sizeof(counts));
  }
}

double median_of(std::vector<double> values) {
  if (values.empty()) {
    return 0.0;
  }
  const std::size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid),
                   values.end());
  if (values.size() % 2 == 1) {
    return values[mid];
  }
  const double upper = values[mid];
  std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid - 1),
                   values.end());
  return 0.5 * (values[mid - 1] + upper);
}

void lattice_from_cloud(const Cloud &cloud, double lattice[9]) {
  lattice[0] = 1.0;
  lattice[1] = 0.0;
  lattice[2] = 0.0;
  lattice[3] = 0.0;
  lattice[4] = 1.0;
  lattice[5] = 0.0;
  lattice[6] = 0.0;
  lattice[7] = 0.0;
  lattice[8] = 1.0;
  const auto &box = cloud.box;
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
}

void tally_cages(seams_cage_set *out, const Cloud &work, int type, const std::vector<char> &hc,
                 const std::vector<char> &ddc) {
  const int n = work.nop;
  out->n_atoms = n;
  out->per_atom.assign(static_cast<std::size_t>(std::max(n, 0)), SEAMS_CAGE_WATER);
  out->n_hex = 0;
  out->n_cubic = 0;
  out->n_water = 0;
  out->n_both = 0;
  for (int i = 0; i < n; ++i) {
    const auto &point = work.pts[static_cast<std::size_t>(i)];
    if (type > 0 && point.type != type) {
      continue;
    }
    const bool hexagonal = static_cast<std::size_t>(i) < hc.size() && hc[static_cast<std::size_t>(i)] != 0;
    const bool cubic = static_cast<std::size_t>(i) < ddc.size() && ddc[static_cast<std::size_t>(i)] != 0;
    int code = SEAMS_CAGE_WATER;
    if (hexagonal && cubic) {
      code = SEAMS_CAGE_BOTH;
      out->n_both += 1;
    } else if (hexagonal) {
      code = SEAMS_CAGE_HEX;
      out->n_hex += 1;
    } else if (cubic) {
      code = SEAMS_CAGE_CUBIC;
      out->n_cubic += 1;
    } else {
      out->n_water += 1;
    }
    out->per_atom[static_cast<std::size_t>(i)] = code;
  }
}

void flags_from_affiliation(int n, const std::vector<std::vector<int>> &six,
                            const ring::CageAffiliation &aff, std::vector<char> &hc,
                            std::vector<char> &ddc) {
  for (std::size_t ring_i = 0; ring_i < six.size(); ++ring_i) {
    const bool hexagonal = ring_i < aff.hc.size() && aff.hc[ring_i];
    const bool cubic = ring_i < aff.ddc.size() && aff.ddc[ring_i];
    for (int atom : six[ring_i]) {
      if (atom < 0 || atom >= n) {
        continue;
      }
      if (hexagonal) {
        hc[static_cast<std::size_t>(atom)] = 1;
      }
      if (cubic) {
        ddc[static_cast<std::size_t>(atom)] = 1;
      }
    }
  }
}

} // namespace

extern "C" {

int seams_abi(void) { return SEAMS_ABI; }
int seams_version(void) { return SEAMS_VERSION_NUMBER; }
int seams_version_major(void) { return SEAMS_VERSION_MAJOR; }
int seams_version_minor(void) { return SEAMS_VERSION_MINOR; }
int seams_version_patch(void) { return SEAMS_VERSION_PATCH; }
const char *seams_version_string(void) { return SEAMS_VERSION_STRING; }

const char *seams_status_string(seams_status status) {
  switch (status) {
  case SEAMS_OK:
    return "ok";
  case SEAMS_ERR_NULL:
    return "null pointer";
  case SEAMS_ERR_RANGE:
    return "value out of range";
  case SEAMS_ERR_BOX:
    return "box is not a positive cell";
  case SEAMS_ERR_ARG:
    return "bad argument";
  case SEAMS_ERR_STATE:
    return "system has no frame";
  case SEAMS_ERR_NOMEM:
    return "out of memory";
  case SEAMS_ERR_INTERNAL:
    return "internal error";
  }
  return "unknown status";
}

const char *seams_last_error(void) { return g_error.c_str(); }

const char *seams_handle_error(const void *handle) {
  if (handle == nullptr) {
    return g_error.c_str();
  }
  try {
    return static_cast<const seams_object *>(handle)->error.c_str();
  } catch (...) {
    return "";
  }
}

const char *seams_label_name(int label) {
  switch (label) {
  case SEAMS_LABEL_CUBIC:
    return "cubic";
  case SEAMS_LABEL_HEXAGONAL:
    return "hexagonal";
  case SEAMS_LABEL_WATER:
    return "water";
  case SEAMS_LABEL_INTERFACIAL:
    return "interfacial";
  case SEAMS_LABEL_CLATHRATE:
    return "clathrate";
  case SEAMS_LABEL_INTERCLATHRATE:
    return "interClathrate";
  case SEAMS_LABEL_UNCLASSIFIED:
    return "unclassified";
  case SEAMS_LABEL_RECUBIC:
    return "reCubic";
  case SEAMS_LABEL_REHEX:
    return "reHex";
  default:
    return "unknown";
  }
}

const char *seams_ring_symmetry_name(int symmetry) {
  switch (symmetry) {
  case SEAMS_RING_UNCLASSIFIED:
    return "unclassified";
  case SEAMS_RING_DDC:
    return "ddc";
  case SEAMS_RING_HC_BASAL:
    return "hc-basal";
  case SEAMS_RING_HC_PRISMATIC:
    return "hc-prismatic";
  case SEAMS_RING_BOTH_BASAL:
    return "both-basal";
  case SEAMS_RING_BOTH_PRISMATIC:
    return "both-prismatic";
  case SEAMS_RING_PRISM:
    return "prism";
  case SEAMS_RING_DEFORMED_PRISM:
    return "deformed-prism";
  case SEAMS_RING_MIXED_PRISM:
    return "mixed-prism";
  case SEAMS_RING_PENTAGON:
    return "pentagon";
  case SEAMS_RING_HEXAGON:
    return "hexagon";
  case SEAMS_RING_OTHER:
    return "other";
  default:
    return "unknown";
  }
}

seams_status seams_set_num_threads(int n) {
  try {
  if (n < 0) {
    return fail(SEAMS_ERR_RANGE, "thread count is negative");
  }
  g_threads = n;
#ifdef SEAMS_HAS_OPENMP
  const int saved = saved_omp_threads();
  omp_set_num_threads(n > 0 ? n : saved);
#endif
#ifdef SEAMS_HAS_LINKCELL
  int width = n;
  if (width <= 0) {
#ifdef SEAMS_HAS_OPENMP
    width = omp_get_max_threads();
#else
    width = 1;
#endif
  }
  lc_set_threads(width);
#else
  (void)n;
#endif
  clear_error();
  return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

int seams_get_num_threads(void) {
  if (g_threads > 0) {
    return g_threads;
  }
#ifdef SEAMS_HAS_OPENMP
  const int n = omp_get_max_threads();
  return n > 0 ? n : 1;
#else
  return 1;
#endif
}

seams_system *seams_system_create(void) {
  try {
    clear_error();
    return new seams_system();
  } catch (const std::bad_alloc &) {
    g_error = "out of memory";
    return nullptr;
  }
}

void seams_system_destroy(seams_system *sys) {
  try {
    delete sys;
  } catch (...) {
  }
}

seams_status seams_system_set_atoms(seams_system *sys, const double *xyz, size_t n,
                                    const double *box9, size_t box_len, const int *types,
                                    const int *ids) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || xyz == nullptr || box9 == nullptr) {
    return fail(SEAMS_ERR_NULL, "system, coordinates, or box is null");
  }
  if (box_len != 9) {
    return fail(SEAMS_ERR_RANGE, "box length is 9");
  }
  if (n == 0) {
    return fail(SEAMS_ERR_RANGE, "frame has no atoms");
  }
  if (n > static_cast<std::size_t>(INT_MAX)) {
    return fail(SEAMS_ERR_RANGE, "atom count does not fit in the frame");
  }
  const int count = static_cast<int>(n);
  if (!finite_box(box9)) {
    return fail(SEAMS_ERR_BOX, "box contains a non-finite number");
  }
  try {
    Cloud cloud;
    if (!nneigh::setBoxFromLammpsEdges(cloud.box, cloud.boxLow, box9[0], box9[1],
                                      box9[2], box9[3], box9[4], box9[5], box9[6],
                                      box9[7], box9[8])) {
      return fail(SEAMS_ERR_BOX, "box edge lengths must be positive");
    }
    cloud.nop = count;
    cloud.currentFrame = 1;
    cloud.pts.resize(n);
    for (int i = 0; i < count; ++i) {
      const std::size_t base = static_cast<std::size_t>(i) * 3;
      const double x = xyz[base];
      const double y = xyz[base + 1];
      const double z = xyz[base + 2];
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
        return fail(SEAMS_ERR_ARG, "a coordinate is not finite");
      }
      auto &point = cloud.pts[static_cast<std::size_t>(i)];
      point.x = x;
      point.y = y;
      point.z = z;
      point.type = types == nullptr ? 1 : types[i];
      point.atomID = ids == nullptr ? i + 1 : ids[i];
      point.molID = point.atomID;
      point.iceType = molSys::atom_state_type::unclassified;
      if (!cloud.idIndexMap.emplace(point.atomID, i).second) {
        return fail(SEAMS_ERR_ARG, "atom ids are not unique");
      }
    }
    sys->cloud = std::move(cloud);
    sys->ready = true;
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return g_error == "out of memory" ? SEAMS_ERR_NOMEM : SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_system_set_frame(seams_system *sys, const double *xyz, int n,
                                    const double *box9, const int *types,
                                    const int *ids) {
  handle_guard guard(sys);
  clear_error();
  if (n < 0) {
    return fail(SEAMS_ERR_RANGE, "frame has no atoms");
  }
  return seams_system_set_atoms(sys, xyz, static_cast<std::size_t>(n), box9, 9, types,
                                ids);
}

seams_status seams_system_set_ortho(seams_system *sys, const double *xyz, int n,
                                    double ox, double oy, double oz, double lx,
                                    double ly, double lz, const int *types,
                                    const int *ids) {
  const double box9[9] = {ox, oy, oz, lx, ly, lz, 0.0, 0.0, 0.0};
  return seams_system_set_frame(sys, xyz, n, box9, types, ids);
}

int seams_system_count(const seams_system *sys) {
  if (sys == nullptr || !sys->ready) {
    return 0;
  }
  return sys->cloud.nop;
}

seams_status seams_neighbours_build(seams_system *sys, double cutoff, int type,
                                    int graph, seams_neighbours **out) {
  handle_guard guard(sys);
  clear_error();
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "neighbours output is null");
  }
  *out = nullptr;
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  if (graph != SEAMS_GRAPH_CUTOFF && graph != SEAMS_GRAPH_MUTUAL) {
    return fail(SEAMS_ERR_ARG, "graph is cutoff or mutual");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    auto rows = bond_rows(work, cutoff, typ, graph);
    auto idx = nneigh::neighbourListByIndex(work, rows);
    auto nb = std::make_unique<seams_neighbours>();
    const int n = work.nop;
    if (static_cast<int>(idx.size()) < n) {
      idx.resize(static_cast<std::size_t>(n));
    }
    nb->n = n;
    nb->offsets.assign(static_cast<std::size_t>(n + 1), 0);
    for (int i = 0; i < n; ++i) {
      const auto &row = idx[static_cast<std::size_t>(i)];
      int degree = 0;
      for (std::size_t k = 0; k < row.size(); ++k) {
        if (row[k] != i) {
          ++degree;
        }
      }
      nb->offsets[static_cast<std::size_t>(i + 1)] =
          nb->offsets[static_cast<std::size_t>(i)] + degree;
    }
    nb->cols.resize(static_cast<std::size_t>(nb->offsets[static_cast<std::size_t>(n)]));
    for (int i = 0; i < n; ++i) {
      int cursor = nb->offsets[static_cast<std::size_t>(i)];
      const auto &row = idx[static_cast<std::size_t>(i)];
      for (int id : row) {
        if (id == i) {
          continue;
        }
        nb->cols[static_cast<std::size_t>(cursor++)] = id;
      }
    }
    *out = nb.release();
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_neighbours_destroy(seams_neighbours *nb) {
  try {
    delete nb;
  } catch (...) {
  }
}

int seams_neighbours_count(const seams_neighbours *nb) {
  return nb == nullptr ? 0 : nb->n;
}

seams_status seams_neighbours_csr(const seams_neighbours *nb, const int **offsets,
                                  int *n_offsets, const int **cols, int *nnz) {
  handle_guard guard(const_cast<seams_neighbours *>(nb));
  clear_error();
  if (nb == nullptr || offsets == nullptr || n_offsets == nullptr || cols == nullptr ||
      nnz == nullptr) {
    return fail(SEAMS_ERR_NULL, "csr output is null");
  }
  *offsets = nb->offsets.data();
  *n_offsets = static_cast<int>(nb->offsets.size());
  *cols = nb->cols.empty() ? nullptr : nb->cols.data();
  *nnz = static_cast<int>(nb->cols.size());
  return SEAMS_OK;
}

seams_status seams_chill_plus(seams_system *sys, double cutoff, int type, int graph,
                              int rule, int *labels, int counts9[9]) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || labels == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or labels is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  if (rule != SEAMS_CHILL_PUBLISHED && rule != SEAMS_CHILL_OVITO) {
    return fail(SEAMS_ERR_ARG, "rule is published or ovito");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    if (rule == SEAMS_CHILL_PUBLISHED && graph == SEAMS_GRAPH_CUTOFF) {
      std::vector<int> subset;
      const std::vector<int> *sub = nullptr;
      if (type > 0) {
        subset.reserve(static_cast<std::size_t>(work.nop));
        for (int i = 0; i < work.nop; ++i) {
          if (work.pts[static_cast<std::size_t>(i)].type == typ) {
            subset.push_back(i);
          }
        }
        sub = &subset;
      }
      std::vector<int> offsets;
      std::vector<int> cols;
      if (nneigh::orthoFullCsr(work, sub, cutoff, offsets, cols, false)) {
        chill::chillPlusFromCsr(work, offsets, cols);
        const int n = work.nop;
        int counts[9] = {};
        for (int i = 0; i < n; ++i) {
          const auto &point = work.pts[static_cast<std::size_t>(i)];
          int label = static_cast<int>(point.iceType);
          if (type > 0 && point.type != type) {
            label = SEAMS_LABEL_UNCLASSIFIED;
          }
          if (label < 0 || label > 8) {
            label = SEAMS_LABEL_UNCLASSIFIED;
          }
          labels[i] = label;
          counts[label] += 1;
        }
        if (counts9 != nullptr) {
          std::memcpy(counts9, counts, sizeof(counts));
        }
        return SEAMS_OK;
      }
    }
    auto rows = bond_rows(work, cutoff, typ, graph);
    const char *spec = rule == SEAMS_CHILL_OVITO ? "ovito-chill+" : "chill+";
    chill::applyChillPlus(work, rows, chill::chillPlusRule(spec), false);
    const int n = work.nop;
    int counts[9] = {};
    for (int i = 0; i < n; ++i) {
      const auto &point = work.pts[static_cast<std::size_t>(i)];
      int label = static_cast<int>(point.iceType);
      if (type > 0 && point.type != type) {
        label = SEAMS_LABEL_UNCLASSIFIED;
      }
      if (label < 0 || label > 8) {
        label = SEAMS_LABEL_UNCLASSIFIED;
      }
      labels[i] = label;
      counts[label] += 1;
    }
    if (counts9 != nullptr) {
      std::memcpy(counts9, counts, sizeof(counts));
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_q4q6(seams_system *sys, double cutoff, int type, double *q4,
                        double *q6) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (q4 == nullptr && q6 == nullptr) {
    return fail(SEAMS_ERR_NULL, "q4 and q6 are both null");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    const int n = work.nop;
    if (q4 != nullptr) {
      const auto result = chill::steinhardtQlCutoff(work, cutoff, 4, typ, false);
      if (static_cast<int>(result.ql.size()) < n) {
        return fail(SEAMS_ERR_INTERNAL, "q4 length is short");
      }
      std::memcpy(q4, result.ql.data(), static_cast<std::size_t>(n) * sizeof(double));
    }
    if (q6 != nullptr) {
      const auto result = chill::steinhardtQlCutoff(work, cutoff, 6, typ, false);
      if (static_cast<int>(result.ql.size()) < n) {
        return fail(SEAMS_ERR_INTERNAL, "q6 length is short");
      }
      std::memcpy(q6, result.ql.data(), static_cast<std::size_t>(n) * sizeof(double));
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_rings_find(seams_system *sys, double cutoff, int type, int graph,
                              int max_depth, seams_rings **out) {
  handle_guard guard(sys);
  clear_error();
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "rings output is null");
  }
  *out = nullptr;
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (max_depth < 3 || max_depth > 12) {
    return fail(SEAMS_ERR_RANGE, "max_depth must be from 3 to 12");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    std::vector<std::vector<int>> idx;
    std::vector<std::vector<int>> found;
    if (graph == SEAMS_GRAPH_CUTOFF) {
      std::vector<int> subset;
      const std::vector<int> *sub = nullptr;
      if (type > 0) {
        subset.reserve(static_cast<std::size_t>(work.nop));
        for (int i = 0; i < work.nop; ++i) {
          if (work.pts[static_cast<std::size_t>(i)].type == typ) {
            subset.push_back(i);
          }
        }
        sub = &subset;
      }
      std::vector<int> offsets;
      std::vector<int> cols;
      if (nneigh::orthoFullCsr(work, sub, cutoff, offsets, cols, false)) {
        found = primitive::ringNetwork(offsets, cols, max_depth);
        idx.assign(static_cast<std::size_t>(work.nop), {});
        for (int i = 0; i < work.nop && i + 1 < static_cast<int>(offsets.size()); ++i) {
          const int a = offsets[static_cast<std::size_t>(i)];
          const int b = offsets[static_cast<std::size_t>(i) + 1];
          auto &row = idx[static_cast<std::size_t>(i)];
          row.reserve(static_cast<std::size_t>(std::max(0, b - a)));
          for (int k = a; k < b; ++k) {
            row.push_back(cols[static_cast<std::size_t>(k)]);
          }
        }
      }
    }
    if (found.empty() && idx.empty()) {
      auto rows = bond_rows(work, cutoff, typ, graph);
      idx = nneigh::neighbourListByIndex(work, rows);
      if (static_cast<int>(idx.size()) < work.nop) {
        idx.resize(static_cast<std::size_t>(work.nop));
      }
      found = primitive::ringNetwork(idx, max_depth);
    }
    if (static_cast<int>(idx.size()) < work.nop) {
      idx.resize(static_cast<std::size_t>(work.nop));
    }
    auto rings = std::make_unique<seams_rings>();
    rings->n_atoms = work.nop;
    rings->rings = std::move(found);
    rings->symmetry.assign(rings->rings.size(), SEAMS_RING_OTHER);
    std::vector<std::vector<int>> six;
    std::vector<int> six_at;
    for (std::size_t i = 0; i < rings->rings.size(); ++i) {
      const int size = static_cast<int>(rings->rings[i].size());
      if (size == 5) {
        rings->symmetry[i] = SEAMS_RING_PENTAGON;
      } else if (size == 6) {
        rings->symmetry[i] = SEAMS_RING_HEXAGON;
        six.push_back(rings->rings[i]);
        six_at.push_back(static_cast<int>(i));
      }
    }
    if (!six.empty()) {
      try {
        const auto aff = ring::cageAffiliation(six, idx);
        const auto planes = ring::stackingPlanes(six, idx);
        for (std::size_t k = 0; k < six.size(); ++k) {
          const bool hc = k < aff.hc.size() && aff.hc[k];
          const bool ddc = k < aff.ddc.size() && aff.ddc[k];
          const bool basal = k < planes.basal.size() && planes.basal[k];
          const bool equatorial = k < planes.equatorial.size() && planes.equatorial[k];
          int code = SEAMS_RING_HEXAGON;
          if (basal && ddc) {
            code = SEAMS_RING_BOTH_BASAL;
          } else if (basal) {
            code = SEAMS_RING_HC_BASAL;
          } else if (equatorial) {
            code = SEAMS_RING_DDC;
          } else if (hc && ddc) {
            code = SEAMS_RING_BOTH_PRISMATIC;
          } else if (hc) {
            code = SEAMS_RING_HC_PRISMATIC;
          } else if (ddc) {
            code = SEAMS_RING_DDC;
          }
          rings->symmetry[static_cast<std::size_t>(six_at[k])] = code;
        }
      } catch (const std::exception &) {
        /* Size labels stay. A cage pass that rejects the graph is not a failed find. */
      }
    }
    try {
      std::vector<std::vector<int>> bySize[7];
      std::vector<int> at[7];
      for (std::size_t i = 0; i < rings->rings.size(); ++i) {
        const int size = static_cast<int>(rings->rings[i].size());
        if (size >= 3 && size <= 6) {
          bySize[size].push_back(rings->rings[i]);
          at[size].push_back(static_cast<int>(i));
        }
      }
      for (int size = 3; size <= 6; ++size) {
        if (bySize[size].empty()) {
          continue;
        }
        std::vector<ring::strucType> kind(bySize[size].size(), ring::strucType::unclassified);
        std::vector<double> rmsd(static_cast<std::size_t>(std::max(work.nop, 0)), -1.0);
        int perfect = 0;
        int imperfect = 0;
        ring::findPrisms(bySize[size], kind, perfect, imperfect, idx, work, rmsd, false);
        for (std::size_t k = 0; k < kind.size(); ++k) {
          int code = -1;
          if (kind[k] == ring::strucType::Prism) {
            code = SEAMS_RING_PRISM;
          } else if (kind[k] == ring::strucType::deformedPrism) {
            code = SEAMS_RING_DEFORMED_PRISM;
          } else if (kind[k] == ring::strucType::mixedPrismRing) {
            code = SEAMS_RING_MIXED_PRISM;
          }
          if (code < 0) {
            continue;
          }
          const int slot = at[size][k];
          const int now = rings->symmetry[static_cast<std::size_t>(slot)];
          if (now == SEAMS_RING_OTHER || now == SEAMS_RING_PENTAGON ||
              now == SEAMS_RING_HEXAGON) {
            rings->symmetry[static_cast<std::size_t>(slot)] = code;
          }
        }
      }
    } catch (const std::exception &) {
      /* Prism labels stay at the size or cage code already written. */
    }
    *out = rings.release();
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_rings_destroy(seams_rings *rings) {
  try {
    delete rings;
  } catch (...) {
  }
}

int seams_rings_count(const seams_rings *rings) {
  return rings == nullptr ? 0 : static_cast<int>(rings->rings.size());
}

seams_status seams_rings_histogram(const seams_rings *rings, int *counts,
                                   int counts_len) {
  handle_guard guard(const_cast<seams_rings *>(rings));
  clear_error();
  if (rings == nullptr || counts == nullptr) {
    return fail(SEAMS_ERR_NULL, "histogram output is null");
  }
  if (counts_len <= 0) {
    return fail(SEAMS_ERR_RANGE, "histogram length is empty");
  }
  for (int i = 0; i < counts_len; ++i) {
    counts[i] = 0;
  }
  for (const auto &ring : rings->rings) {
    const int size = static_cast<int>(ring.size());
    if (size >= 0 && size < counts_len) {
      counts[size] += 1;
    }
  }
  return SEAMS_OK;
}

seams_status seams_rings_labels(const seams_rings *rings, int *sizes, int *symmetry) {
  handle_guard guard(const_cast<seams_rings *>(rings));
  clear_error();
  if (rings == nullptr) {
    return fail(SEAMS_ERR_NULL, "rings handle is null");
  }
  if (sizes == nullptr && symmetry == nullptr) {
    return fail(SEAMS_ERR_NULL, "sizes and symmetry are both null");
  }
  for (std::size_t i = 0; i < rings->rings.size(); ++i) {
    if (sizes != nullptr) {
      sizes[i] = static_cast<int>(rings->rings[i].size());
    }
    if (symmetry != nullptr) {
      symmetry[i] = rings->symmetry[i];
    }
  }
  return SEAMS_OK;
}

seams_status seams_rings_fraction(const seams_rings *rings, int size, double *fraction) {
  handle_guard guard(const_cast<seams_rings *>(rings));
  clear_error();
  if (rings == nullptr || fraction == nullptr) {
    return fail(SEAMS_ERR_NULL, "fraction output is null");
  }
  if (rings->n_atoms <= 0) {
    *fraction = 0.0;
    return SEAMS_OK;
  }
  std::vector<char> hit(static_cast<std::size_t>(rings->n_atoms), 0);
  for (const auto &ring : rings->rings) {
    if (size > 0 && static_cast<int>(ring.size()) != size) {
      continue;
    }
    for (int vertex : ring) {
      if (vertex >= 0 && vertex < rings->n_atoms) {
        hit[static_cast<std::size_t>(vertex)] = 1;
      }
    }
  }
  int marked = 0;
  for (char bit : hit) {
    marked += bit != 0;
  }
  *fraction = static_cast<double>(marked) / static_cast<double>(rings->n_atoms);
  return SEAMS_OK;
}

seams_status seams_neighbours_copy(const seams_neighbours *nb, int *offsets,
                                   size_t offsets_cap, size_t *n_offsets, int *cols,
                                   size_t cols_cap, size_t *nnz) {
  handle_guard guard(const_cast<seams_neighbours *>(nb));
  clear_error();
  if (nb == nullptr || n_offsets == nullptr || nnz == nullptr) {
    return fail(SEAMS_ERR_NULL, "csr output is null");
  }
  const std::size_t need_off = nb->offsets.size();
  const std::size_t need_col = nb->cols.size();
  *n_offsets = need_off;
  *nnz = need_col;
  if ((offsets != nullptr && offsets_cap < need_off) ||
      (cols != nullptr && cols_cap < need_col)) {
    return fail(SEAMS_ERR_RANGE, "buffer is short");
  }
  if (offsets != nullptr && need_off > 0) {
    std::memcpy(offsets, nb->offsets.data(), need_off * sizeof(int));
  }
  if (cols != nullptr && need_col > 0) {
    std::memcpy(cols, nb->cols.data(), need_col * sizeof(int));
  }
  return SEAMS_OK;
}

seams_status seams_rings_vertices(const seams_rings *rings, int index, int *atoms,
                                  size_t cap, size_t *n_out) {
  handle_guard guard(const_cast<seams_rings *>(rings));
  clear_error();
  if (rings == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "ring vertex output is null");
  }
  if (index < 0 || static_cast<std::size_t>(index) >= rings->rings.size()) {
    return fail(SEAMS_ERR_RANGE, "ring index is out of range");
  }
  const auto &ring = rings->rings[static_cast<std::size_t>(index)];
  *n_out = ring.size();
  if (atoms == nullptr) {
    return SEAMS_OK;
  }
  if (cap < ring.size()) {
    return fail(SEAMS_ERR_RANGE, "buffer is short");
  }
  for (std::size_t i = 0; i < ring.size(); ++i) {
    atoms[i] = ring[i];
  }
  return SEAMS_OK;
}

seams_status seams_path_nframes(const char *path, size_t *nframes) {
  clear_error();
  if (path == nullptr || nframes == nullptr) {
    return fail(SEAMS_ERR_NULL, "path or frame count is null");
  }
  try {
    const std::string file(path);
    if (!gen::file_exists(file)) {
      return fail(SEAMS_ERR_ARG, "cannot read path");
    }
    const std::string ext = extension_of(file);
    if (ext == "lammpstrj" || ext == "dump") {
      const int frames = sinp::nLammpsFrames(file);
      if (frames <= 0) {
        return fail(SEAMS_ERR_ARG, "unreadable dump");
      }
      *nframes = static_cast<std::size_t>(frames);
      return SEAMS_OK;
    }
#ifdef SEAMS_HAS_READCON
    if (ext == "con") {
      const int frames = count_con_frames(file);
      if (frames <= 0) {
        return fail(SEAMS_ERR_ARG, "unreadable trajectory");
      }
      *nframes = static_cast<std::size_t>(frames);
      return SEAMS_OK;
    }
#endif
    const int steps = sinp::chemfilesSteps(file);
    if (steps > 0) {
      *nframes = static_cast<std::size_t>(steps);
      return SEAMS_OK;
    }
    if (steps == 0) {
      return fail(SEAMS_ERR_RANGE, "file has no frames");
    }
    int native = -1;
    if (ext == "gro") {
      native = count_gro_frames(file);
    } else if (ext == "pdb") {
      native = count_pdb_frames(file);
    } else if (ext == "xyz") {
      native = count_xyz_frames(file);
    }
    if (native <= 0) {
      return fail(SEAMS_ERR_ARG, "unreadable trajectory");
    }
    *nframes = static_cast<std::size_t>(native);
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::exception &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_system_read(seams_system *sys, const char *path, int frame,
                               int type_filter) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || path == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or path is null");
  }
  if (frame < 1) {
    return fail(SEAMS_ERR_RANGE, "frame index starts at 1");
  }
  try {
    Cloud cloud = sinp::openFrame(path, frame, type_filter);
    if (cloud.nop <= 0) {
      return fail(SEAMS_ERR_RANGE, "frame has no atoms");
    }
    if (cloud.currentFrame <= 0) {
      cloud.currentFrame = frame;
    }
    sys->cloud = std::move(cloud);
    sys->ready = true;
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::exception &ex) {
    const std::string message = ex.what();
    const bool past = message.find("past the last") != std::string::npos;
    g_error = message.empty() ? "cannot read frame" : message;
    note_error();
    return past ? SEAMS_ERR_RANGE : SEAMS_ERR_ARG;
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_chill(seams_system *sys, double cutoff, int type, int graph, int rule,
                         int *labels, size_t labels_len, int counts9[9]) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || labels == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or labels is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!buffer_fits(labels_len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "label buffer is short");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  if (graph != SEAMS_GRAPH_CUTOFF && graph != SEAMS_GRAPH_MUTUAL) {
    return fail(SEAMS_ERR_ARG, "graph is cutoff or mutual");
  }
  if (rule != SEAMS_CHILL_PUBLISHED && rule != SEAMS_CHILL_OVITO) {
    return fail(SEAMS_ERR_ARG, "rule is published or ovito");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    auto rows = bond_rows(work, cutoff, typ, graph);
    if (rule == SEAMS_CHILL_OVITO) {
      chill::BondClassifier classifier = chill::chillRule();
      chill::classifyBonds(work, rows, classifier, false, false);
    } else {
      chill::getCorrel(work, rows, false);
    }
    chill::getIceTypeNoPrint(work, rows, false);
    write_ice_labels(work, type, labels, counts9);
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_ql(seams_system *sys, double cutoff, int type, int order_l, int graph,
                      double *ql, size_t ql_len, double *qlbar, size_t qlbar_len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (ql == nullptr && qlbar == nullptr) {
    return fail(SEAMS_ERR_NULL, "ql and qlbar are both null");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  if (order_l != 3 && order_l != 4 && order_l != 6 && order_l != 8 && order_l != 12) {
    return fail(SEAMS_ERR_RANGE, "order is 3, 4, 6, 8, or 12");
  }
  if (graph != SEAMS_GRAPH_CUTOFF && graph != SEAMS_GRAPH_MUTUAL) {
    return fail(SEAMS_ERR_ARG, "graph is cutoff or mutual");
  }
  if (ql != nullptr && !buffer_fits(ql_len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "ql buffer is short");
  }
  if (qlbar != nullptr && !buffer_fits(qlbar_len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "qlbar buffer is short");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    chill::SteinhardtQl result;
    if (graph == SEAMS_GRAPH_CUTOFF) {
      result = chill::steinhardtQlCutoff(work, cutoff, order_l, typ, qlbar != nullptr);
    } else if (work.box.size() == 3 || work.box.size() >= 6) {
      result = chill::steinhardtQlMutual(work, 4, cutoff + 1.5, typ, order_l);
    } else {
      auto rows = bond_rows(work, cutoff, typ, graph);
      result = chill::steinhardtQl(work, rows, order_l);
    }
    const int n = work.nop;
    if (ql != nullptr) {
      if (static_cast<int>(result.ql.size()) < n) {
        return fail(SEAMS_ERR_INTERNAL, "ql length is short");
      }
      if (n > 0) {
        std::memcpy(ql, result.ql.data(), static_cast<std::size_t>(n) * sizeof(double));
      }
    }
    if (qlbar != nullptr) {
      if (static_cast<int>(result.qlBar.size()) < n) {
        return fail(SEAMS_ERR_INTERNAL, "qlbar length is short");
      }
      if (n > 0) {
        std::memcpy(qlbar, result.qlBar.data(), static_cast<std::size_t>(n) * sizeof(double));
      }
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_q6bar(seams_system *sys, double cutoff, int type, double ice_cut,
                         double *q6bar, size_t q6bar_len, double *mean, double *median,
                         int *n_ice) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (q6bar == nullptr && mean == nullptr && median == nullptr && n_ice == nullptr) {
    return fail(SEAMS_ERR_NULL, "q6bar outputs are null");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff) || !std::isfinite(ice_cut)) {
    return fail(SEAMS_ERR_RANGE, "cutoff and ice cut must be finite");
  }
  if (q6bar != nullptr && !buffer_fits(q6bar_len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "q6bar buffer is short");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    const auto result = chill::steinhardtQlCutoff(work, cutoff, 6, typ, true);
    const int n = work.nop;
    if (static_cast<int>(result.qlBar.size()) < n) {
      return fail(SEAMS_ERR_INTERNAL, "q6bar length is short");
    }
    if (q6bar != nullptr && n > 0) {
      std::memcpy(q6bar, result.qlBar.data(), static_cast<std::size_t>(n) * sizeof(double));
    }
    std::vector<double> finite;
    finite.reserve(static_cast<std::size_t>(std::max(n, 0)));
    double sum = 0.0;
    int ice = 0;
    for (int i = 0; i < n; ++i) {
      const double value = result.qlBar[static_cast<std::size_t>(i)];
      if (!std::isfinite(value)) {
        continue;
      }
      sum += value;
      finite.push_back(value);
      if (value > ice_cut) {
        ice += 1;
      }
    }
    if (mean != nullptr) {
      *mean = finite.empty() ? 0.0 : sum / static_cast<double>(finite.size());
    }
    if (median != nullptr) {
      *median = median_of(std::move(finite));
    }
    if (n_ice != nullptr) {
      *n_ice = ice;
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_cages_find(seams_system *sys, double cutoff, int type, int graph, int k,
                              int complete, seams_cage_set **out) {
  handle_guard guard(sys);
  clear_error();
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "cages output is null");
  }
  *out = nullptr;
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  if (graph != SEAMS_GRAPH_CUTOFF && graph != SEAMS_GRAPH_MUTUAL &&
      graph != SEAMS_GRAPH_SEEDED) {
    return fail(SEAMS_ERR_ARG, "graph is cutoff, mutual, or seeded");
  }
  if (k < 0) {
    return fail(SEAMS_ERR_RANGE, "neighbour count is negative");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    const int n = work.nop;
    const int neighbours = k > 0 ? k : 4;
    std::vector<char> hc(static_cast<std::size_t>(std::max(n, 0)), 0);
    std::vector<char> ddc(hc.size(), 0);
    if (graph == SEAMS_GRAPH_SEEDED) {
      const double candidate = cutoff + 1.5;
      auto graphs = nneigh::kNearestNeighbourPair(work, neighbours, candidate, typ);
      auto idx_strict = nneigh::neighbourListByIndex(work, graphs.first);
      auto idx_union = nneigh::neighbourListByIndex(work, graphs.second);
      const bool same = graphs.first == graphs.second ||
                        seams::provenEqual(seams::edgeSetOf(idx_strict),
                                           seams::edgeSetOf(idx_union));
      ring::SeededAtomLabels aff;
      if (same) {
        auto six = six_of(ring::contractibleRings(primitive::ringNetwork(idx_strict, 6), work));
        const auto once = ring::cageAffiliation(six, idx_strict);
        aff = ring::seededCageAffiliation(six, once, six, once, idx_strict, complete != 0);
      } else {
        auto six_strict =
            six_of(ring::contractibleRings(primitive::ringNetwork(idx_strict, 6), work));
        auto six_union =
            six_of(ring::contractibleRings(primitive::ringNetwork(idx_union, 6), work));
        aff = ring::seededCageAffiliation(six_strict, idx_strict, six_union, idx_union,
                                          complete != 0);
      }
      const int m = std::min(n, static_cast<int>(aff.hc.size()));
      for (int i = 0; i < m; ++i) {
        hc[static_cast<std::size_t>(i)] = aff.hc[static_cast<std::size_t>(i)];
        ddc[static_cast<std::size_t>(i)] = aff.ddc[static_cast<std::size_t>(i)];
      }
    } else {
      auto rows = bond_rows(work, cutoff, typ, graph);
      auto idx = nneigh::neighbourListByIndex(work, rows);
      if (static_cast<int>(idx.size()) < n) {
        idx.resize(static_cast<std::size_t>(n));
      }
      if (tum::preferOffload()) {
        const auto counts = tum::cageCounts(idx);
        const int m = std::min(n, static_cast<int>(counts.atomHc.size()));
        for (int i = 0; i < m; ++i) {
          hc[static_cast<std::size_t>(i)] = counts.atomHc[static_cast<std::size_t>(i)] != 0;
          ddc[static_cast<std::size_t>(i)] = counts.atomDdc[static_cast<std::size_t>(i)] != 0;
        }
      } else {
        auto six = six_of(ring::contractibleRings(primitive::ringNetwork(idx, 6), work));
        flags_from_affiliation(n, six, ring::cageAffiliation(six, idx), hc, ddc);
      }
    }
    auto cages = std::make_unique<seams_cage_set>();
    tally_cages(cages.get(), work, type, hc, ddc);
    *out = cages.release();
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_cages_destroy(seams_cage_set *cages) {
  try {
    delete cages;
  } catch (...) {
  }
}

seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex, int *n_cubic,
                                int *n_water, int *n_both) {
  handle_guard guard(const_cast<seams_cage_set *>(cages));
  clear_error();
  if (cages == nullptr) {
    return fail(SEAMS_ERR_NULL, "cages handle is null");
  }
  if (n_hex == nullptr && n_cubic == nullptr && n_water == nullptr && n_both == nullptr) {
    return fail(SEAMS_ERR_NULL, "cage counts are null");
  }
  if (n_hex != nullptr) {
    *n_hex = cages->n_hex;
  }
  if (n_cubic != nullptr) {
    *n_cubic = cages->n_cubic;
  }
  if (n_water != nullptr) {
    *n_water = cages->n_water;
  }
  if (n_both != nullptr) {
    *n_both = cages->n_both;
  }
  return SEAMS_OK;
}

seams_status seams_cages_labels(const seams_cage_set *cages, int *per_atom, size_t len) {
  handle_guard guard(const_cast<seams_cage_set *>(cages));
  clear_error();
  if (cages == nullptr || per_atom == nullptr) {
    return fail(SEAMS_ERR_NULL, "cage labels are null");
  }
  if (!buffer_fits(len, cages->n_atoms)) {
    return fail(SEAMS_ERR_RANGE, "label buffer is short");
  }
  if (cages->n_atoms > 0) {
    std::memcpy(per_atom, cages->per_atom.data(),
                static_cast<std::size_t>(cages->n_atoms) * sizeof(int));
  }
  return SEAMS_OK;
}

seams_status seams_cages_signature(seams_system *sys, double cutoff, int type, int graph,
                                   const char *spec, size_t *n_cages, size_t *n_atoms) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || spec == nullptr || n_cages == nullptr || n_atoms == nullptr) {
    return fail(SEAMS_ERR_NULL, "signature needs a frame, a census, and counts");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  if (graph != SEAMS_GRAPH_CUTOFF && graph != SEAMS_GRAPH_MUTUAL) {
    return fail(SEAMS_ERR_ARG, "graph is cutoff or mutual");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    auto rows = bond_rows(work, cutoff, typ, graph);
    auto idx = nneigh::neighbourListByIndex(work, rows);
    const auto signature = cage::Signature::parse(spec);
    const int depth = std::max(signature.maxRingSize(), 3);
    const auto rings = ring::cageRings(primitive::ringNetwork(idx, depth), work, cutoff);
    const auto window = cage::parseEdgeWindow("");
    std::vector<std::array<double, 3>> xyz(static_cast<std::size_t>(std::max(work.nop, 0)));
    for (int i = 0; i < work.nop; ++i) {
      const auto &point = work.pts[static_cast<std::size_t>(i)];
      xyz[static_cast<std::size_t>(i)] = {point.x, point.y, point.z};
    }
    const auto found_of = cage::findBySignatures(rings, idx, std::vector<cage::Signature>{signature},
                                                 {}, &xyz, &work.box, window);
    const auto &found = found_of.at(0);
    std::set<int> atoms;
    for (const auto &cage_i : found) {
      atoms.insert(cage_i.vertices.begin(), cage_i.vertices.end());
    }
    *n_cages = found.size();
    *n_atoms = atoms.size();
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (const std::out_of_range &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd, int *labels,
                                   size_t labels_len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || labels == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or labels is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!buffer_fits(labels_len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "label buffer is short");
  }
  if (shell_rmsd > 0.0 && !std::isfinite(shell_rmsd)) {
    return fail(SEAMS_ERR_RANGE, "shell residual must be finite");
  }
  const int n = sys->cloud.nop;
  if (n <= 0) {
    return SEAMS_OK;
  }
  try {
    use_threads();
    std::vector<double> xyz(static_cast<std::size_t>(n) * 3);
    for (int i = 0; i < n; ++i) {
      const auto &point = sys->cloud.pts[static_cast<std::size_t>(i)];
      const std::size_t base = static_cast<std::size_t>(i) * 3;
      xyz[base] = point.x;
      xyz[base + 1] = point.y;
      xyz[base + 2] = point.z;
    }
    double lattice[9];
    lattice_from_cloud(sys->cloud, lattice);
    const double ceiling =
        shell_rmsd > 0.0 ? shell_rmsd : symm::kShellRmsd;
    const int code = seams_symmetry_classify_tol(xyz.data(), n, lattice, ceiling, labels);
    if (code != 0) {
      const char *message = seams_c_error();
      return fail(SEAMS_ERR_INTERNAL,
                  message != nullptr && message[0] != '\0' ? message : "symmetry classify failed");
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_brine_analyse(seams_system *sys, const int *ice, size_t n_ice,
                                 int water_type, int ion_type, int n_bins, int axis,
                                 seams_brine **out) {
  handle_guard guard(sys);
  clear_error();
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "brine output is null");
  }
  *out = nullptr;
  if (sys == nullptr || ice == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or ice flags are null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!buffer_fits(n_ice, sys->cloud.nop) || n_ice != static_cast<std::size_t>(sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "ice flags must cover the frame");
  }
  if (n_bins < 0) {
    return fail(SEAMS_ERR_RANGE, "bin count is negative");
  }
  if (axis < 0 || axis > 2) {
    return fail(SEAMS_ERR_RANGE, "axis is 0, 1, or 2");
  }
  try {
    use_threads();
    const int n = sys->cloud.nop;
    std::vector<bool> ice_flag(static_cast<std::size_t>(n), false);
    std::vector<int> waters;
    std::vector<int> ions;
    waters.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      ice_flag[static_cast<std::size_t>(i)] = ice[static_cast<std::size_t>(i)] != 0;
      const int kind = sys->cloud.pts[static_cast<std::size_t>(i)].type;
      if (ion_type > 0 && kind == ion_type) {
        ions.push_back(i);
      } else if (water_type <= 0 || kind == water_type) {
        waters.push_back(i);
      }
    }
    site::FrontOptions options;
    if (n_bins > 0) {
      options.nBins = n_bins;
    }
    options.axis = axis;
    auto front = std::make_unique<seams_brine>();
    front->snapshot = site::analyseFront(sys->cloud, ice_flag, waters, ions, options);
    *out = front.release();
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_brine_destroy(seams_brine *front) {
  try {
    delete front;
  } catch (...) {
  }
}

seams_status seams_brine_fronts(const seams_brine *front, double *z_enter, double *z_exit,
                                double *thickness, int *found) {
  handle_guard guard(const_cast<seams_brine *>(front));
  clear_error();
  if (front == nullptr) {
    return fail(SEAMS_ERR_NULL, "brine handle is null");
  }
  if (z_enter == nullptr && z_exit == nullptr && thickness == nullptr && found == nullptr) {
    return fail(SEAMS_ERR_NULL, "front outputs are null");
  }
  const auto &pair = front->snapshot.fronts;
  if (z_enter != nullptr) {
    *z_enter = pair.zEnter;
  }
  if (z_exit != nullptr) {
    *z_exit = pair.zExit;
  }
  if (thickness != nullptr) {
    *thickness = pair.thickness;
  }
  if (found != nullptr) {
    *found = pair.found ? 1 : 0;
  }
  return SEAMS_OK;
}

seams_status seams_brine_rejection(const seams_brine *front, double *c_ice, double *c_liquid,
                                   double *k, int *defined) {
  handle_guard guard(const_cast<seams_brine *>(front));
  clear_error();
  if (front == nullptr) {
    return fail(SEAMS_ERR_NULL, "brine handle is null");
  }
  if (c_ice == nullptr && c_liquid == nullptr && k == nullptr && defined == nullptr) {
    return fail(SEAMS_ERR_NULL, "rejection outputs are null");
  }
  const auto &rejection = front->snapshot.rejection;
  if (c_ice != nullptr) {
    *c_ice = rejection.cIce;
  }
  if (c_liquid != nullptr) {
    *c_liquid = rejection.cLiquid;
  }
  if (defined != nullptr) {
    *defined = rejection.defined ? 1 : 0;
  }
  if (k != nullptr && rejection.defined) {
    *k = rejection.k;
  }
  return SEAMS_OK;
}

seams_status seams_brine_profile(const seams_brine *front, double *z_mid, double *molality,
                                 double *ice_fraction, size_t cap, size_t *n_out) {
  handle_guard guard(const_cast<seams_brine *>(front));
  clear_error();
  if (front == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "profile output is null");
  }
  const auto &bins = front->snapshot.profile.bins;
  *n_out = bins.size();
  if (z_mid == nullptr && molality == nullptr && ice_fraction == nullptr) {
    return SEAMS_OK;
  }
  if (cap < bins.size()) {
    return fail(SEAMS_ERR_RANGE, "buffer is short");
  }
  for (std::size_t i = 0; i < bins.size(); ++i) {
    if (z_mid != nullptr) {
      z_mid[i] = bins[i].zMid;
    }
    if (molality != nullptr) {
      molality[i] = bins[i].molality;
    }
    if (ice_fraction != nullptr) {
      ice_fraction[i] = bins[i].iceFraction;
    }
  }
  return SEAMS_OK;
}

seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph, int k,
                         int complete, seams_cage_set **out) {
  return seams_cages_find(sys, cutoff, type, graph, k, complete, out);
}

seams_status seams_signature_cages(seams_system *sys, double cutoff, int type, int graph,
                                   const char *spec, size_t *n_cages, size_t *n_atoms) {
  return seams_cages_signature(sys, cutoff, type, graph, spec, n_cages, n_atoms);
}

seams_status seams_steinhardt(seams_system *sys, double cutoff, int type, int order_l,
                              double *ql, size_t ql_len, double *qlbar, size_t qlbar_len) {
  return seams_ql(sys, cutoff, type, order_l, SEAMS_GRAPH_CUTOFF, ql, ql_len, qlbar,
                  qlbar_len);
}

namespace {

seams_status copy_text(const std::string &text, char *buf, size_t cap, size_t *n_out) {
  if (n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "text length is null");
  }
  *n_out = text.size();
  if (buf == nullptr) {
    return SEAMS_OK;
  }
  if (cap < text.size() + 1) {
    return fail(SEAMS_ERR_RANGE, "buffer is short");
  }
  if (!text.empty()) {
    std::memcpy(buf, text.data(), text.size());
  }
  buf[text.size()] = '\0';
  return SEAMS_OK;
}

rdf::PartialRdf partial_of(const Cloud &cloud, double rmax, int bins, int type_i, int type_j) {
  int bins_used = bins;
  if (bins_used <= 0) {
    bins_used = std::max(1, static_cast<int>(std::lround(rmax / 0.1)));
  }
  int typ = 1;
  Cloud work = unified_copy(cloud, type_i, typ);
  const int j = type_j > 0 ? type_j : typ;
  return rdf::partialRdf(work, typ, j, rmax, bins_used);
}

} // namespace

seams_status seams_neighbours_knn(seams_system *sys, int k, double cutoff, int type,
                                  int mutual, seams_neighbours **out) {
  handle_guard guard(sys);
  clear_error();
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "neighbours output is null");
  }
  *out = nullptr;
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (k <= 0) {
    return fail(SEAMS_ERR_RANGE, "neighbour count must be positive");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  try {
    use_threads();
    int typ = 1;
    const Cloud *src = &sys->cloud;
    Cloud owned;
    if (type > 0) {
      typ = type;
      bool uniform = sys->cloud.nop > 0;
      for (int i = 0; i < sys->cloud.nop && uniform; ++i) {
        if (sys->cloud.pts[static_cast<std::size_t>(i)].type != type) {
          uniform = false;
        }
      }
      if (!uniform) {
        owned = unified_copy(sys->cloud, type, typ);
        src = &owned;
      }
    } else {
      owned = unified_copy(sys->cloud, type, typ);
      src = &owned;
    }
    auto nb = std::make_unique<seams_neighbours>();
    nb->n = src->nop;
    nneigh::kNearestCsr(*src, k, cutoff + 1.5, typ, mutual != 0, nb->offsets, nb->cols);
    *out = nb.release();
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i, int type_j,
                       double *r, double *g, int *count, size_t cap, size_t *n_out,
                       double *volume) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "rdf output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(rmax > 0.0) || !std::isfinite(rmax)) {
    return fail(SEAMS_ERR_RANGE, "rmax must be positive");
  }
  try {
    use_threads();
    const auto gr = partial_of(sys->cloud, rmax, bins, type_i, type_j);
    *n_out = gr.r.size();
    if (volume != nullptr) {
      *volume = gr.volume;
    }
    if (r == nullptr && g == nullptr && count == nullptr) {
      return SEAMS_OK;
    }
    if (cap < gr.r.size()) {
      return fail(SEAMS_ERR_RANGE, "buffer is short");
    }
    for (std::size_t i = 0; i < gr.r.size(); ++i) {
      if (r != nullptr) {
        r[i] = gr.r[i];
      }
      if (g != nullptr) {
        g[i] = gr.g[i];
      }
      if (count != nullptr) {
        count[i] = gr.count[i];
      }
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i,
                              int type_j, double *r, double *cn, size_t cap, size_t *n_out) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "running coordination output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(rmax > 0.0) || !std::isfinite(rmax)) {
    return fail(SEAMS_ERR_RANGE, "rmax must be positive");
  }
  try {
    use_threads();
    const auto gr = partial_of(sys->cloud, rmax, bins, type_i, type_j);
    const double rho = gr.volume > 0.0 ? static_cast<double>(gr.nJ) / gr.volume : 0.0;
    const auto series = rdf::runningCN(gr, rho);
    *n_out = series.size();
    if (r == nullptr && cn == nullptr) {
      return SEAMS_OK;
    }
    if (cap < series.size() || gr.r.size() < series.size()) {
      return fail(SEAMS_ERR_RANGE, "buffer is short");
    }
    for (std::size_t i = 0; i < series.size(); ++i) {
      if (r != nullptr) {
        r[i] = gr.r[i];
      }
      if (cn != nullptr) {
        cn[i] = series[i];
      }
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_coordination(seams_system *sys, double rmax, int bins, int type_i,
                                int type_j, double *cn) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || cn == nullptr) {
    return fail(SEAMS_ERR_NULL, "coordination output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(rmax > 0.0) || !std::isfinite(rmax)) {
    return fail(SEAMS_ERR_RANGE, "rmax must be positive");
  }
  try {
    use_threads();
    const auto gr = partial_of(sys->cloud, rmax, bins, type_i, type_j);
    const double rho = gr.volume > 0.0 ? static_cast<double>(gr.nJ) / gr.volume : 0.0;
    *cn = rdf::coordinationNumber(gr, rmax, rho);
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type,
                          int hydrogen_type, double dist_cutoff, double angle_cutoff,
                          int *n_bonds) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_bonds == nullptr) {
    return fail(SEAMS_ERR_NULL, "hbond output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !(dist_cutoff > 0.0) || !(angle_cutoff > 0.0)) {
    return fail(SEAMS_ERR_RANGE, "cutoffs must be positive");
  }
  if (hydrogen_type <= 0) {
    return fail(SEAMS_ERR_ARG, "hydrogen type must be positive");
  }
  try {
    use_threads();
    Cloud heavy;
    Cloud hydro;
    heavy.box = sys->cloud.box;
    heavy.boxLow = sys->cloud.boxLow;
    hydro.box = sys->cloud.box;
    hydro.boxLow = sys->cloud.boxLow;
    for (int i = 0; i < sys->cloud.nop; ++i) {
      const auto &point = sys->cloud.pts[static_cast<std::size_t>(i)];
      if (point.type == hydrogen_type) {
        hydro.pts.push_back(point);
      } else if (oxygen_type <= 0 || point.type == oxygen_type) {
        heavy.pts.push_back(point);
      }
    }
    heavy.nop = static_cast<int>(heavy.pts.size());
    hydro.nop = static_cast<int>(hydro.pts.size());
    for (int i = 0; i < heavy.nop; ++i) {
      heavy.idIndexMap[heavy.pts[static_cast<std::size_t>(i)].atomID] = i;
    }
    for (int i = 0; i < hydro.nop; ++i) {
      hydro.idIndexMap[hydro.pts[static_cast<std::size_t>(i)].atomID] = i;
    }
    if (heavy.nop == 0 || hydro.nop == 0) {
      *n_bonds = 0;
      return SEAMS_OK;
    }
    const int typ = heavy.pts[0].type;
    auto rows = nneigh::neighListO(cutoff, heavy, typ);
    auto net = bond::populateHbondsWithInputClouds(heavy, hydro, rows, dist_cutoff,
                                                   angle_cutoff);
    int edges = 0;
    for (const auto &row : net) {
      if (row.size() > 1) {
        edges += static_cast<int>(row.size()) - 1;
      }
    }
    *n_bonds = edges / 2;
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis,
                                double *z, double *rho, size_t cap, size_t *n_out) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "density output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (axis < 0 || axis > 2) {
    return fail(SEAMS_ERR_RANGE, "axis is 0, 1, or 2");
  }
  try {
    int bins_used = bins;
    if (bins_used <= 0) {
      const double span = static_cast<int>(sys->cloud.box.size()) > axis
                              ? sys->cloud.box[static_cast<std::size_t>(axis)]
                              : 0.0;
      bins_used = std::max(1, static_cast<int>(std::lround(span / 0.1)));
    }
    const auto density = site::densityZ(sys->cloud, type, bins_used, axis);
    *n_out = density.z.size();
    if (z == nullptr && rho == nullptr) {
      return SEAMS_OK;
    }
    if (cap < density.z.size()) {
      return fail(SEAMS_ERR_RANGE, "buffer is short");
    }
    for (std::size_t i = 0; i < density.z.size(); ++i) {
      if (z != nullptr) {
        z[i] = density.z[i];
      }
      if (rho != nullptr) {
        rho[i] = density.rho[i];
      }
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k, int graph,
                               int hops, int colour_types, seams_topo **out) {
  handle_guard guard(sys);
  clear_error();
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "fingerprint output is null");
  }
  *out = nullptr;
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  if (graph != SEAMS_GRAPH_CUTOFF && graph != SEAMS_GRAPH_MUTUAL) {
    return fail(SEAMS_ERR_ARG, "graph is cutoff or mutual");
  }
  if (hops <= 0) {
    hops = 2;
  }
  if (k <= 0) {
    k = 4;
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    std::vector<std::vector<int>> rows;
    if (graph == SEAMS_GRAPH_CUTOFF) {
      rows = bond_rows(work, cutoff, typ, graph);
    } else {
      rows = nneigh::kNearestNeighbourList(work, k, cutoff + 1.5, typ, true);
    }
    auto idx = nneigh::neighbourListByIndex(work, rows);
    std::vector<int> colours;
    if (colour_types != 0) {
      colours.reserve(static_cast<std::size_t>(work.nop));
      for (int i = 0; i < work.nop; ++i) {
        colours.push_back(work.pts[static_cast<std::size_t>(i)].type);
      }
    }
    const auto fp = topo::fingerprint(idx, hops, 7, colours);
    auto handle = std::make_unique<seams_topo>();
    handle->frame_key = fp.key;
    handle->method = fp.method;
    handle->atom_keys = fp.atomKeys;
    handle->ring_census = fp.ringCensus;
    handle->hops = fp.hops;
    *out = handle.release();
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_fingerprint_destroy(seams_topo *fp) {
  try {
    delete fp;
  } catch (...) {
  }
}

seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap, size_t *n_out) {
  handle_guard guard(const_cast<seams_topo *>(fp));
  clear_error();
  if (fp == nullptr) {
    return fail(SEAMS_ERR_NULL, "fingerprint handle is null");
  }
  return copy_text(fp->frame_key, buf, cap, n_out);
}

seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf, size_t cap,
                                    size_t *n_out) {
  handle_guard guard(const_cast<seams_topo *>(fp));
  clear_error();
  if (fp == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "fingerprint atom output is null");
  }
  if (index < 0 || static_cast<std::size_t>(index) >= fp->atom_keys.size()) {
    return fail(SEAMS_ERR_RANGE, "atom index is out of range");
  }
  return copy_text(fp->atom_keys[static_cast<std::size_t>(index)], buf, cap, n_out);
}

seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts, size_t counts_len) {
  handle_guard guard(const_cast<seams_topo *>(fp));
  clear_error();
  if (fp == nullptr || counts == nullptr) {
    return fail(SEAMS_ERR_NULL, "fingerprint ring output is null");
  }
  for (std::size_t i = 0; i < counts_len; ++i) {
    counts[i] = i < fp->ring_census.size() ? fp->ring_census[i] : 0;
  }
  return SEAMS_OK;
}

seams_status seams_contact_pairs(seams_system *sys, const char *site_spec, size_t *n_pairs,
                                 int *n_cation, int *n_anion) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || site_spec == nullptr) {
    return fail(SEAMS_ERR_NULL, "contact pairs need a frame and a site spec");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (n_pairs == nullptr && n_cation == nullptr && n_anion == nullptr) {
    return fail(SEAMS_ERR_NULL, "contact pair outputs are null");
  }
  try {
    const auto table = site::parseSiteSpec(site_spec);
    const auto ions = site::ionCloud(sys->cloud, table);
    const auto pairs = nneigh::mutualNearestUnlike(ions, 1, 2);
    int cations = 0;
    int anions = 0;
    for (const auto &point : ions.pts) {
      if (point.type == 1) {
        ++cations;
      } else if (point.type == 2) {
        ++anions;
      }
    }
    if (n_pairs != nullptr) {
      *n_pairs = pairs.size();
    }
    if (n_cation != nullptr) {
      *n_cation = cations;
    }
    if (n_anion != nullptr) {
      *n_anion = anions;
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec, int subset,
                           int *n_subset, int *largest, double *percolation) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || site_spec == nullptr) {
    return fail(SEAMS_ERR_NULL, "domains need a frame and a site spec");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (n_subset == nullptr && largest == nullptr && percolation == nullptr) {
    return fail(SEAMS_ERR_NULL, "domain outputs are null");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  site::Kind kind = site::Kind::apolar;
  if (subset == 1) {
    kind = site::Kind::polar;
  } else if (subset != 0) {
    return fail(SEAMS_ERR_ARG, "subset is 0 (apolar) or 1 (polar)");
  }
  try {
    use_threads();
    const auto table = site::parseSiteSpec(site_spec);
    std::vector<bool> mask(static_cast<std::size_t>(sys->cloud.nop), false);
    for (const int index : site::indicesOf(sys->cloud, table, kind)) {
      if (index >= 0 && index < sys->cloud.nop) {
        mask[static_cast<std::size_t>(index)] = true;
      }
    }
    const auto idx = nneigh::getNewNeighbourListByIndex(sys->cloud, cutoff);
    std::vector<std::vector<int>> ids(idx.size());
    for (std::size_t i = 0; i < idx.size(); ++i) {
      for (const int index : idx[i]) {
        if (index >= 0 && index < sys->cloud.nop) {
          ids[i].push_back(sys->cloud.pts[static_cast<std::size_t>(index)].atomID);
        }
      }
    }
    const auto domain = clump::largestDomain(sys->cloud, ids, mask);
    if (n_subset != nullptr) {
      *n_subset = domain.subset;
    }
    if (largest != nullptr) {
      *largest = domain.largest;
    }
    if (percolation != nullptr) {
      *percolation = domain.percolation;
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect,
                          int *imperfect, size_t len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || perfect == nullptr || imperfect == nullptr) {
    return fail(SEAMS_ERR_NULL, "prism output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (len < 7) {
    return fail(SEAMS_ERR_RANGE, "prism buffers need length 7");
  }
  if (!(cutoff > 0.0) || !std::isfinite(cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoff must be positive");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    auto rows = nneigh::kNearestNeighbourList(work, 4, cutoff + 1.5, typ, true);
    auto idx = nneigh::neighbourListByIndex(work, rows);
    const auto rings = primitive::ringNetwork(idx, 6);
    for (std::size_t i = 0; i < len; ++i) {
      perfect[i] = 0;
      imperfect[i] = 0;
    }
    for (int size = 3; size <= 6; ++size) {
      std::vector<std::vector<int>> one;
      for (const auto &ring : rings) {
        if (static_cast<int>(ring.size()) == size) {
          one.push_back(ring);
        }
      }
      std::vector<ring::strucType> kind(one.size());
      std::vector<double> rmsd(static_cast<std::size_t>(std::max(work.nop, 0)), -1.0);
      int n_perfect = 0;
      int n_imperfect = 0;
      if (!one.empty()) {
        ring::findPrisms(one, kind, n_perfect, n_imperfect, idx, work, rmsd, false);
      }
      perfect[size] = n_perfect;
      imperfect[size] = n_imperfect;
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_monolayer(seams_system *sys, double cutoff, int type, int max_depth,
                             double sheet_area, int *counts, size_t counts_len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || counts == nullptr) {
    return fail(SEAMS_ERR_NULL, "monolayer output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (max_depth < 3 || counts_len < static_cast<std::size_t>(max_depth) + 1) {
    return fail(SEAMS_ERR_RANGE, "monolayer needs max_depth >= 3 and a long counts buffer");
  }
  if (!(cutoff > 0.0) || !(sheet_area > 0.0)) {
    return fail(SEAMS_ERR_RANGE, "cutoff and sheet area must be positive");
  }
  try {
    use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    auto rows = bond_rows(work, cutoff, typ, SEAMS_GRAPH_CUTOFF);
    auto idx = nneigh::neighbourListByIndex(work, rows);
    const auto rings = primitive::ringNetwork(idx, max_depth);
    for (std::size_t i = 0; i < counts_len; ++i) {
      counts[i] = 0;
    }
    for (const auto &ring : rings) {
      const std::size_t size = ring.size();
      if (size < counts_len) {
        counts[size] += 1;
      }
    }
    (void)sheet_area;
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds, size_t len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || spec == nullptr || kinds == nullptr) {
    return fail(SEAMS_ERR_NULL, "site table output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!buffer_fits(len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "kind buffer is short");
  }
  try {
    const auto table = site::parseSiteSpec(spec);
    for (int i = 0; i < sys->cloud.nop; ++i) {
      kinds[i] = static_cast<int>(table.of(sys->cloud.pts[static_cast<std::size_t>(i)]));
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::invalid_argument &ex) {
    return fail(SEAMS_ERR_ARG, ex.what());
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_brine_pockets(const seams_brine *front, int *n_pockets,
                                 int *reservoir_waters) {
  handle_guard guard(const_cast<seams_brine *>(front));
  clear_error();
  if (front == nullptr) {
    return fail(SEAMS_ERR_NULL, "brine handle is null");
  }
  if (n_pockets == nullptr && reservoir_waters == nullptr) {
    return fail(SEAMS_ERR_NULL, "pocket outputs are null");
  }
  if (n_pockets != nullptr) {
    *n_pockets = front->snapshot.pockets.nPockets;
  }
  if (reservoir_waters != nullptr) {
    *reservoir_waters = front->snapshot.pockets.reservoirWaters;
  }
  return SEAMS_OK;
}

} // extern "C"
