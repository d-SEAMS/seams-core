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
#include <seams_input.hpp>
#include <topo_bulk.hpp>
#include <topo_one_dim.hpp>
#include <topo_two_dim.hpp>
#include <tum_offload.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <exception>
#include <fstream>
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

std::string extension_of(const std::string &path) {
  const auto dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return {};
  }
  return path.substr(dot + 1);
}

Cloud load_path(const std::string &path, int frame, int typeI) {
  Cloud cloud;
  const auto ext = extension_of(path);
  if (ext == "xyz") {
    return sinp::readXYZ(path);
  }
  if (ext == "con") {
#ifdef SEAMS_HAS_READCON
    // A missing path is an error. readCon clears the cloud and returns
    // it, which would look like an empty frame.
    if (!gen::file_exists(path)) {
      throw std::runtime_error("cannot read " + path);
    }
    return sinp::readCon(path, frame, cloud);
#else
    throw std::runtime_error(
        "readcon-core is not in this build; CON input needs it");
#endif
  }
  if (ext == "pdb" || ext == "gro" || ext == "dcd") {
#ifdef SEAMS_HAS_CHEMFILES
    const int filter = typeI > 0 ? typeI : -1;
    return sinp::readChemfiles(path, frame, cloud, filter);
#else
    throw std::runtime_error("chemfiles is not in this build; " + ext +
                             " input needs it");
#endif
  }
  if (typeI < 0) {
    return sinp::readLammpsTrj(path, frame, cloud);
  }
  if (typeI > 0) {
    return sinp::readLammpsTrjO(path, frame, cloud, typeI);
  }
  cloud = sinp::readLammpsTrjO(path, frame, cloud, 2);
  if (cloud.nop > 0) {
    return cloud;
  }
  Cloud again;
  return sinp::readLammpsTrjO(path, frame, again, 1);
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
    const auto ext = extension_of(path);
    if (ext == "xyz") {
      std::ifstream in(path);
      if (!in) {
        set_error("cannot read " + std::string(path));
        return -1;
      }
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

int classify_chill(seams_frame *frame, double cutoff, int type,
                   const char *graph, bool plus, const char *rule, int *labels,
                   int counts9[9]) {
  clear_error();
  if (frame == nullptr) {
    set_error("frame is null");
    return 1;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    const std::string which =
        (graph == nullptr || graph[0] == '\0') ? "cutoff" : graph;
    if (which != "cutoff" && which != "mutual") {
      throw std::invalid_argument("unknown CHILL graph '" + which +
                                  "' (use cutoff or mutual)");
    }
    const bool publishedRule =
        rule == nullptr || rule[0] == '\0' || std::string(rule) == "chill+";
    // Published cutoff CHILL+ reads its own pair rows. Other graphs and
    // the OVITO rule keep the neighbour list the caller asked for.
    if (!(plus && which == "cutoff" && publishedRule &&
          chill::chillPlusPublishedCutoff(frame->cloud, cutoff, typ))) {
      auto nList = bond_list(frame->cloud, cutoff, typ, which, 4);
      if (plus) {
        chill::applyChillPlus(frame->cloud, nList,
                              chill::chillPlusRule(rule == nullptr ? "" : rule),
                              false);
      } else {
        chill::getCorrel(frame->cloud, nList, false);
        chill::getIceTypeNoPrint(frame->cloud, nList, false);
      }
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

} // namespace

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
  try {
    const int typ = resolve_type(frame->cloud, type);
    const std::string which = graph == nullptr ? "mutual" : graph;
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
  clear_error();
  if (frame == nullptr) {
    set_error("frame is null");
    return 1;
  }
  try {
    const int typ = resolve_type(frame->cloud, type);
    // Orthorhombic and restricted-triclinic frames use the cutoff shell,
    // the same neighbours OVITO's cutoff Steinhardt reads, built by the
    // shared CSR. Anything else keeps the mutual four-nearest graph.
    chill::SteinhardtQl result;
    if (frame->cloud.box.size() == 3 || frame->cloud.box.size() >= 6) {
      result = chill::steinhardtQlCutoff(frame->cloud, cutoff, order_l, typ, true);
    } else {
      const double cand = cutoff + 1.5;
      result = chill::steinhardtQlMutual(frame->cloud, 4, cand, typ, order_l);
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
