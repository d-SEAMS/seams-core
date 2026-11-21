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
#include <cpu_features.hpp>
#include <bop.hpp>
#include <brine_front.hpp>
#include <chordless.hpp>
#include <cluster.hpp>
#include <density.hpp>
#include <cage_affiliation.hpp>
#include <cage_enum.hpp>
#include <chill_plus_scheme.hpp>
#include <edge_set.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <index_range.hpp>
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
#include <trajectory.hpp>
#include <tum_offload.hpp>

#include <array>
#include <algorithm>
#include <charconv>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

/* The error slot is the only state a const-handle call may write, and only
 * when it fails. Several threads may run const calls on one handle, so
 * those writes take the lock. */
struct seams_object {
  mutable std::string error;
  mutable std::mutex error_mu;
};

thread_local std::string g_error;
thread_local const seams_object *g_current = nullptr;
int g_threads = 0;
std::mutex g_live_mu;
std::unordered_set<void *> g_live;

void *track(void *ptr) {
  if (ptr == nullptr) {
    return nullptr;
  }
  std::lock_guard<std::mutex> lock(g_live_mu);
  g_live.insert(ptr);
  return ptr;
}

/* Storing a message must not throw out of a C entry point. On a failed
 * allocation the message is dropped and the status code still returns. */
void set_error(const char *message) noexcept {
  try {
    g_error = message == nullptr ? "" : message;
  } catch (...) {
    g_error.clear();
  }
}

/* NULL is a no-op. A second destroy records an error and does not free. */
bool claim(void *ptr) noexcept {
  if (ptr == nullptr) {
    return false;
  }
  try {
    std::lock_guard<std::mutex> lock(g_live_mu);
    if (g_live.erase(ptr) == 0) {
      set_error("double destroy");
      return false;
    }
  } catch (...) {
    set_error("handle registry is unavailable");
    return false;
  }
  return true;
}

/* Calls on a mutable handle own it, so the slot is cleared on entry and
 * seams_handle_error reports this call. */
struct handle_guard {
  const seams_object *previous;
  explicit handle_guard(seams_object *object) : previous(g_current) {
    g_current = object;
    if (object != nullptr) {
      std::lock_guard<std::mutex> lock(object->error_mu);
      object->error.clear();
    }
  }
  ~handle_guard() { g_current = previous; }
  handle_guard(const handle_guard &) = delete;
  handle_guard &operator=(const handle_guard &) = delete;
};

/* Calls on a const handle write nothing on success. A failure still lands
 * in the slot, so bindings that read seams_handle_error after a failed
 * status keep their message. */
struct read_guard {
  const seams_object *previous;
  explicit read_guard(const seams_object *object) : previous(g_current) {
    g_current = object;
  }
  ~read_guard() { g_current = previous; }
  read_guard(const read_guard &) = delete;
  read_guard &operator=(const read_guard &) = delete;
};

void note_error() noexcept {
  if (g_current == nullptr) {
    return;
  }
  try {
    std::lock_guard<std::mutex> lock(g_current->error_mu);
    try {
      g_current->error = g_error;
    } catch (...) {
      g_current->error.clear();
    }
  } catch (...) {
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

void clear_error() noexcept { g_error.clear(); }

seams_status fail(seams_status code, const char *message) noexcept {
  set_error(message);
  note_error();
  return code;
}

/* Call only from inside a catch handler. */
void remember() noexcept {
  try {
    throw;
  } catch (const std::bad_alloc &) {
    set_error("out of memory");
  } catch (const std::exception &ex) {
    set_error(ex.what());
  } catch (...) {
    set_error("unknown error");
  }
  note_error();
}

seams::FairWidth use_threads() {
#ifdef SEAMS_HAS_OPENMP
  if (g_threads > 0) {
    omp_set_num_threads(g_threads);
  }
#endif
  return seams::FairWidth(true);
}

bool finite_box(const double *box9) {
  for (int i = 0; i < 9; ++i) {
    if (!std::isfinite(box9[i])) {
      return false;
    }
  }
  return true;
}

// A double is NaN or infinite exactly when its exponent bits are all set,
// and adding one to an all-ones exponent carries into the sign bit, so the
// OR of the sums has its top bit set exactly when some value is not finite.
// Add, AND and OR vectorise at any width, where a 64-bit compare needs
// SSE4.1 and std::isfinite does not vectorise. The scan runs from the end,
// which leaves the start of the input in cache for a write pass that walks
// forward. The explicit prefetch 8 KiB ahead cut the reload of a 21 MB
// input by about 4%.
template <std::size_t Lanes> inline bool finite_scan(const double *v, std::size_t count) {
  constexpr std::uint64_t exponent = 0x7ff0000000000000ULL;
  constexpr std::uint64_t carry = 0x0010000000000000ULL;
  constexpr std::size_t ahead = 1024;
  std::uint64_t acc[Lanes] = {};
  std::size_t k = count;
  while (k >= Lanes) {
    k -= Lanes;
#if defined(__GNUC__) || defined(__clang__)
    if (k > ahead) {
      __builtin_prefetch(v + (k - ahead));
    }
#endif
    std::uint64_t bits[Lanes];
    std::memcpy(bits, v + k, sizeof bits);
    for (std::size_t l = 0; l < Lanes; ++l) {
      acc[l] |= (bits[l] & exponent) + carry;
    }
  }
  std::uint64_t bad = 0;
  for (std::size_t l = 0; l < Lanes; ++l) {
    bad |= acc[l];
  }
  while (k > 0) {
    --k;
    std::uint64_t bits;
    std::memcpy(&bits, v + k, sizeof bits);
    bad |= (bits & exponent) + carry;
  }
  return (bad >> 63) == 0;
}

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
__attribute__((target("avx512f"))) bool finite_scan_avx512(const double *v, std::size_t count) {
  return finite_scan<32>(v, count);
}

__attribute__((target("avx2"))) bool finite_scan_avx2(const double *v, std::size_t count) {
  return finite_scan<16>(v, count);
}
#endif

bool finite_coordinates(const double *xyz, std::size_t n_atoms) {
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
  // __builtin_cpu_supports needs libgcc's __cpu_model, which the macOS
  // and Windows wheels do not link. The CPUID helper is the same check.
  const auto &cpu = seams::cpu::features();
  if (cpu.avx512f) {
    return finite_scan_avx512(xyz, 3 * n_atoms);
  }
  if (cpu.avx2) {
    return finite_scan_avx2(xyz, 3 * n_atoms);
  }
#endif
  return finite_scan<8>(xyz, 3 * n_atoms);
}

// The mutual k-nearest fast path searches `typ` on the frame as held.
// With `type <= 0` that matches `unified_copy` only when every atom has
// one positive type; a mixed frame takes the copied path instead.
bool mutual_fast_type(const Cloud &cloud, int type, int &typ) {
  if (type > 0) {
    typ = type;
    return true;
  }
  if (cloud.nop <= 0) {
    return false;
  }
  const int first = cloud.pts[0].type;
  if (first <= 0) {
    return false;
  }
  for (int i = 1; i < cloud.nop; ++i) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != first) {
      return false;
    }
  }
  typ = first;
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
  // cloud.idIndexMap holds atom id i + 1 -> i for this many atoms, so a
  // reload with default or sequential ids keeps the map. Zero after any
  // other load.
  std::size_t identity_ids = 0;
  // A load is built here first and swapped in only once it has passed
  // every check, so a failed load leaves the frame as it was.
  std::unordered_map<int, int> spare_ids;
  std::vector<double> spare_box;
  std::vector<double> spare_low;
  // With stage set, set_positions writes the new coordinates here in the
  // same pass that checks them, and swaps this with cloud.pts only when
  // every value was finite. While staged_current is set, types and ids
  // here match the frame; any load clears it.
  bool stage = false;
  decltype(Cloud::pts) staged;
  bool staged_current = false;
};

struct seams_neighbours : seams_object {
  std::vector<int> offsets;
  std::vector<int> cols;
  int n = 0;
};

struct seams_rings : seams_object {
  std::vector<std::vector<int>> rings;
  // Filled lazily by const label reads, under refine_mu.
  mutable std::vector<int> symmetry;
  int n_atoms = 0;
  // Kept so prism and cage labels can be filled when a caller asks for
  // them. The ring census itself does not need that pass. The prism
  // search takes the frame by non-const reference, so it is mutable too
  // and only touched under refine_mu.
  mutable Cloud cloud;
  std::vector<int> offsets;
  std::vector<int> cols;
  mutable std::vector<std::vector<int>> idx;
  mutable bool refined = false;
  mutable std::mutex refine_mu;
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

struct seams_msd : seams_object {
  traj::MsdAccumulator acc;
  seams_msd(std::size_t natoms, std::size_t max_lag, std::size_t origin_stride,
            double dt, const double axis[3])
      : acc(natoms, max_lag, origin_stride, dt, axis) {}
};

struct seams_hbond_tcf : seams_object {
  traj::HbondLifetimeAccumulator acc;
  seams_hbond_tcf(std::size_t max_lag, std::size_t origin_stride, double dt,
                  bool directed)
      : acc(max_lag, origin_stride, dt, directed) {}
};

struct seams_reorient : seams_object {
  traj::ReorientationAccumulator acc;
  seams_reorient(std::size_t n_vectors, std::size_t max_lag,
                 std::size_t origin_stride, double dt)
      : acc(n_vectors, max_lag, origin_stride, dt) {}
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

struct hbond_edge_internal {
  int donor = -1;
  int acceptor = -1;
  int hydrogen = -1;
  double d_DA = NAN;
  double angle = NAN;
};

std::vector<hbond_edge_internal>
hbond_edges_detailed(const Cloud &cloud, double cutoff, int oxygen_type, int hydrogen_type,
                     double dist_cutoff, double angle_cutoff, bool with_h,
                     bool with_geometry) {
  Cloud heavy;
  Cloud hydro;
  heavy.box = cloud.box;
  heavy.boxLow = cloud.boxLow;
  hydro.box = cloud.box;
  hydro.boxLow = cloud.boxLow;
  std::vector<int> heavy_to_frame;
  std::vector<int> hydro_to_frame;
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &point = cloud.pts[static_cast<std::size_t>(i)];
    if (point.type == hydrogen_type) {
      hydro.pts.push_back(point);
      hydro_to_frame.push_back(i);
    } else if (oxygen_type <= 0 || point.type == oxygen_type) {
      heavy.pts.push_back(point);
      if (oxygen_type <= 0) {
        heavy.pts.back().type = 1;
      }
      heavy_to_frame.push_back(i);
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
  std::vector<hbond_edge_internal> out;
  if (heavy.nop == 0 || hydro.nop == 0) {
    return out;
  }
  const int typ = oxygen_type > 0 ? oxygen_type : 1;
  auto rows = nneigh::neighListO(cutoff, heavy, typ);
  auto molIDlist = molSys::hAtomMolList(hydro, heavy);
  std::unordered_map<int, std::vector<int>> donor_h_by_mol;
  for (const auto &row : molIDlist) {
    if (row.size() < 3) {
      continue;
    }
    donor_h_by_mol[row[0]].push_back(row[1]);
    donor_h_by_mol[row[0]].push_back(row[2]);
  }
  for (const auto &row : rows) {
    if (row.size() < 2) {
      continue;
    }
    const auto it_i = heavy.idIndexMap.find(row[0]);
    if (it_i == heavy.idIndexMap.end()) {
      continue;
    }
    const int iatom = it_i->second;
    for (std::size_t j = 1; j < row.size(); ++j) {
      const auto it_j = heavy.idIndexMap.find(row[j]);
      if (it_j == heavy.idIndexMap.end()) {
        continue;
      }
      const int jatom = it_j->second;
      if (iatom >= jatom) {
        continue;
      }
      const int donors[2] = {iatom, jatom};
      const int acceptors[2] = {jatom, iatom};
      bool accepted = false;
      for (int d = 0; d < 2 && !accepted; ++d) {
        const int donor = donors[d];
        const int acceptor = acceptors[d];
        const auto donor_it = donor_h_by_mol.find(heavy.pts[static_cast<std::size_t>(donor)].molID);
        if (donor_it == donor_h_by_mol.end()) {
          continue;
        }
        const auto oo = gen::relDist(heavy, acceptor, donor);
        const std::vector<double> oo_vec{oo[0], oo[1], oo[2]};
        const double oo_len = std::sqrt(oo[0] * oo[0] + oo[1] * oo[1] + oo[2] * oo[2]);
        const double dist2 = dist_cutoff * dist_cutoff;
        for (const int hindex : donor_it->second) {
          if (hindex < 0 || hindex >= hydro.nop) {
            continue;
          }
          const auto oh = gen::relDistFromPoint(
              heavy, acceptor, hydro.pts[static_cast<std::size_t>(hindex)].x,
              hydro.pts[static_cast<std::size_t>(hindex)].y,
              hydro.pts[static_cast<std::size_t>(hindex)].z);
          const double oh2 = oh[0] * oh[0] + oh[1] * oh[1] + oh[2] * oh[2];
          if (oh2 >= dist2) {
            continue;
          }
          const std::vector<double> oh_vec{oh[0], oh[1], oh[2]};
          const double angle = gen::radDeg(gen::eigenVecAngle(oo_vec, oh_vec));
          if (angle > angle_cutoff) {
            continue;
          }
          hbond_edge_internal edge;
          edge.donor = heavy_to_frame[static_cast<std::size_t>(donor)];
          edge.acceptor = heavy_to_frame[static_cast<std::size_t>(acceptor)];
          edge.hydrogen = with_h ? hydro_to_frame[static_cast<std::size_t>(hindex)] : -1;
          if (with_geometry) {
            edge.d_DA = oo_len;
            edge.angle = angle;
          }
          out.push_back(edge);
          accepted = true;
          break;
        }
      }
    }
  }
  std::sort(out.begin(), out.end(), [](const hbond_edge_internal &a, const hbond_edge_internal &b) {
    if (a.donor != b.donor) {
      return a.donor < b.donor;
    }
    if (a.acceptor != b.acceptor) {
      return a.acceptor < b.acceptor;
    }
    return a.hydrogen < b.hydrogen;
  });
  return out;
}

} // namespace

/* mols NULL stores molID = atomID, the historical set_atoms behaviour. */
static seams_status load_frame(seams_system *sys, const double *xyz, size_t n,
                               const double *box9, size_t box_len, const int *types,
                               const int *ids, const int *mols) {
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
    if (!nneigh::setBoxFromLammpsEdges(sys->spare_box, sys->spare_low, box9[0], box9[1],
                                      box9[2], box9[3], box9[4], box9[5], box9[6],
                                      box9[7], box9[8])) {
      return fail(SEAMS_ERR_BOX, "box edge lengths must be positive");
    }
    if (!finite_coordinates(xyz, n)) {
      return fail(SEAMS_ERR_ARG, "a coordinate is not finite");
    }
    bool sequential = true;
    if (ids != nullptr) {
      for (int i = 0; i < count && sequential; ++i) {
        sequential = ids[i] == i + 1;
      }
    }
    const bool rebuild = !sequential || !sys->ready || sys->identity_ids != n;
    if (rebuild) {
      auto &map = sys->spare_ids;
      map.clear();
      map.reserve(n);
      for (int i = 0; i < count; ++i) {
        if (!map.emplace(ids == nullptr ? i + 1 : ids[i], i).second) {
          return fail(SEAMS_ERR_ARG, "atom ids are not unique");
        }
      }
    }
    Cloud &cloud = sys->cloud;
    sys->staged_current = false;
    cloud.pts.resize(n);
    // Nothing below allocates, so the frame is either the old one or the
    // whole new one.
    for (int i = 0; i < count; ++i) {
      const std::size_t base = static_cast<std::size_t>(i) * 3;
      auto &point = cloud.pts[static_cast<std::size_t>(i)];
      point.x = xyz[base];
      point.y = xyz[base + 1];
      point.z = xyz[base + 2];
      point.type = types == nullptr ? 1 : types[i];
      point.atomID = ids == nullptr ? i + 1 : ids[i];
      point.molID = mols == nullptr ? point.atomID : mols[i];
      point.iceType = molSys::atom_state_type::unclassified;
      point.inSlice = true;
      point.c_ij.clear();
    }
    if (rebuild) {
      cloud.idIndexMap.swap(sys->spare_ids);
    }
    cloud.box.swap(sys->spare_box);
    cloud.boxLow.swap(sys->spare_low);
    cloud.nop = count;
    cloud.currentFrame = 1;
    sys->identity_ids = sequential ? n : 0;
    sys->ready = true;
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return g_error == "out of memory" ? SEAMS_ERR_NOMEM : SEAMS_ERR_INTERNAL;
  }
}

extern "C" {

int seams_abi(void) { return SEAMS_ABI; }
int seams_version(void) { return SEAMS_VERSION_NUMBER; }
int seams_version_major(void) { return SEAMS_VERSION_MAJOR; }
int seams_version_minor(void) { return SEAMS_VERSION_MINOR; }
int seams_version_patch(void) { return SEAMS_VERSION_PATCH; }
const char *seams_version_string(void) { return SEAMS_VERSION_STRING; }

const char *seams_status_string(seams_status status) {
  // A C caller may pass any int. In C++ a value past the enumerators is not
  // a seams_status, so the switch reads the bits as an int.
  static_assert(sizeof(seams_status) == sizeof(int), "seams_status is an int");
  int code = 0;
  std::memcpy(&code, &status, sizeof code);
  switch (code) {
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
    auto *sys = new seams_system();
    track(sys);
    return sys;
  } catch (...) {
    set_error("out of memory");
    return nullptr;
  }
}

void seams_system_destroy(seams_system *sys) {
  if (!claim(sys)) {
    return;
  }
  delete sys;
}

seams_status seams_system_set_atoms(seams_system *sys, const double *xyz, size_t n,
                                    const double *box9, size_t box_len, const int *types,
                                    const int *ids) {
  handle_guard guard(sys);
  clear_error();
  return load_frame(sys, xyz, n, box9, box_len, types, ids, nullptr);
}

seams_status seams_system_set_hydrogen_frame(
    seams_system *sys, const double *heavy_xyz, size_t n_heavy, const int *heavy_types,
    const int *heavy_ids, const int *heavy_mols, const double *hydrogen_xyz,
    size_t n_hydrogen, const int *hydrogen_types, const int *hydrogen_ids,
    const int *hydrogen_mols, const double *box9, size_t box_len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || heavy_xyz == nullptr || box9 == nullptr) {
    return fail(SEAMS_ERR_NULL, "system, heavy coordinates, or box is null");
  }
  if (n_heavy == 0) {
    return fail(SEAMS_ERR_RANGE, "frame has no heavy atoms");
  }
  if (heavy_types == nullptr) {
    return fail(SEAMS_ERR_NULL, "heavy types are null");
  }
  if (n_hydrogen > 0 && (hydrogen_xyz == nullptr || hydrogen_types == nullptr)) {
    return fail(SEAMS_ERR_NULL, "hydrogen coordinates or types are null");
  }
  if (n_heavy > static_cast<std::size_t>(INT_MAX) ||
      n_hydrogen > static_cast<std::size_t>(INT_MAX) - n_heavy) {
    return fail(SEAMS_ERR_RANGE, "atom count does not fit in the frame");
  }
  try {
    const std::size_t n = n_heavy + n_hydrogen;
    std::vector<double> xyz(n * 3);
    std::vector<int> types(n);
    std::vector<int> ids(n);
    std::vector<int> mols(n);
    for (std::size_t i = 0; i < n_heavy; ++i) {
      const std::size_t base = i * 3;
      xyz[base] = heavy_xyz[base];
      xyz[base + 1] = heavy_xyz[base + 1];
      xyz[base + 2] = heavy_xyz[base + 2];
      types[i] = heavy_types[i];
      ids[i] = heavy_ids == nullptr ? static_cast<int>(i + 1) : heavy_ids[i];
      mols[i] = heavy_mols == nullptr ? ids[i] : heavy_mols[i];
    }
    for (std::size_t i = 0; i < n_hydrogen; ++i) {
      const std::size_t k = n_heavy + i;
      const std::size_t src = i * 3;
      const std::size_t base = k * 3;
      xyz[base] = hydrogen_xyz[src];
      xyz[base + 1] = hydrogen_xyz[src + 1];
      xyz[base + 2] = hydrogen_xyz[src + 2];
      types[k] = hydrogen_types[i];
      ids[k] = hydrogen_ids == nullptr ? static_cast<int>(k + 1) : hydrogen_ids[i];
      mols[k] = hydrogen_mols == nullptr ? ids[k] : hydrogen_mols[i];
    }
    return load_frame(sys, xyz.data(), n, box9, box_len, types.data(), ids.data(), mols.data());
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
    const auto fair = use_threads();
    int typ = 1;
    if (graph == SEAMS_GRAPH_CUTOFF) {
      // The cutoff list only reads coordinates. Copying the frame first
      // duplicated every point, including its bond vectors, and that copy
      // was slower than the search.
      const Cloud &work = sys->cloud;
      const std::vector<int> *sub = nullptr;
      std::vector<int> subset;
      if (type > 0) {
        typ = type;
        bool uniform = work.nop > 0;
        for (int i = 0; i < work.nop && uniform; ++i) {
          if (work.pts[static_cast<std::size_t>(i)].type != typ) {
            uniform = false;
          }
        }
        if (!uniform) {
          subset.reserve(static_cast<std::size_t>(work.nop));
          for (int i = 0; i < work.nop; ++i) {
            if (work.pts[static_cast<std::size_t>(i)].type == typ) {
              subset.push_back(i);
            }
          }
          sub = &subset;
        }
      }
      auto direct = std::make_unique<seams_neighbours>();
      direct->n = work.nop;
      if (nneigh::orthoFullCsr(work, sub, cutoff, direct->offsets, direct->cols, false)) {
        *out = static_cast<seams_neighbours *>(track(direct.release()));
        return SEAMS_OK;
      }
    }
    Cloud work = unified_copy(sys->cloud, type, typ);
    auto rows = bond_rows(work, cutoff, typ, graph);
    auto idx = nneigh::neighbourListByIndex(work, rows);
    auto nb = std::make_unique<seams_neighbours>();
    const int n = work.nop;
    if (static_cast<int>(idx.size()) < n) {
      idx.resize(static_cast<std::size_t>(n));
    }
    nb->n = n;
    std::vector<int> degrees(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
      const auto &row = idx[static_cast<std::size_t>(i)];
      int degree = 0;
      for (std::size_t k = 0; k < row.size(); ++k) {
        if (row[k] != i) {
          ++degree;
        }
      }
      degrees[static_cast<std::size_t>(i)] = degree;
    }
    if (!seams::prefixCounts(degrees.data(), n, nb->offsets)) {
      return fail(SEAMS_ERR_RANGE, "neighbour list does not fit");
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
    *out = static_cast<seams_neighbours *>(track(nb.release()));
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

seams_status seams_rings_chordless(const seams_neighbours *nb, int size, int *rings,
                                   size_t rings_len, size_t *n_rings, int *per_atom,
                                   size_t per_atom_len) {
  handle_guard guard(const_cast<seams_neighbours *>(nb));
  clear_error();
  if (nb == nullptr) {
    return fail(SEAMS_ERR_NULL, "neighbours handle is null");
  }
  if (rings != nullptr && n_rings == nullptr) {
    return fail(SEAMS_ERR_NULL, "rings needs n_rings");
  }
  if (rings == nullptr && n_rings == nullptr && per_atom == nullptr) {
    return fail(SEAMS_ERR_NULL, "rings, n_rings and per_atom are all null");
  }
  if (size < 3 || size > chordless::kMaxSize) {
    return fail(SEAMS_ERR_RANGE, "size must be from 3 to 12");
  }
  const std::size_t n = nb->offsets.empty() ? 0 : nb->offsets.size() - 1;
  if (per_atom != nullptr && per_atom_len < n) {
    return fail(SEAMS_ERR_RANGE, "per_atom is short");
  }
  try {
    const auto fair = use_threads();
    std::vector<int> found;
    const std::size_t count = chordless::cycles(nb->offsets, nb->cols, size,
                                                rings != nullptr ? &found : nullptr, per_atom);
    if (n_rings != nullptr) {
      *n_rings = count;
    }
    if (rings != nullptr) {
      if (rings_len < found.size()) {
        return fail(SEAMS_ERR_RANGE, "rings buffer is short");
      }
      if (!found.empty()) {
        std::memcpy(rings, found.data(), found.size() * sizeof(int));
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

void seams_neighbours_destroy(seams_neighbours *nb) {
  if (!claim(nb)) {
    return;
  }
  delete nb;
}

int seams_neighbours_count(const seams_neighbours *nb) {
  return nb == nullptr ? 0 : nb->n;
}

seams_status seams_neighbours_csr(const seams_neighbours *nb, const int **offsets,
                                  int *n_offsets, const int **cols, int *nnz) {
  read_guard guard(nb);
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
    const auto fair = use_threads();
    int typ = 1;
    if ((rule == SEAMS_CHILL_PUBLISHED || rule == SEAMS_CHILL_OVITO) &&
        (graph == SEAMS_GRAPH_CUTOFF ||
         (graph == SEAMS_GRAPH_MUTUAL && mutual_fast_type(sys->cloud, type, typ)))) {
      // Labels are written onto the frame. The bond vectors are not: the
      // caller reads the label buffer. Skipping the frame copy keeps the
      // same iceType values the copied frame used to carry out.
      Cloud &work = sys->cloud;
      if (type > 0) {
        typ = type;
      }
      std::vector<int> subset;
      const std::vector<int> *sub = nullptr;
      if (graph == SEAMS_GRAPH_CUTOFF && type > 0) {
        bool uniform = work.nop > 0;
        for (int i = 0; i < work.nop && uniform; ++i) {
          if (work.pts[static_cast<std::size_t>(i)].type != typ) {
            uniform = false;
          }
        }
        if (!uniform) {
          subset.reserve(static_cast<std::size_t>(work.nop));
          for (int i = 0; i < work.nop; ++i) {
            if (work.pts[static_cast<std::size_t>(i)].type == typ) {
              subset.push_back(i);
            }
          }
          sub = &subset;
        }
      }
      std::vector<int> offsets;
      std::vector<int> cols;
      bool listed = false;
      if (graph == SEAMS_GRAPH_MUTUAL) {
        nneigh::kNearestCsr(work, 4, cutoff + 1.5, typ, true, offsets, cols, false);
        listed = static_cast<int>(offsets.size()) == work.nop + 1;
      } else {
        listed = nneigh::orthoFullCsr(work, sub, cutoff, offsets, cols, false);
      }
      if (listed) {
        if (rule == SEAMS_CHILL_OVITO) {
          chill::chillPlusOvitoFromCsr(work, offsets, cols);
        } else {
          chill::chillPlusFromCsr(work, offsets, cols, false);
        }
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
    Cloud work = unified_copy(sys->cloud, type, typ);
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
    const auto fair = use_threads();
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
    const auto fair = use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    std::vector<std::vector<int>> idx;
    std::vector<std::vector<int>> found;
    std::vector<int> offsets;
    std::vector<int> cols;
    bool listed = false;
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
      if (nneigh::orthoFullCsr(work, sub, cutoff, offsets, cols, false)) {
        listed = true;
        found = primitive::ringNetwork(offsets, cols, max_depth);
      }
    }
    if (!listed) {
      auto rows = bond_rows(work, cutoff, typ, graph);
      idx = nneigh::neighbourListByIndex(work, rows);
      if (static_cast<int>(idx.size()) < work.nop) {
        idx.resize(static_cast<std::size_t>(work.nop));
      }
      found = primitive::ringNetwork(idx, max_depth);
    }
    auto rings = std::make_unique<seams_rings>();
    rings->n_atoms = work.nop;
    rings->cloud = std::move(work);
    rings->offsets = std::move(offsets);
    rings->cols = std::move(cols);
    rings->idx = std::move(idx);
    rings->rings = std::move(found);
    rings->symmetry.assign(rings->rings.size(), SEAMS_RING_OTHER);
    for (std::size_t i = 0; i < rings->rings.size(); ++i) {
      const int size = static_cast<int>(rings->rings[i].size());
      if (size == 5) {
        rings->symmetry[i] = SEAMS_RING_PENTAGON;
      } else if (size == 6) {
        rings->symmetry[i] = SEAMS_RING_HEXAGON;
      }
    }
    *out = static_cast<seams_rings *>(track(rings.release()));
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

// The census is the ring list. Cage and prism codes are filled the first
// time a caller reads symmetry, from the cloud and neighbour rows kept
// on the handle.
void refine_ring_labels(const seams_rings *rings) {
  if (rings == nullptr) {
    return;
  }
  std::lock_guard<std::mutex> lock(rings->refine_mu);
  if (rings->refined) {
    return;
  }
  if (rings->idx.empty() && rings->offsets.size() >= 2) {
    const int nop = rings->n_atoms;
    rings->idx.assign(static_cast<std::size_t>(std::max(nop, 0)), {});
    for (int i = 0; i < nop && i + 1 < static_cast<int>(rings->offsets.size()); ++i) {
      const int a = rings->offsets[static_cast<std::size_t>(i)];
      const int b = rings->offsets[static_cast<std::size_t>(i) + 1];
      // cageAffiliation and findPrisms read a row from its second entry,
      // as neighbourListByIndex writes it.
      auto &row = rings->idx[static_cast<std::size_t>(i)];
      row.reserve(static_cast<std::size_t>(std::max(0, b - a)) + 1);
      row.push_back(i);
      for (int k = a; k < b; ++k) {
        row.push_back(rings->cols[static_cast<std::size_t>(k)]);
      }
    }
  }
  const auto &idx = rings->idx;
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
        std::vector<double> rmsd(static_cast<std::size_t>(std::max(rings->cloud.nop, 0)), -1.0);
        int perfect = 0;
        int imperfect = 0;
        ring::findPrisms(bySize[size], kind, perfect, imperfect, idx, rings->cloud, rmsd,
                         false);
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
  rings->refined = true;
}

void seams_rings_destroy(seams_rings *rings) {
  if (!claim(rings)) {
    return;
  }
  delete rings;
}

int seams_rings_count(const seams_rings *rings) {
  return rings == nullptr ? 0 : static_cast<int>(rings->rings.size());
}

seams_status seams_rings_histogram(const seams_rings *rings, int *counts,
                                   int counts_len) {
  read_guard guard(rings);
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
  read_guard guard(rings);
  clear_error();
  if (rings == nullptr) {
    return fail(SEAMS_ERR_NULL, "rings handle is null");
  }
  if (sizes == nullptr && symmetry == nullptr) {
    return fail(SEAMS_ERR_NULL, "sizes and symmetry are both null");
  }
  if (symmetry != nullptr) {
    try {
      refine_ring_labels(rings);
    } catch (const std::bad_alloc &) {
      return fail(SEAMS_ERR_NOMEM, "out of memory");
    } catch (...) {
      remember();
      return SEAMS_ERR_INTERNAL;
    }
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
  read_guard guard(rings);
  clear_error();
  if (rings == nullptr || fraction == nullptr) {
    return fail(SEAMS_ERR_NULL, "fraction output is null");
  }
  if (rings->n_atoms <= 0) {
    *fraction = 0.0;
    return SEAMS_OK;
  }
  std::vector<char> hit;
  try {
    hit.assign(static_cast<std::size_t>(rings->n_atoms), 0);
  } catch (...) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  }
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
  read_guard guard(nb);
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
  read_guard guard(rings);
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
    if (!gen::file_exists(path)) {
      return fail(SEAMS_ERR_ARG, "cannot read path");
    }
    Cloud cloud = sinp::openFrame(path, frame, type_filter);
    if (cloud.nop <= 0) {
      return fail(SEAMS_ERR_RANGE, "frame has no atoms");
    }
    if (cloud.currentFrame <= 0) {
      cloud.currentFrame = frame;
    }
    sys->cloud = std::move(cloud);
    sys->identity_ids = 0;
    sys->staged_current = false;
    sys->ready = true;
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (const std::exception &ex) {
    const char *what = ex.what();
    const bool empty = what == nullptr || what[0] == '\0';
    const bool past = !empty && std::strstr(what, "past the last") != nullptr;
    return fail(past ? SEAMS_ERR_RANGE : SEAMS_ERR_ARG, empty ? "cannot read frame" : what);
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
    const auto fair = use_threads();
    int typ = 1;
    if (rule == SEAMS_CHILL_PUBLISHED &&
        (graph == SEAMS_GRAPH_CUTOFF ||
         (graph == SEAMS_GRAPH_MUTUAL && mutual_fast_type(sys->cloud, type, typ)))) {
      // Correlations stay in a flat table. The per-atom bond vectors are
      // not built. iceType is written on the frame the caller already holds.
      Cloud &work = sys->cloud;
      if (type > 0) {
        typ = type;
      }
      std::vector<int> offsets;
      std::vector<int> cols;
      bool listed = false;
      if (graph == SEAMS_GRAPH_CUTOFF) {
        const std::vector<int> *sub = nullptr;
        std::vector<int> subset;
        if (type > 0) {
          bool uniform = work.nop > 0;
          for (int i = 0; i < work.nop && uniform; ++i) {
            if (work.pts[static_cast<std::size_t>(i)].type != typ) {
              uniform = false;
            }
          }
          if (!uniform) {
            subset.reserve(static_cast<std::size_t>(work.nop));
            for (int i = 0; i < work.nop; ++i) {
              if (work.pts[static_cast<std::size_t>(i)].type == typ) {
                subset.push_back(i);
              }
            }
            sub = &subset;
          }
        }
        listed = nneigh::orthoFullCsr(work, sub, cutoff, offsets, cols, false);
      } else {
        nneigh::kNearestCsr(work, 4, cutoff + 1.5, typ, true, offsets, cols, false);
        listed = static_cast<int>(offsets.size()) == work.nop + 1;
      }
      if (listed) {
        chill::chillFromCsr(work, offsets, cols);
        write_ice_labels(work, type, labels, counts9);
        return SEAMS_OK;
      }
    }
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
  if (sys->cloud.box.size() >= 3) {
    const double edge =
        std::min(sys->cloud.box[0], std::min(sys->cloud.box[1], sys->cloud.box[2]));
    double absurd = edge;
    if (sys->cloud.box.size() >= 6) {
      absurd = std::max(absurd, std::max(std::abs(sys->cloud.box[3]),
                                         std::max(std::abs(sys->cloud.box[4]),
                                                  std::abs(sys->cloud.box[5]))));
      absurd = std::max(absurd, std::max(sys->cloud.box[0], std::max(sys->cloud.box[1],
                                                                    sys->cloud.box[2])));
    }
    double shear = 0.0;
    if (sys->cloud.box.size() >= 6 && edge > 0.0) {
      shear = std::max(std::abs(sys->cloud.box[3]),
                       std::max(std::abs(sys->cloud.box[4]), std::abs(sys->cloud.box[5]))) /
              edge;
    }
    const bool many = edge > 0.0 && (cutoff / edge > 64.0 || shear > 64.0);
    if (many || absurd > 1.0e8) {
      return fail(SEAMS_ERR_BOX, "cutoff spans too many periodic images");
    }
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
    const auto fair = use_threads();
    int typ = 1;
    // A positive type already names the atoms. Copying the cloud also
    // copies the atom-id map, which is the bulk of a q6 call at 1e5
    // molecules and does not change the harmonics.
    Cloud owned;
    const Cloud *view = &sys->cloud;
    if (type > 0) {
      typ = type;
    } else {
      owned = unified_copy(sys->cloud, type, typ);
      view = &owned;
    }
    const Cloud &work = *view;
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
    const auto fair = use_threads();
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

// Cutoff graph as the CSR the ring search already uses, plus the index
// rows the cage walk wants (self first, then neighbours). Skips the
// atom-ID round trip in neighListO. Returns false when the box cannot
// use that list or the frame mixes types.
bool cutoff_csr_graph(const Cloud &work, double cutoff, int typ,
                      std::vector<std::vector<int>> &idx,
                      std::vector<std::vector<int>> &rings, int depth) {
  if (work.nop <= 0) {
    idx.clear();
    rings.clear();
    return true;
  }
  // The tilted cell list does not feed the same pairs as neighListO.
  // Keep that builder for triclinic frames.
  if (work.box.size() != 3) {
    return false;
  }
  for (int i = 0; i < work.nop; ++i) {
    if (work.pts[static_cast<std::size_t>(i)].type != typ) {
      return false;
    }
  }
  std::vector<int> offsets;
  std::vector<int> cols;
  if (!nneigh::orthoFullCsr(work, nullptr, cutoff, offsets, cols, false)) {
    return false;
  }
  idx.assign(static_cast<std::size_t>(work.nop), {});
  const int rows = static_cast<int>(offsets.size()) - 1;
  for (int i = 0; i < work.nop && i < rows; ++i) {
    const int a = offsets[static_cast<std::size_t>(i)];
    const int b = offsets[static_cast<std::size_t>(i + 1)];
    auto &row = idx[static_cast<std::size_t>(i)];
    row.reserve(static_cast<std::size_t>(std::max(0, b - a) + 1));
    row.push_back(i);
    for (int k = a; k < b; ++k) {
      row.push_back(cols[static_cast<std::size_t>(k)]);
    }
  }
  // One thread is faster on the packed rows. A team is faster on the
  // per-atom vectors, which is the layout the command uses.
#ifdef SEAMS_HAS_OPENMP
  const bool packed = omp_in_parallel() || omp_get_max_threads() <= 1;
#else
  const bool packed = true;
#endif
  rings = packed ? primitive::ringNetwork(offsets, cols, depth)
                 : primitive::ringNetwork(idx, depth);
  return true;
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
    const auto fair = use_threads();
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
      std::vector<std::vector<int>> idx;
      std::vector<std::vector<int>> found_rings;
      const bool csr = graph == SEAMS_GRAPH_CUTOFF &&
                       cutoff_csr_graph(work, cutoff, typ, idx, found_rings, 6);
      if (!csr) {
        auto rows = bond_rows(work, cutoff, typ, graph);
        idx = nneigh::neighbourListByIndex(work, rows);
      }
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
        if (!csr) {
          found_rings = primitive::ringNetwork(idx, 6);
        }
        auto six = six_of(ring::contractibleRings(found_rings, work));
        flags_from_affiliation(n, six, ring::cageAffiliation(six, idx), hc, ddc);
      }
    }
    auto cages = std::make_unique<seams_cage_set>();
    tally_cages(cages.get(), work, type, hc, ddc);
    *out = static_cast<seams_cage_set *>(track(cages.release()));
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
  if (!claim(cages)) {
    return;
  }
  delete cages;
}

seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex, int *n_cubic,
                                int *n_water, int *n_both) {
  read_guard guard(cages);
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
  read_guard guard(cages);
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

static seams_status signature_with_edge(seams_system *sys, double cutoff, int type, int graph,
                                        const char *spec, const char *edge, size_t *n_cages,
                                        size_t *n_atoms) {
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
    const auto fair = use_threads();
    // A type filter that keeps every atom is the frame itself. Copying the
    // cloud copies its atom-id map, and that copy is what made the library
    // call slower than the command on a wide run.
    int typ = type > 0 ? type : 1;
    Cloud copied;
    const Cloud *held = &sys->cloud;
    if (type > 0) {
      copied = unified_copy(sys->cloud, type, typ);
      held = &copied;
    }
    const Cloud &work = *held;
    const auto signature = cage::Signature::parse(spec);
    const int depth = std::max(signature.maxRingSize(), 3);
    std::vector<std::vector<int>> idx;
    std::vector<std::vector<int>> found_rings;
    const bool csr = graph == SEAMS_GRAPH_CUTOFF &&
                     cutoff_csr_graph(work, cutoff, typ, idx, found_rings, depth);
    if (!csr) {
      auto rows = bond_rows(work, cutoff, typ, graph);
      idx = nneigh::neighbourListByIndex(work, rows);
      found_rings = primitive::ringNetwork(idx, depth);
    }
    const auto rings = ring::cageRings(found_rings, work, cutoff);
    const auto window = cage::parseEdgeWindow(edge == nullptr ? "" : edge);
    std::vector<std::array<double, 3>> xyz(static_cast<std::size_t>(std::max(work.nop, 0)));
    for (int i = 0; i < work.nop; ++i) {
      const auto &point = work.pts[static_cast<std::size_t>(i)];
      xyz[static_cast<std::size_t>(i)] = {point.x, point.y, point.z};
    }
    const auto found_of = cage::findBySignatures(rings, idx, std::vector<cage::Signature>{signature},
                                                 {}, &xyz, &work.box, window, &work.boxLow);
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
    const auto fair = use_threads();
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
    const auto fair = use_threads();
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
    *out = static_cast<seams_brine *>(track(front.release()));
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
  if (!claim(front)) {
    return;
  }
  delete front;
}

seams_status seams_brine_fronts(const seams_brine *front, double *z_enter, double *z_exit,
                                double *thickness, int *found) {
  read_guard guard(front);
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
  read_guard guard(front);
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
  read_guard guard(front);
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

seams_status seams_cages_signature(seams_system *sys, double cutoff, int type, int graph,
                                   const char *spec, size_t *n_cages, size_t *n_atoms) {
  return signature_with_edge(sys, cutoff, type, graph, spec, nullptr, n_cages, n_atoms);
}

seams_status seams_signature_cages(seams_system *sys, double cutoff, int type, int graph,
                                   const char *spec, size_t *n_cages, size_t *n_atoms) {
  return signature_with_edge(sys, cutoff, type, graph, spec, nullptr, n_cages, n_atoms);
}

seams_status seams_signature_cages_edge(seams_system *sys, double cutoff, int type, int graph,
                                        const char *spec, const char *edge, size_t *n_cages,
                                        size_t *n_atoms) {
  return signature_with_edge(sys, cutoff, type, graph, spec, edge, n_cages, n_atoms);
}

seams_status seams_steinhardt(seams_system *sys, double cutoff, int type, int order_l,
                              double *ql, size_t ql_len, double *qlbar, size_t qlbar_len) {
  return seams_ql(sys, cutoff, type, order_l, SEAMS_GRAPH_CUTOFF, ql, ql_len, qlbar,
                  qlbar_len);
}

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
    const auto fair = use_threads();
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
    // Parking the OpenMP team at the start of a k-nearest call made the
    // next search late by a fraction of a millisecond, which is the whole
    // of a few-thousand-atom triclinic query. The linkcell pool sleeps
    // on its own. The handle stays in the live set so destroy frees it.
    nneigh::kNearestCsr(*src, k, cutoff + 1.5, typ, mutual != 0, nb->offsets, nb->cols,
                        false);
    *out = static_cast<seams_neighbours *>(track(nb.release()));
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
    const auto fair = use_threads();
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
    const auto fair = use_threads();
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
    const auto fair = use_threads();
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
  if (!(cutoff > 0.0) || !(dist_cutoff > 0.0) || !(angle_cutoff > 0.0) ||
      !std::isfinite(cutoff) || !std::isfinite(dist_cutoff) || !std::isfinite(angle_cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoffs must be positive");
  }
  if (hydrogen_type <= 0) {
    return fail(SEAMS_ERR_ARG, "hydrogen type must be positive");
  }
  try {
    const auto fair = use_threads();
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
        if (oxygen_type <= 0) {
          heavy.pts.back().type = 1;
        }
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
    const int typ = oxygen_type > 0 ? oxygen_type : 1;
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

seams_status seams_hbond_edges(seams_system *sys, double cutoff, int oxygen_type,
                               int hydrogen_type, double dist_cutoff, double angle_cutoff,
                               int fields, seams_hbond_edge *edges, size_t cap, size_t *n_out) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "hbond edge output is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !(dist_cutoff > 0.0) || !(angle_cutoff > 0.0) ||
      !std::isfinite(cutoff) || !std::isfinite(dist_cutoff) || !std::isfinite(angle_cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoffs must be positive");
  }
  if (hydrogen_type <= 0) {
    return fail(SEAMS_ERR_ARG, "hydrogen type must be positive");
  }
  if ((fields & ~(SEAMS_HBOND_EDGE_WITH_H | SEAMS_HBOND_EDGE_WITH_GEOMETRY)) != 0) {
    return fail(SEAMS_ERR_ARG, "unknown hbond edge fields");
  }
  try {
    const auto fair = use_threads();
    const bool with_h = (fields & SEAMS_HBOND_EDGE_WITH_H) != 0;
    const bool with_geometry = (fields & SEAMS_HBOND_EDGE_WITH_GEOMETRY) != 0;
    const auto detailed = hbond_edges_detailed(sys->cloud, cutoff, oxygen_type, hydrogen_type,
                                               dist_cutoff, angle_cutoff, with_h,
                                               with_geometry);
    *n_out = detailed.size();
    if (edges == nullptr) {
      return SEAMS_OK;
    }
    if (cap > 0 && edges[0].size < sizeof(uint32_t)) {
      return fail(SEAMS_ERR_ARG, "hbond edge struct size is invalid");
    }
    if (cap < detailed.size()) {
      return fail(SEAMS_ERR_RANGE, "buffer is short");
    }
    if (detailed.empty()) {
      return SEAMS_OK;
    }
    const std::size_t elem_size = edges[0].size;
    auto *raw_out = reinterpret_cast<unsigned char *>(edges);
    for (std::size_t i = 0; i < detailed.size(); ++i) {
      uint32_t row_size = 0;
      std::memcpy(&row_size, raw_out + i * elem_size, sizeof(row_size));
      if (row_size != elem_size) {
        return fail(SEAMS_ERR_ARG, "hbond edge struct size mismatch");
      }
    }
    const std::size_t copy_bytes = std::min(elem_size, sizeof(seams_hbond_edge));
    for (std::size_t i = 0; i < detailed.size(); ++i) {
      seams_hbond_edge row{};
      row.size = static_cast<uint32_t>(elem_size);
      row.donor = static_cast<uint32_t>(detailed[i].donor);
      row.acceptor = static_cast<uint32_t>(detailed[i].acceptor);
      row.hydrogen = detailed[i].hydrogen >= 0
                         ? static_cast<uint32_t>(detailed[i].hydrogen)
                         : UINT32_MAX;
      row.d_DA = detailed[i].d_DA;
      row.angle = detailed[i].angle;
      std::memcpy(raw_out + i * elem_size, &row, copy_bytes);
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

seams_status seams_directed_hbonds(
    seams_system *sys, double pair_cutoff, double dist_cutoff, double angle_cutoff,
    const int *donor_types, size_t n_donor_types, const int *acceptor_types,
    size_t n_acceptor_types, int hydrogen_type, int exclude_same_mol,
    seams_directed_hbond *edges, size_t cap, size_t *n_out) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_out == nullptr) {
    return fail(SEAMS_ERR_NULL, "directed hbond output is null");
  }
  if ((n_donor_types > 0 && donor_types == nullptr) ||
      (n_acceptor_types > 0 && acceptor_types == nullptr)) {
    return fail(SEAMS_ERR_NULL, "type set is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(pair_cutoff > 0.0) || !(dist_cutoff > 0.0) || !(angle_cutoff > 0.0) ||
      !std::isfinite(pair_cutoff) || !std::isfinite(dist_cutoff) ||
      !std::isfinite(angle_cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoffs must be positive");
  }
  if (hydrogen_type <= 0) {
    return fail(SEAMS_ERR_ARG, "hydrogen type must be positive");
  }
  if (exclude_same_mol != 0 && exclude_same_mol != 1) {
    return fail(SEAMS_ERR_ARG, "exclude_same_mol is 0 or 1");
  }
  try {
    for (size_t i = 0; i < n_donor_types; ++i) {
      if (donor_types[i] == hydrogen_type) {
        return fail(SEAMS_ERR_ARG, "hydrogen type is distinct from the donor types");
      }
    }
    for (size_t i = 0; i < n_acceptor_types; ++i) {
      if (acceptor_types[i] == hydrogen_type) {
        return fail(SEAMS_ERR_ARG, "hydrogen type is distinct from the acceptor types");
      }
    }
    const std::vector<int> donors(donor_types, donor_types + n_donor_types);
    const std::vector<int> acceptors(acceptor_types, acceptor_types + n_acceptor_types);
    const auto detailed = bond::directedHydrogenBonds(
        sys->cloud, pair_cutoff, dist_cutoff, angle_cutoff, donors, acceptors, hydrogen_type,
        exclude_same_mol != 0);
    *n_out = detailed.size();
    if (edges == nullptr) {
      return SEAMS_OK;
    }
    if (cap > 0 && edges[0].size < sizeof(seams_directed_hbond)) {
      return fail(SEAMS_ERR_ARG, "directed hbond struct size is short");
    }
    if (cap < detailed.size()) {
      return fail(SEAMS_ERR_RANGE, "buffer is short");
    }
    if (detailed.empty()) {
      return SEAMS_OK;
    }
    const std::size_t elem_size = edges[0].size;
    auto *raw_out = reinterpret_cast<unsigned char *>(edges);
    for (std::size_t i = 0; i < detailed.size(); ++i) {
      uint32_t row_size = 0;
      std::memcpy(&row_size, raw_out + i * elem_size, sizeof(row_size));
      if (row_size != elem_size) {
        return fail(SEAMS_ERR_ARG, "directed hbond struct size mismatch");
      }
    }
    const std::size_t copy_bytes = std::min(elem_size, sizeof(seams_directed_hbond));
    for (std::size_t i = 0; i < detailed.size(); ++i) {
      seams_directed_hbond row{};
      row.size = static_cast<uint32_t>(elem_size);
      row.donor = static_cast<uint32_t>(detailed[i].donor);
      row.hydrogen = static_cast<uint32_t>(detailed[i].hydrogen);
      row.acceptor = static_cast<uint32_t>(detailed[i].acceptor);
      row.d_AH = detailed[i].d_AH;
      row.d_DA = detailed[i].d_DA;
      row.angle = detailed[i].angle;
      std::memcpy(raw_out + i * elem_size, &row, copy_bytes);
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

seams_status seams_quenched_analysis(
    seams_system *sys, seams_quench_callback callback, void *user_data, double cutoff, int type,
    int graph, int max_depth, int oxygen_type, int hydrogen_type, double dist_cutoff,
    double angle_cutoff, int *raw_labels, size_t raw_labels_len, int *quenched_labels,
    size_t quenched_labels_len, int counts_raw[9], int counts_quenched[9], int *rings_raw,
    int *rings_quenched, size_t ring_counts_len, int *hbond_raw, int *hbond_quenched,
    int *hbond_edge_delta, int *ring_count_delta, int *class_change_count) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || callback == nullptr || raw_labels == nullptr || quenched_labels == nullptr ||
      rings_raw == nullptr || rings_quenched == nullptr || hbond_raw == nullptr ||
      hbond_quenched == nullptr || hbond_edge_delta == nullptr || ring_count_delta == nullptr ||
      class_change_count == nullptr) {
    return fail(SEAMS_ERR_NULL, "quenched analysis inputs are null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!(cutoff > 0.0) || !(dist_cutoff > 0.0) || !(angle_cutoff > 0.0) ||
      !std::isfinite(cutoff) || !std::isfinite(dist_cutoff) || !std::isfinite(angle_cutoff)) {
    return fail(SEAMS_ERR_RANGE, "cutoffs must be positive");
  }
  if (max_depth < 3 || max_depth > 12) {
    return fail(SEAMS_ERR_RANGE, "max_depth must be from 3 to 12");
  }
  if (ring_counts_len < static_cast<size_t>(max_depth) + 1) {
    return fail(SEAMS_ERR_RANGE, "ring buffer is short");
  }
  if (hydrogen_type <= 0) {
    return fail(SEAMS_ERR_ARG, "hydrogen type must be positive");
  }
  if (graph != SEAMS_GRAPH_CUTOFF && graph != SEAMS_GRAPH_MUTUAL) {
    return fail(SEAMS_ERR_ARG, "graph is cutoff or mutual");
  }
  const size_t n = static_cast<size_t>(sys->cloud.nop);
  if (raw_labels_len < n || quenched_labels_len < n) {
    return fail(SEAMS_ERR_RANGE, "label buffer is short");
  }
  try {
    const auto fair = use_threads();
    auto analyse_one = [&](Cloud frame_cloud, int *labels_out, int counts9[9], int *rings_out,
                           int *hbond_out, std::vector<hbond_edge_internal> &edges_out) -> seams_status {
      int typ = 0;
      Cloud work = unified_copy(frame_cloud, type, typ);
      auto nlist = bond_rows(work, cutoff, typ, graph);
      chill::applyChillPlus(work, nlist, chill::chillPlusRule(""), false);
      write_ice_labels(work, type, labels_out, counts9);
      auto idx = nneigh::neighbourListByIndex(work, nlist);
      auto rings = primitive::ringNetwork(idx, max_depth);
      std::fill(rings_out, rings_out + ring_counts_len, 0);
      for (const auto &ring : rings) {
        const int size = static_cast<int>(ring.size());
        if (size >= 0 && static_cast<size_t>(size) < ring_counts_len) {
          rings_out[size] += 1;
        }
      }
      edges_out = hbond_edges_detailed(frame_cloud, cutoff, oxygen_type, hydrogen_type, dist_cutoff,
                                       angle_cutoff, true, true);
      *hbond_out = static_cast<int>(edges_out.size());
      return SEAMS_OK;
    };

    Cloud raw_cloud = sys->cloud;
    std::vector<hbond_edge_internal> raw_edges;
    seams_status rc = analyse_one(raw_cloud, raw_labels, counts_raw, rings_raw, hbond_raw, raw_edges);
    if (rc != SEAMS_OK) {
      return rc;
    }

    double raw_box9[9];
    {
      double H[3][3];
      double origin[3];
      nneigh::dumpBoundsToH(sys->cloud.box, sys->cloud.boxLow, H, origin);
      raw_box9[0] = origin[0];
      raw_box9[1] = origin[1];
      raw_box9[2] = origin[2];
      raw_box9[3] = H[0][0];
      raw_box9[4] = H[1][1];
      raw_box9[5] = H[2][2];
      raw_box9[6] = H[1][0];
      raw_box9[7] = H[2][0];
      raw_box9[8] = H[2][1];
    }
    std::vector<double> xyz_in(3 * n);
    std::vector<double> xyz_out(3 * n);
    for (size_t i = 0; i < n; ++i) {
      xyz_in[3 * i] = sys->cloud.pts[i].x;
      xyz_in[3 * i + 1] = sys->cloud.pts[i].y;
      xyz_in[3 * i + 2] = sys->cloud.pts[i].z;
    }
    std::memcpy(xyz_out.data(), xyz_in.data(), xyz_in.size() * sizeof(double));
    double box_out[9];
    std::memcpy(box_out, raw_box9, sizeof(box_out));
    rc = callback(xyz_in.data(), n, raw_box9, 9, xyz_out.data(), xyz_out.size(), box_out, 9, user_data);
    if (rc != SEAMS_OK) {
      return fail(rc, "quenched callback failed");
    }
    Cloud quenched_cloud = sys->cloud;
    for (size_t i = 0; i < n; ++i) {
      const double x = xyz_out[3 * i];
      const double y = xyz_out[3 * i + 1];
      const double z = xyz_out[3 * i + 2];
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
        return fail(SEAMS_ERR_ARG, "quenched callback returned non-finite coordinates");
      }
      quenched_cloud.pts[i].x = x;
      quenched_cloud.pts[i].y = y;
      quenched_cloud.pts[i].z = z;
    }
    if (!finite_box(box_out)) {
      return fail(SEAMS_ERR_BOX, "quenched box contains a non-finite number");
    }
    if (!nneigh::setBoxFromLammpsEdges(quenched_cloud.box, quenched_cloud.boxLow, box_out[0],
                                       box_out[1], box_out[2], box_out[3], box_out[4], box_out[5],
                                       box_out[6], box_out[7], box_out[8])) {
      return fail(SEAMS_ERR_BOX, "box edge lengths must be positive");
    }
    std::vector<hbond_edge_internal> quenched_edges;
    rc = analyse_one(quenched_cloud, quenched_labels, counts_quenched, rings_quenched, hbond_quenched,
                     quenched_edges);
    if (rc != SEAMS_OK) {
      return rc;
    }
    std::set<std::tuple<int, int, int>> raw_set;
    std::set<std::tuple<int, int, int>> quenched_set;
    for (const auto &edge : raw_edges) {
      raw_set.emplace(edge.donor, edge.acceptor, edge.hydrogen);
    }
    for (const auto &edge : quenched_edges) {
      quenched_set.emplace(edge.donor, edge.acceptor, edge.hydrogen);
    }
    int shared = 0;
    auto it_r = raw_set.begin();
    auto it_q = quenched_set.begin();
    while (it_r != raw_set.end() && it_q != quenched_set.end()) {
      if (*it_r == *it_q) {
        ++shared;
        ++it_r;
        ++it_q;
      } else if (*it_r < *it_q) {
        ++it_r;
      } else {
        ++it_q;
      }
    }
    *hbond_edge_delta =
        static_cast<int>(raw_set.size() + quenched_set.size() - static_cast<size_t>(2 * shared));
    int ring_delta = 0;
    for (size_t i = 0; i < ring_counts_len; ++i) {
      ring_delta += std::abs(rings_raw[i] - rings_quenched[i]);
    }
    *ring_count_delta = ring_delta;
    int changed = 0;
    for (size_t i = 0; i < n; ++i) {
      if (raw_labels[i] != quenched_labels[i]) {
        ++changed;
      }
    }
    *class_change_count = changed;
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
                               int hops, int color_types, seams_topo **out) {
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
    const auto fair = use_threads();
    int typ = 1;
    Cloud work = unified_copy(sys->cloud, type, typ);
    std::vector<std::vector<int>> rows;
    if (graph == SEAMS_GRAPH_CUTOFF) {
      rows = bond_rows(work, cutoff, typ, graph);
    } else {
      rows = nneigh::kNearestNeighbourList(work, k, cutoff + 1.5, typ, true);
    }
    auto idx = nneigh::neighbourListByIndex(work, rows);
    std::vector<int> colors;
    if (color_types != 0) {
      colors.reserve(static_cast<std::size_t>(work.nop));
      for (int i = 0; i < work.nop; ++i) {
        colors.push_back(work.pts[static_cast<std::size_t>(i)].type);
      }
    }
    const auto fp = topo::fingerprint(idx, hops, 7, colors);
    auto handle = std::make_unique<seams_topo>();
    handle->frame_key = fp.key;
    handle->method = fp.method;
    handle->atom_keys = fp.atomKeys;
    handle->ring_census = fp.ringCensus;
    handle->hops = fp.hops;
    *out = static_cast<seams_topo *>(track(handle.release()));
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
  if (!claim(fp)) {
    return;
  }
  delete fp;
}

seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap, size_t *n_out) {
  read_guard guard(fp);
  clear_error();
  if (fp == nullptr) {
    return fail(SEAMS_ERR_NULL, "fingerprint handle is null");
  }
  return copy_text(fp->frame_key, buf, cap, n_out);
}

seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf, size_t cap,
                                    size_t *n_out) {
  read_guard guard(fp);
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
  read_guard guard(fp);
  clear_error();
  if (fp == nullptr || counts == nullptr) {
    return fail(SEAMS_ERR_NULL, "fingerprint ring output is null");
  }
  if (counts_len == 0) {
    return fail(SEAMS_ERR_RANGE, "ring census length is empty");
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
    const auto fair = use_threads();
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
    const auto fair = use_threads();
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
  if (!(cutoff > 0.0) || !(sheet_area > 0.0) || !std::isfinite(cutoff) ||
      !std::isfinite(sheet_area)) {
    return fail(SEAMS_ERR_RANGE, "cutoff and sheet area must be positive");
  }
  try {
    const auto fair = use_threads();
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
  read_guard guard(front);
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

const char *seams_last_error_message(void) { return g_error.c_str(); }

seams_status seams_system_natoms(const seams_system *sys, size_t *n) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || n == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or count is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  *n = static_cast<size_t>(sys->cloud.nop);
  return SEAMS_OK;
}

seams_status seams_neighbours_nrows(const seams_neighbours *nb, size_t *n) {
  read_guard guard(nb);
  clear_error();
  if (nb == nullptr || n == nullptr) {
    return fail(SEAMS_ERR_NULL, "neighbours or count is null");
  }
  *n = static_cast<size_t>(nb->n);
  return SEAMS_OK;
}

seams_status seams_rings_nrings(const seams_rings *rings, size_t *n) {
  read_guard guard(rings);
  clear_error();
  if (rings == nullptr || n == nullptr) {
    return fail(SEAMS_ERR_NULL, "rings or count is null");
  }
  *n = rings->rings.size();
  return SEAMS_OK;
}

seams_status seams_chill_plus_n(seams_system *sys, double cutoff, int type, int graph,
                                int rule, int *labels, size_t labels_len, int counts9[9]) {
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
  return seams_chill_plus(sys, cutoff, type, graph, rule, labels, counts9);
}

seams_status seams_q4q6_n(seams_system *sys, double cutoff, int type, double *q4,
                          size_t q4_len, double *q6, size_t q6_len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  if (q4 == nullptr && q6 == nullptr) {
    return fail(SEAMS_ERR_NULL, "q4 and q6 are both null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if ((q4 != nullptr && !buffer_fits(q4_len, sys->cloud.nop)) ||
      (q6 != nullptr && !buffer_fits(q6_len, sys->cloud.nop))) {
    return fail(SEAMS_ERR_RANGE, "q buffer is short");
  }
  return seams_q4q6(sys, cutoff, type, q4, q6);
}

seams_status seams_rings_labels_n(const seams_rings *rings, int *sizes, size_t sizes_len,
                                  int *symmetry, size_t symmetry_len) {
  read_guard guard(rings);
  clear_error();
  if (rings == nullptr) {
    return fail(SEAMS_ERR_NULL, "rings handle is null");
  }
  if (sizes == nullptr && symmetry == nullptr) {
    return fail(SEAMS_ERR_NULL, "sizes and symmetry are both null");
  }
  const size_t need = rings->rings.size();
  if ((sizes != nullptr && sizes_len < need) ||
      (symmetry != nullptr && symmetry_len < need)) {
    return fail(SEAMS_ERR_RANGE, "ring label buffer is short");
  }
  return seams_rings_labels(rings, sizes, symmetry);
}

seams_status seams_rings_histogram_n(const seams_rings *rings, int *counts,
                                     size_t counts_len) {
  read_guard guard(rings);
  clear_error();
  if (rings == nullptr || counts == nullptr) {
    return fail(SEAMS_ERR_NULL, "histogram output is null");
  }
  if (counts_len == 0) {
    return fail(SEAMS_ERR_RANGE, "histogram length is empty");
  }
  for (size_t i = 0; i < counts_len; ++i) {
    counts[i] = 0;
  }
  for (const auto &ring : rings->rings) {
    if (ring.size() < counts_len) {
      counts[ring.size()] += 1;
    }
  }
  return SEAMS_OK;
}

seams_status seams_system_set_positions(seams_system *sys, const double *xyz, size_t n,
                                        const double *box9, size_t box_len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || xyz == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or coordinates is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (n != static_cast<std::size_t>(sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "atom count differs from the loaded frame");
  }
  if (box9 != nullptr) {
    if (box_len != 9) {
      return fail(SEAMS_ERR_RANGE, "box length is 9");
    }
    if (!finite_box(box9)) {
      return fail(SEAMS_ERR_BOX, "box contains a non-finite number");
    }
  }
  try {
    if (box9 != nullptr &&
        !nneigh::setBoxFromLammpsEdges(sys->spare_box, sys->spare_low, box9[0], box9[1],
                                       box9[2], box9[3], box9[4], box9[5], box9[6],
                                       box9[7], box9[8])) {
      return fail(SEAMS_ERR_BOX, "box edge lengths must be positive");
    }
    if (sys->stage) {
      auto &next = sys->staged;
      if (!sys->staged_current) {
        next = sys->cloud.pts;
        sys->staged_current = true;
      }
      constexpr std::uint64_t exponent = 0x7ff0000000000000ULL;
      std::uint64_t bad = 0;
      for (std::size_t i = 0; i < n; ++i) {
        std::uint64_t bits[3];
        std::memcpy(bits, xyz + 3 * i, sizeof bits);
        bad |= static_cast<std::uint64_t>((bits[0] & exponent) == exponent) |
               static_cast<std::uint64_t>((bits[1] & exponent) == exponent) |
               static_cast<std::uint64_t>((bits[2] & exponent) == exponent);
        auto &point = next[i];
        point.x = xyz[3 * i];
        point.y = xyz[3 * i + 1];
        point.z = xyz[3 * i + 2];
        point.iceType = molSys::atom_state_type::unclassified;
        point.inSlice = true;
        point.c_ij.clear();
      }
      if (bad != 0) {
        return fail(SEAMS_ERR_ARG, "a coordinate is not finite");
      }
      sys->cloud.pts.swap(next);
    } else {
      if (!finite_coordinates(xyz, n)) {
        return fail(SEAMS_ERR_ARG, "a coordinate is not finite");
      }
      for (std::size_t i = 0; i < n; ++i) {
        auto &point = sys->cloud.pts[i];
        point.x = xyz[3 * i];
        point.y = xyz[3 * i + 1];
        point.z = xyz[3 * i + 2];
        point.iceType = molSys::atom_state_type::unclassified;
        point.inSlice = true;
        point.c_ij.clear();
      }
    }
    if (box9 != nullptr) {
      sys->cloud.box.swap(sys->spare_box);
      sys->cloud.boxLow.swap(sys->spare_low);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_system_stage_positions(seams_system *sys, int enable) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr) {
    return fail(SEAMS_ERR_NULL, "system is null");
  }
  sys->stage = enable != 0;
  if (!sys->stage) {
    decltype(Cloud::pts)().swap(sys->staged);
    sys->staged_current = false;
  }
  return SEAMS_OK;
}

seams_status seams_system_box(const seams_system *sys, double *box9, size_t box_len) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || box9 == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or box is null");
  }
  if (box_len != 9) {
    return fail(SEAMS_ERR_RANGE, "box length is 9");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  double H[3][3];
  double origin[3];
  nneigh::dumpBoundsToH(sys->cloud.box, sys->cloud.boxLow, H, origin);
  box9[0] = origin[0];
  box9[1] = origin[1];
  box9[2] = origin[2];
  box9[3] = H[0][0];
  box9[4] = H[1][1];
  box9[5] = H[2][2];
  box9[6] = H[1][0];
  box9[7] = H[2][0];
  box9[8] = H[2][1];
  return SEAMS_OK;
}

seams_status seams_system_positions(const seams_system *sys, double *xyz, size_t len) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || xyz == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or coordinates is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  const std::size_t n = static_cast<std::size_t>(sys->cloud.nop);
  if (len / 3 < n) {
    return fail(SEAMS_ERR_RANGE, "coordinate buffer is short");
  }
  for (std::size_t i = 0; i < n; ++i) {
    const auto &point = sys->cloud.pts[i];
    xyz[3 * i] = point.x;
    xyz[3 * i + 1] = point.y;
    xyz[3 * i + 2] = point.z;
  }
  return SEAMS_OK;
}

seams_status seams_system_types(const seams_system *sys, int *types, size_t len) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || types == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or types is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!buffer_fits(len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "type buffer is short");
  }
  for (int i = 0; i < sys->cloud.nop; ++i) {
    types[i] = sys->cloud.pts[static_cast<std::size_t>(i)].type;
  }
  return SEAMS_OK;
}

seams_status seams_system_ids(const seams_system *sys, int *ids, size_t len) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || ids == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or ids is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!buffer_fits(len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "id buffer is short");
  }
  for (int i = 0; i < sys->cloud.nop; ++i) {
    ids[i] = sys->cloud.pts[static_cast<std::size_t>(i)].atomID;
  }
  return SEAMS_OK;
}

seams_status seams_mols(const seams_system *sys, int *mols, size_t len) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || mols == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or mols is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (!buffer_fits(len, sys->cloud.nop)) {
    return fail(SEAMS_ERR_RANGE, "molecule id buffer is short");
  }
  for (int i = 0; i < sys->cloud.nop; ++i) {
    mols[i] = sys->cloud.pts[static_cast<std::size_t>(i)].molID;
  }
  return SEAMS_OK;
}

seams_status seams_set_mols(seams_system *sys, const int *mols, size_t len) {
  handle_guard guard(sys);
  clear_error();
  if (sys == nullptr || mols == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or mols is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  const auto n = static_cast<std::size_t>(sys->cloud.nop);
  if (len != n) {
    return fail(SEAMS_ERR_RANGE, "molecule id length does not match atom count");
  }
  for (std::size_t i = 0; i < n; ++i) {
    sys->cloud.pts[i].molID = mols[i];
  }
  return SEAMS_OK;
}

seams_status seams_molecule_centres(const seams_system *sys, const int *mols,
                                    size_t mols_len, const double *masses,
                                    size_t masses_len, int *mol_ids_out,
                                    double *centres_out, size_t cap,
                                    size_t *n_mols) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_mols == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or molecule count is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  const auto n = static_cast<std::size_t>(sys->cloud.nop);
  if ((mols == nullptr && mols_len != 0) ||
      (mols != nullptr && mols_len != n)) {
    return fail(SEAMS_ERR_RANGE, "mols length must match atom count");
  }
  if ((masses == nullptr && masses_len != 0) ||
      (masses != nullptr && masses_len != n)) {
    return fail(SEAMS_ERR_RANGE, "masses length must match atom count");
  }
  if (masses != nullptr) {
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(masses[i]) || !(masses[i] > 0.0)) {
        return fail(SEAMS_ERR_ARG, "masses must be finite and positive");
      }
    }
  }
  try {
    const traj::BoxCell cell = traj::BoxCell::fromCloud(sys->cloud);
    if (!cell.ok) {
      return fail(SEAMS_ERR_BOX, "stored cell is invalid");
    }
    std::vector<double> xyz(3 * n);
    std::vector<int> stored_mols;
    if (mols == nullptr) {
      stored_mols.resize(n);
    }
    for (std::size_t i = 0; i < n; ++i) {
      xyz[3 * i] = sys->cloud.pts[i].x;
      xyz[3 * i + 1] = sys->cloud.pts[i].y;
      xyz[3 * i + 2] = sys->cloud.pts[i].z;
      if (mols == nullptr) {
        stored_mols[i] = sys->cloud.pts[i].molID;
      }
    }
    std::vector<int> unique_mols;
    std::vector<double> centres;
    const int *use_mols = mols != nullptr ? mols : stored_mols.data();
    if (!traj::moleculeCentres(xyz.data(), n, cell, use_mols, masses,
                               unique_mols, centres)) {
      return fail(SEAMS_ERR_ARG, "molecule centre assembly failed");
    }
    const std::size_t nm = unique_mols.size();
    *n_mols = nm;
    if (mol_ids_out == nullptr && centres_out == nullptr) {
      return SEAMS_OK;
    }
    if (cap < nm) {
      return fail(SEAMS_ERR_RANGE, "molecule output buffer is short");
    }
    if (mol_ids_out != nullptr && nm > 0) {
      std::copy_n(unique_mols.data(), nm, mol_ids_out);
    }
    if (centres_out != nullptr && nm > 0) {
      std::copy_n(centres.data(), 3 * nm, centres_out);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_traj_unwrap_positions(const double *xyz, size_t nframes,
                                         size_t natoms, const double *boxes,
                                         size_t nboxes, double *unwrapped,
                                         size_t unwrapped_len) {
  clear_error();
  if (xyz == nullptr || boxes == nullptr || unwrapped == nullptr) {
    return fail(SEAMS_ERR_NULL, "coordinates, boxes, or output is null");
  }
  if (nframes == 0 || natoms == 0) {
    return fail(SEAMS_ERR_RANGE, "nframes and natoms must be positive");
  }
  if (nboxes != 1 && nboxes != nframes) {
    return fail(SEAMS_ERR_RANGE, "nboxes must be 1 or nframes");
  }
  if (natoms > SIZE_MAX / 3 || nframes > SIZE_MAX / (3 * natoms)) {
    return fail(SEAMS_ERR_RANGE, "coordinate count overflows size_t");
  }
  const std::size_t total = nframes * 3 * natoms;
  if (unwrapped_len < total) {
    return fail(SEAMS_ERR_RANGE, "unwrapped buffer is short");
  }
  for (std::size_t b = 0; b < nboxes; ++b) {
    const traj::BoxCell cell = traj::BoxCell::fromBox9(boxes + 9 * b);
    if (!cell.ok) {
      return fail(SEAMS_ERR_BOX, "periodic box is invalid or non-finite");
    }
  }
  if (!finite_coordinates(xyz, nframes * natoms)) {
    return fail(SEAMS_ERR_ARG, "a coordinate is not finite");
  }
  try {
    if (!traj::unwrapPositions(xyz, nframes, natoms, boxes, nboxes,
                               unwrapped)) {
      return fail(SEAMS_ERR_ARG, "failed to unwrap coordinates");
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_traj_unwrap_centres(
    const double *xyz, size_t nframes, size_t natoms, const double *boxes,
    size_t nboxes, const int *mols, size_t mols_len, const double *masses,
    size_t masses_len, int *mol_ids_out, size_t mol_ids_cap,
    double *centres_out, size_t centres_len, size_t *n_mols) {
  clear_error();
  if (xyz == nullptr || boxes == nullptr || mols == nullptr ||
      n_mols == nullptr) {
    return fail(SEAMS_ERR_NULL, "coordinates, boxes, mols, or n_mols is null");
  }
  if (nframes == 0 || natoms == 0) {
    return fail(SEAMS_ERR_RANGE, "nframes and natoms must be positive");
  }
  if (nboxes != 1 && nboxes != nframes) {
    return fail(SEAMS_ERR_RANGE, "nboxes must be 1 or nframes");
  }
  if (mols_len != natoms) {
    return fail(SEAMS_ERR_RANGE, "mols_len must equal natoms");
  }
  if ((masses == nullptr && masses_len != 0) ||
      (masses != nullptr && masses_len != natoms)) {
    return fail(SEAMS_ERR_RANGE, "masses_len must equal natoms or 0");
  }
  if (natoms > SIZE_MAX / 3 || nframes > SIZE_MAX / (3 * natoms)) {
    return fail(SEAMS_ERR_RANGE, "coordinate count overflows size_t");
  }
  for (std::size_t b = 0; b < nboxes; ++b) {
    const traj::BoxCell cell = traj::BoxCell::fromBox9(boxes + 9 * b);
    if (!cell.ok) {
      return fail(SEAMS_ERR_BOX, "periodic box is invalid or non-finite");
    }
  }
  if (!finite_coordinates(xyz, nframes * natoms)) {
    return fail(SEAMS_ERR_ARG, "a coordinate is not finite");
  }
  if (masses != nullptr) {
    for (std::size_t i = 0; i < natoms; ++i) {
      if (!std::isfinite(masses[i]) || !(masses[i] > 0.0)) {
        return fail(SEAMS_ERR_ARG, "masses must be finite and positive");
      }
    }
  }
  try {
    std::vector<int> unique_mols;
    std::vector<double> centres;
    if (!traj::unwrapMoleculeCentres(xyz, nframes, natoms, boxes, nboxes, mols,
                                     masses, unique_mols, centres)) {
      return fail(SEAMS_ERR_ARG, "failed to unwrap molecule centres");
    }
    const std::size_t nm = unique_mols.size();
    *n_mols = nm;
    if (mol_ids_out == nullptr && centres_out == nullptr) {
      return SEAMS_OK;
    }
    if (mol_ids_out != nullptr && mol_ids_cap < nm) {
      return fail(SEAMS_ERR_RANGE, "mol_ids_out buffer is short");
    }
    if (centres_out != nullptr && centres_len < nframes * 3 * nm) {
      return fail(SEAMS_ERR_RANGE, "centres_out buffer is short");
    }
    if (mol_ids_out != nullptr && nm > 0) {
      std::copy_n(unique_mols.data(), nm, mol_ids_out);
    }
    if (centres_out != nullptr && !centres.empty()) {
      std::copy(centres.begin(), centres.end(), centres_out);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_water_vectors(const seams_system *sys, int oxygen_type,
                                 int hydrogen_type, int *mol_ids,
                                 double *oh_mid, double *bisector, double *hh,
                                 double *hh_unit, size_t cap,
                                 size_t *n_waters) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_waters == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or n_waters is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  if (oxygen_type == hydrogen_type) {
    return fail(SEAMS_ERR_ARG, "oxygen_type and hydrogen_type must differ");
  }
  try {
    traj::WaterVectorsResult res;
    if (!traj::waterOrientationVectors(sys->cloud, oxygen_type, hydrogen_type,
                                       res)) {
      return fail(SEAMS_ERR_BOX, "stored cell is invalid");
    }
    const std::size_t nw = res.molIds.size();
    *n_waters = nw;
    if (mol_ids == nullptr && oh_mid == nullptr && bisector == nullptr &&
        hh == nullptr && hh_unit == nullptr) {
      return SEAMS_OK;
    }
    if (cap < nw) {
      return fail(SEAMS_ERR_RANGE, "water vector buffer is short");
    }
    if (mol_ids != nullptr && nw > 0) {
      std::copy_n(res.molIds.data(), nw, mol_ids);
    }
    if (oh_mid != nullptr && nw > 0) {
      std::copy_n(res.ohMidVectors.data(), 3 * nw, oh_mid);
    }
    if (bisector != nullptr && nw > 0) {
      std::copy_n(res.bisectorUnit.data(), 3 * nw, bisector);
    }
    if (hh != nullptr && nw > 0) {
      std::copy_n(res.hhVectors.data(), 3 * nw, hh);
    }
    if (hh_unit != nullptr && nw > 0) {
      std::copy_n(res.hhUnit.data(), 3 * nw, hh_unit);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_molecule_dipoles(
    const seams_system *sys, const int *mols, size_t mols_len,
    const double *charges, size_t charges_len, int *mol_ids_out,
    double *dipoles, double *unit_vectors, size_t cap, size_t *n_mols) {
  read_guard guard(sys);
  clear_error();
  if (sys == nullptr || n_mols == nullptr) {
    return fail(SEAMS_ERR_NULL, "system or n_mols is null");
  }
  if (!sys->ready) {
    return fail(SEAMS_ERR_STATE, "system has no frame");
  }
  const auto n = static_cast<std::size_t>(sys->cloud.nop);
  if ((mols == nullptr && mols_len != 0) ||
      (mols != nullptr && mols_len != n)) {
    return fail(SEAMS_ERR_RANGE, "mols_len must equal atom count or 0");
  }
  if ((charges == nullptr && charges_len != 0) ||
      (charges != nullptr && charges_len != n)) {
    return fail(SEAMS_ERR_RANGE, "charges_len must equal atom count or 0");
  }
  if (charges != nullptr && !finite_scan<4>(charges, n)) {
    return fail(SEAMS_ERR_ARG, "charges contain a non-finite value");
  }
  try {
    traj::DipoleVectorsResult res;
    if (!traj::moleculeDipoleVectors(sys->cloud, mols, charges, res)) {
      return fail(SEAMS_ERR_BOX, "stored cell is invalid");
    }
    const std::size_t nm = res.molIds.size();
    *n_mols = nm;
    if (mol_ids_out == nullptr && dipoles == nullptr &&
        unit_vectors == nullptr) {
      return SEAMS_OK;
    }
    if (cap < nm) {
      return fail(SEAMS_ERR_RANGE, "dipole output buffer is short");
    }
    if (mol_ids_out != nullptr && nm > 0) {
      std::copy_n(res.molIds.data(), nm, mol_ids_out);
    }
    if (dipoles != nullptr && nm > 0) {
      std::copy_n(res.dipoles.data(), 3 * nm, dipoles);
    }
    if (unit_vectors != nullptr && nm > 0) {
      std::copy_n(res.unitVectors.data(), 3 * nm, unit_vectors);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_orientation_order(const double *vectors, size_t n_vectors,
                                     const double axis[3], double *cos_theta,
                                     double *p2, size_t cap, double *mean_p1,
                                     double *mean_p2) {
  clear_error();
  if (axis == nullptr) {
    return fail(SEAMS_ERR_NULL, "axis is null");
  }
  if (n_vectors > 0 && vectors == nullptr) {
    return fail(SEAMS_ERR_NULL, "vectors is null");
  }
  if (cos_theta == nullptr && p2 == nullptr && mean_p1 == nullptr &&
      mean_p2 == nullptr) {
    return fail(SEAMS_ERR_NULL, "at least one output pointer is required");
  }
  if ((cos_theta != nullptr || p2 != nullptr) && cap < n_vectors) {
    return fail(SEAMS_ERR_RANGE, "per-vector output buffer is short");
  }
  if (!std::isfinite(axis[0]) || !std::isfinite(axis[1]) ||
      !std::isfinite(axis[2])) {
    return fail(SEAMS_ERR_ARG, "axis contains a non-finite number");
  }
  const double anorm =
      std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  if (!(anorm > 0.0) || !std::isfinite(anorm)) {
    return fail(SEAMS_ERR_ARG, "axis must be non-zero");
  }
  if (n_vectors > 0 && !finite_coordinates(vectors, n_vectors)) {
    return fail(SEAMS_ERR_ARG, "vectors contain a non-finite number");
  }
  try {
    traj::OrientationOrderResult res;
    if (!traj::axisOrderParameters(vectors, n_vectors, axis, res)) {
      return fail(SEAMS_ERR_ARG, "orientation order computation failed");
    }
    if (cos_theta != nullptr && n_vectors > 0) {
      std::copy_n(res.cosTheta.data(), n_vectors, cos_theta);
    }
    if (p2 != nullptr && n_vectors > 0) {
      std::copy_n(res.p2.data(), n_vectors, p2);
    }
    if (mean_p1 != nullptr) {
      *mean_p1 = res.meanP1;
    }
    if (mean_p2 != nullptr) {
      *mean_p2 = res.meanP2;
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_msd_create(size_t natoms, size_t max_lag,
                              size_t origin_stride, double dt,
                              const double axis[3], seams_msd **out) {
  clear_error();
  if (out != nullptr) {
    *out = nullptr;
  }
  if (axis == nullptr || out == nullptr) {
    return fail(SEAMS_ERR_NULL, "axis or output pointer is null");
  }
  if (natoms == 0 || origin_stride == 0) {
    return fail(SEAMS_ERR_RANGE, "natoms and origin_stride must be positive");
  }
  if (!std::isfinite(dt) || !(dt > 0.0)) {
    return fail(SEAMS_ERR_ARG, "dt must be finite and positive");
  }
  if (!std::isfinite(axis[0]) || !std::isfinite(axis[1]) ||
      !std::isfinite(axis[2])) {
    return fail(SEAMS_ERR_ARG, "axis contains a non-finite number");
  }
  const double anorm =
      std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  if (!(anorm > 0.0) || !std::isfinite(anorm)) {
    return fail(SEAMS_ERR_ARG, "axis must be non-zero");
  }
  try {
    auto holder =
        std::make_unique<seams_msd>(natoms, max_lag, origin_stride, dt, axis);
    if (!holder->acc.valid()) {
      return fail(SEAMS_ERR_ARG, "invalid MSD accumulator parameters");
    }
    *out = static_cast<seams_msd *>(track(holder.release()));
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_msd_destroy(seams_msd *msd) {
  clear_error();
  if (claim(msd)) {
    delete msd;
  }
}

seams_status seams_msd_push(seams_msd *msd, const double *xyz, size_t natoms,
                            const double *box9, size_t box_len) {
  handle_guard guard(msd);
  clear_error();
  if (msd == nullptr || xyz == nullptr) {
    return fail(SEAMS_ERR_NULL, "msd handle or coordinates is null");
  }
  if (natoms != msd->acc.numParticles()) {
    return fail(SEAMS_ERR_RANGE, "natoms does not match accumulator");
  }
  if ((box9 == nullptr && box_len != 0) || (box9 != nullptr && box_len != 9)) {
    return fail(SEAMS_ERR_RANGE, "box_len must be 9 with box9 or 0 without");
  }
  if (box9 != nullptr) {
    const traj::BoxCell cell = traj::BoxCell::fromBox9(box9);
    if (!cell.ok) {
      return fail(SEAMS_ERR_BOX, "periodic box is invalid or non-finite");
    }
  }
  if (!finite_coordinates(xyz, natoms)) {
    return fail(SEAMS_ERR_ARG, "a coordinate is not finite");
  }
  try {
    if (!msd->acc.push(xyz, natoms, box9)) {
      return fail(SEAMS_ERR_STATE,
                  "cannot mix wrapped and unwrapped frames on one accumulator");
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_msd_curve(const seams_msd *msd, double *lag_time,
                             double *total, double *axial, double *lateral,
                             size_t *n_origins, size_t cap, size_t *n_lags) {
  read_guard guard(msd);
  clear_error();
  if (msd == nullptr || n_lags == nullptr) {
    return fail(SEAMS_ERR_NULL, "msd handle or n_lags is null");
  }
  if (msd->acc.framesSeen() == 0) {
    return fail(SEAMS_ERR_STATE, "msd accumulator has no frames");
  }
  try {
    std::vector<double> lt;
    std::vector<double> tot;
    std::vector<double> ax;
    std::vector<double> lat;
    std::vector<std::size_t> orig;
    msd->acc.curve(lt, tot, ax, lat, orig);
    const std::size_t nl = lt.size();
    *n_lags = nl;
    if (lag_time == nullptr && total == nullptr && axial == nullptr &&
        lateral == nullptr && n_origins == nullptr) {
      return SEAMS_OK;
    }
    if (cap < nl) {
      return fail(SEAMS_ERR_RANGE, "msd curve output buffer is short");
    }
    if (lag_time != nullptr && nl > 0) {
      std::copy_n(lt.data(), nl, lag_time);
    }
    if (total != nullptr && nl > 0) {
      std::copy_n(tot.data(), nl, total);
    }
    if (axial != nullptr && nl > 0) {
      std::copy_n(ax.data(), nl, axial);
    }
    if (lateral != nullptr && nl > 0) {
      std::copy_n(lat.data(), nl, lateral);
    }
    if (n_origins != nullptr && nl > 0) {
      std::copy_n(orig.data(), nl, n_origins);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_msd_fit(const seams_msd *msd, double t_start, double t_end,
                           double *d_total, double *d_axial, double *d_lateral,
                           double *alpha) {
  read_guard guard(msd);
  clear_error();
  if (msd == nullptr) {
    return fail(SEAMS_ERR_NULL, "msd handle is null");
  }
  if (d_total == nullptr && d_axial == nullptr && d_lateral == nullptr &&
      alpha == nullptr) {
    return fail(SEAMS_ERR_NULL, "at least one fit output pointer is required");
  }
  if (!std::isfinite(t_start) || !std::isfinite(t_end)) {
    return fail(SEAMS_ERR_ARG, "fit window bounds must be finite");
  }
  if (t_start < 0.0 ||
      (!(t_start == 0.0 && t_end == 0.0) && t_end <= t_start)) {
    return fail(SEAMS_ERR_RANGE, "invalid fit time window");
  }
  if (msd->acc.numLags() < 2) {
    return fail(SEAMS_ERR_STATE, "msd fit needs at least two lag points");
  }
  try {
    traj::DiffusionFitResult res;
    if (!msd->acc.fit(t_start, t_end, res)) {
      return fail(SEAMS_ERR_RANGE, "fewer than two lag points in fit window");
    }
    if (d_total != nullptr) {
      *d_total = res.dTotal;
    }
    if (d_axial != nullptr) {
      *d_axial = res.dAxial;
    }
    if (d_lateral != nullptr) {
      *d_lateral = res.dLateral;
    }
    if (alpha != nullptr) {
      *alpha = res.alphaFitted ? res.alpha
                               : std::numeric_limits<double>::quiet_NaN();
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_hbond_tcf_create(size_t max_lag, size_t origin_stride,
                                    double dt, int directed,
                                    seams_hbond_tcf **out) {
  clear_error();
  if (out != nullptr) {
    *out = nullptr;
  }
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "output pointer is null");
  }
  if (origin_stride == 0) {
    return fail(SEAMS_ERR_RANGE, "origin_stride must be positive");
  }
  if (!std::isfinite(dt) || !(dt > 0.0)) {
    return fail(SEAMS_ERR_ARG, "dt must be finite and positive");
  }
  try {
    auto holder = std::make_unique<seams_hbond_tcf>(max_lag, origin_stride, dt,
                                                    directed != 0);
    if (!holder->acc.valid()) {
      return fail(SEAMS_ERR_ARG, "invalid H-bond TCF parameters");
    }
    *out = static_cast<seams_hbond_tcf *>(track(holder.release()));
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_hbond_tcf_destroy(seams_hbond_tcf *tcf) {
  clear_error();
  if (claim(tcf)) {
    delete tcf;
  }
}

seams_status seams_hbond_tcf_push(seams_hbond_tcf *tcf, const int *donors,
                                  const int *acceptors, size_t n_edges) {
  handle_guard guard(tcf);
  clear_error();
  if (tcf == nullptr) {
    return fail(SEAMS_ERR_NULL, "tcf handle is null");
  }
  if (n_edges > 0 && (donors == nullptr || acceptors == nullptr)) {
    return fail(SEAMS_ERR_NULL, "donors or acceptors is null for non-empty edges");
  }
  try {
    if (!tcf->acc.push(donors, acceptors, n_edges)) {
      return fail(SEAMS_ERR_ARG, "failed to push H-bond edges");
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_hbond_tcf_curve(const seams_hbond_tcf *tcf, double *lag_time,
                                   double *c_intermittent,
                                   double *s_continuous, size_t *n_origins,
                                   size_t cap, size_t *n_lags) {
  read_guard guard(tcf);
  clear_error();
  if (tcf == nullptr || n_lags == nullptr) {
    return fail(SEAMS_ERR_NULL, "tcf handle or n_lags is null");
  }
  if (tcf->acc.framesSeen() == 0) {
    return fail(SEAMS_ERR_STATE, "tcf accumulator has no frames");
  }
  try {
    std::vector<double> lt;
    std::vector<double> chb;
    std::vector<double> shb;
    std::vector<std::size_t> orig;
    tcf->acc.curve(lt, chb, shb, orig);
    const std::size_t nl = lt.size();
    *n_lags = nl;
    if (lag_time == nullptr && c_intermittent == nullptr &&
        s_continuous == nullptr && n_origins == nullptr) {
      return SEAMS_OK;
    }
    if (cap < nl) {
      return fail(SEAMS_ERR_RANGE, "tcf curve buffer is short");
    }
    if (lag_time != nullptr && nl > 0) {
      std::copy_n(lt.data(), nl, lag_time);
    }
    if (c_intermittent != nullptr && nl > 0) {
      std::copy_n(chb.data(), nl, c_intermittent);
    }
    if (s_continuous != nullptr && nl > 0) {
      std::copy_n(shb.data(), nl, s_continuous);
    }
    if (n_origins != nullptr && nl > 0) {
      std::copy_n(orig.data(), nl, n_origins);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_hbond_tcf_fit(const seams_hbond_tcf *tcf, double t_start,
                                 double t_end, double *tau_c_int,
                                 double *tau_c_exp, double *tau_s_int,
                                 double *tau_s_exp) {
  read_guard guard(tcf);
  clear_error();
  if (tcf == nullptr) {
    return fail(SEAMS_ERR_NULL, "tcf handle is null");
  }
  if (tau_c_int == nullptr && tau_c_exp == nullptr && tau_s_int == nullptr &&
      tau_s_exp == nullptr) {
    return fail(SEAMS_ERR_NULL, "at least one lifetime output pointer is required");
  }
  if (!std::isfinite(t_start) || !std::isfinite(t_end)) {
    return fail(SEAMS_ERR_ARG, "fit window bounds must be finite");
  }
  if (t_start < 0.0 ||
      (!(t_start == 0.0 && t_end == 0.0) && t_end <= t_start)) {
    return fail(SEAMS_ERR_RANGE, "invalid fit time window");
  }
  if (tcf->acc.numLags() < 2) {
    return fail(SEAMS_ERR_STATE, "tcf fit needs at least two lag points");
  }
  try {
    traj::DecayFitResult inter;
    traj::DecayFitResult cont;
    if (!tcf->acc.fit(t_start, t_end, inter, cont)) {
      return fail(SEAMS_ERR_RANGE, "insufficient lag points in fit window");
    }
    if (tau_c_int != nullptr) {
      *tau_c_int = inter.tauIntegral;
    }
    if (tau_c_exp != nullptr) {
      *tau_c_exp = inter.tauExp;
    }
    if (tau_s_int != nullptr) {
      *tau_s_int = cont.tauIntegral;
    }
    if (tau_s_exp != nullptr) {
      *tau_s_exp = cont.tauExp;
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_reorient_create(size_t n_vectors, size_t max_lag,
                                   size_t origin_stride, double dt,
                                   seams_reorient **out) {
  clear_error();
  if (out != nullptr) {
    *out = nullptr;
  }
  if (out == nullptr) {
    return fail(SEAMS_ERR_NULL, "output pointer is null");
  }
  if (n_vectors == 0 || origin_stride == 0) {
    return fail(SEAMS_ERR_RANGE,
                "n_vectors and origin_stride must be positive");
  }
  if (!std::isfinite(dt) || !(dt > 0.0)) {
    return fail(SEAMS_ERR_ARG, "dt must be finite and positive");
  }
  try {
    auto holder =
        std::make_unique<seams_reorient>(n_vectors, max_lag, origin_stride, dt);
    if (!holder->acc.valid()) {
      return fail(SEAMS_ERR_ARG, "invalid reorientation parameters");
    }
    *out = static_cast<seams_reorient *>(track(holder.release()));
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

void seams_reorient_destroy(seams_reorient *reorient) {
  clear_error();
  if (claim(reorient)) {
    delete reorient;
  }
}

seams_status seams_reorient_push(seams_reorient *reorient,
                                 const double *vectors, size_t n_vectors) {
  handle_guard guard(reorient);
  clear_error();
  if (reorient == nullptr || vectors == nullptr) {
    return fail(SEAMS_ERR_NULL, "reorient handle or vectors is null");
  }
  if (n_vectors != reorient->acc.numVectors()) {
    return fail(SEAMS_ERR_RANGE, "n_vectors does not match accumulator");
  }
  if (!finite_coordinates(vectors, n_vectors)) {
    return fail(SEAMS_ERR_ARG, "vectors contain a non-finite number");
  }
  try {
    if (!reorient->acc.push(vectors, n_vectors)) {
      return fail(SEAMS_ERR_ARG, "failed to push vectors");
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_reorient_curve(const seams_reorient *reorient,
                                  double *lag_time, double *c1, double *c2,
                                  size_t *n_origins, size_t cap,
                                  size_t *n_lags) {
  read_guard guard(reorient);
  clear_error();
  if (reorient == nullptr || n_lags == nullptr) {
    return fail(SEAMS_ERR_NULL, "reorient handle or n_lags is null");
  }
  if (reorient->acc.framesSeen() == 0) {
    return fail(SEAMS_ERR_STATE, "reorient accumulator has no frames");
  }
  try {
    std::vector<double> lt;
    std::vector<double> c1v;
    std::vector<double> c2v;
    std::vector<std::size_t> orig;
    reorient->acc.curve(lt, c1v, c2v, orig);
    const std::size_t nl = lt.size();
    *n_lags = nl;
    if (lag_time == nullptr && c1 == nullptr && c2 == nullptr &&
        n_origins == nullptr) {
      return SEAMS_OK;
    }
    if (cap < nl) {
      return fail(SEAMS_ERR_RANGE, "reorient curve buffer is short");
    }
    if (lag_time != nullptr && nl > 0) {
      std::copy_n(lt.data(), nl, lag_time);
    }
    if (c1 != nullptr && nl > 0) {
      std::copy_n(c1v.data(), nl, c1);
    }
    if (c2 != nullptr && nl > 0) {
      std::copy_n(c2v.data(), nl, c2);
    }
    if (n_origins != nullptr && nl > 0) {
      std::copy_n(orig.data(), nl, n_origins);
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

seams_status seams_reorient_fit(const seams_reorient *reorient, double t_start,
                                double t_end, double *tau1_int,
                                double *tau1_exp, double *tau2_int,
                                double *tau2_exp) {
  read_guard guard(reorient);
  clear_error();
  if (reorient == nullptr) {
    return fail(SEAMS_ERR_NULL, "reorient handle is null");
  }
  if (tau1_int == nullptr && tau1_exp == nullptr && tau2_int == nullptr &&
      tau2_exp == nullptr) {
    return fail(SEAMS_ERR_NULL, "at least one tau output pointer is required");
  }
  if (!std::isfinite(t_start) || !std::isfinite(t_end)) {
    return fail(SEAMS_ERR_ARG, "fit window bounds must be finite");
  }
  if (t_start < 0.0 ||
      (!(t_start == 0.0 && t_end == 0.0) && t_end <= t_start)) {
    return fail(SEAMS_ERR_RANGE, "invalid fit time window");
  }
  if (reorient->acc.numLags() < 2) {
    return fail(SEAMS_ERR_STATE, "reorient fit needs at least two lag points");
  }
  try {
    traj::DecayFitResult f1;
    traj::DecayFitResult f2;
    if (!reorient->acc.fit(t_start, t_end, f1, f2)) {
      return fail(SEAMS_ERR_RANGE, "insufficient lag points in fit window");
    }
    if (tau1_int != nullptr) {
      *tau1_int = f1.tauIntegral;
    }
    if (tau1_exp != nullptr) {
      *tau1_exp = f1.tauExp;
    }
    if (tau2_int != nullptr) {
      *tau2_int = f2.tauIntegral;
    }
    if (tau2_exp != nullptr) {
      *tau2_exp = f2.tauExp;
    }
    return SEAMS_OK;
  } catch (const std::bad_alloc &) {
    return fail(SEAMS_ERR_NOMEM, "out of memory");
  } catch (...) {
    remember();
    return SEAMS_ERR_INTERNAL;
  }
}

} // extern "C"
