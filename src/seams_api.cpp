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

#include "seams.h"

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <chill_plus_scheme.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <thread_team.hpp>

#include <cmath>
#include <cstring>
#include <exception>
#include <memory>
#include <new>
#include <string>
#include <vector>

#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

thread_local std::string g_error;
int g_threads = 0;

void clear_error() { g_error.clear(); }

seams_status fail(seams_status code, const char *message) {
  g_error = message == nullptr ? "" : message;
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

int unify_type(Cloud &cloud, int type) {
  if (type > 0) {
    return type;
  }
  if (cloud.nop <= 0) {
    return 1;
  }
  const int first = cloud.pts[0].type;
  for (int i = 1; i < cloud.nop; ++i) {
    if (cloud.pts[static_cast<std::size_t>(i)].type != first) {
      for (auto &point : cloud.pts) {
        point.type = 1;
      }
      return 1;
    }
  }
  return first > 0 ? first : 1;
}

std::vector<std::vector<int>> bond_rows(Cloud &cloud, double cutoff, int type,
                                        int graph) {
  const int typ = unify_type(cloud, type);
  if (graph == SEAMS_GRAPH_MUTUAL) {
    return nneigh::kNearestNeighbourList(cloud, 4, cutoff + 1.5, typ, true);
  }
  if (graph != SEAMS_GRAPH_CUTOFF) {
    throw std::invalid_argument("graph is cutoff or mutual");
  }
  return nneigh::neighListO(cutoff, cloud, typ);
}

} // namespace

struct seams_system {
  Cloud cloud;
  bool ready = false;
};

struct seams_neighbours {
  std::vector<int> offsets;
  std::vector<int> cols;
  int n = 0;
};

struct seams_rings {
  std::vector<std::vector<int>> rings;
  std::vector<int> symmetry;
  int n_atoms = 0;
};

extern "C" {

int seams_abi(void) { return SEAMS_ABI; }
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
  if (n < 0) {
    return fail(SEAMS_ERR_RANGE, "thread count is negative");
  }
  g_threads = n;
#ifdef SEAMS_HAS_OPENMP
  if (n > 0) {
    omp_set_num_threads(n);
  }
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

void seams_system_destroy(seams_system *sys) { delete sys; }

seams_status seams_system_set_frame(seams_system *sys, const double *xyz, int n,
                                    const double *box9, const int *types,
                                    const int *ids) {
  clear_error();
  if (sys == nullptr || xyz == nullptr || box9 == nullptr) {
    return fail(SEAMS_ERR_NULL, "system, coordinates, or box is null");
  }
  if (n <= 0) {
    return fail(SEAMS_ERR_RANGE, "frame has no atoms");
  }
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
    cloud.nop = n;
    cloud.currentFrame = 1;
    cloud.pts.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const double x = xyz[3 * i];
      const double y = xyz[3 * i + 1];
      const double z = xyz[3 * i + 2];
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
    Cloud work = sys->cloud;
    auto rows = bond_rows(work, cutoff, type, graph);
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

void seams_neighbours_destroy(seams_neighbours *nb) { delete nb; }

int seams_neighbours_count(const seams_neighbours *nb) {
  return nb == nullptr ? 0 : nb->n;
}

seams_status seams_neighbours_csr(const seams_neighbours *nb, const int **offsets,
                                  int *n_offsets, const int **cols, int *nnz) {
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
    Cloud work = sys->cloud;
    auto rows = bond_rows(work, cutoff, type, graph);
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
    Cloud work = sys->cloud;
    auto rows = bond_rows(work, cutoff, type, SEAMS_GRAPH_CUTOFF);
    const int n = work.nop;
    if (q4 != nullptr) {
      const auto result = chill::steinhardtQl(work, rows, 4);
      if (static_cast<int>(result.ql.size()) < n) {
        return fail(SEAMS_ERR_INTERNAL, "q4 length is short");
      }
      std::memcpy(q4, result.ql.data(), static_cast<std::size_t>(n) * sizeof(double));
    }
    if (q6 != nullptr) {
      const auto result = chill::steinhardtQl(work, rows, 6);
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
    Cloud work = sys->cloud;
    auto rows = bond_rows(work, cutoff, type, graph);
    auto idx = nneigh::neighbourListByIndex(work, rows);
    if (static_cast<int>(idx.size()) < work.nop) {
      idx.resize(static_cast<std::size_t>(work.nop));
    }
    auto found = primitive::ringNetwork(idx, max_depth);
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

void seams_rings_destroy(seams_rings *rings) { delete rings; }

int seams_rings_count(const seams_rings *rings) {
  return rings == nullptr ? 0 : static_cast<int>(rings->rings.size());
}

seams_status seams_rings_histogram(const seams_rings *rings, int *counts,
                                   int counts_len) {
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

} // extern "C"
