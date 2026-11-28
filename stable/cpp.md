<!-- https://docs.dseams.info/stable/cpp.md -->

# C++ {#c}

`include/seams.hpp` is the header-only C++ API. It lives in namespace
`seams`, includes `seams.h`, and links `libseams`. Meson installs the
header next to `seams.h`. `pkg-config --cflags seams` and
`find_package(seams)` both put that directory on the include path.
`SEAMS_CXX_HEADER` in `seams-config.cmake` is the installed `seams.hpp`.

## Quickstart {#quickstart}

```bash
c++ -std=c++20 -Iinclude docs/examples/cpp_api.cpp -Lbuild/src -lseams -o cpp_api
./cpp_api "$PWD"
```

The eight-atom diamond is cubic. Ice Ih is hexagonal, 128 atoms, with mean
q6-bar 0.562448843702.

```cpp
/* Copyright (c) d-SEAMS developers
 * Idiomatic seams.hpp on the eight-atom diamond and on ice Ih.
 */
#include "seams.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

std::vector<double> diamond8(double &edge) {
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  std::vector<double> xyz;
  xyz.reserve(24);
  for (const auto &shift : fcc) {
    for (const auto &site : basis) {
      xyz.push_back((shift[0] + site[0]) * a);
      xyz.push_back((shift[1] + site[1]) * a);
      xyz.push_back((shift[2] + site[2]) * a);
    }
  }
  edge = a;
  return xyz;
}

}  // namespace

int main(int argc, char **argv) {
  try {
    seams::set_num_threads(1);
    if (seams::num_threads() != 1 || seams::version().empty()) {
      return 1;
    }
    double edge = 0;
    const auto xyz = diamond8(edge);
    seams::System cell;
    cell.set_ortho(xyz, edge, edge, edge);
    const auto cubic = cell.chill_plus(seams::ChillOptions{});
    const auto orders = cell.q4q6(3.5);
    double mean = 0;
    for (double value : orders.q6) {
      mean += value;
    }
    mean /= static_cast<double>(orders.q6.size());
    if (cubic.counts[0] != 8 || !(mean > 0.3)) {
      return 1;
    }
    std::printf("cpp api cubic %d q6 %.3f %s\n", cubic.counts[0], mean, seams::version().c_str());
    if (argc < 2) {
      return 0;
    }
    seams::System ice;
    ice.read(std::string(argv[1]) + "/input/traj/ice1h.lammpstrj", 1, -1);
    const auto labels = ice.chill_plus(seams::ChillOptions{});
    const auto q6 = ice.q6bar(seams::Q6BarOptions{});
    const auto shell = ice.symmetry(0.0);
    int n_ih = 0;
    for (int label : shell) {
      n_ih += label == 0 ? 1 : 0;
    }
    std::vector<int> flags(static_cast<std::size_t>(ice.count()), 1);
    const auto front = ice.brine(flags, seams::BrineOptions{});
    const auto pockets = front.pockets();
    std::printf("ice Ih atoms %d hexagonal %d q6bar %.12f shell_ih %d pockets %d\n", ice.count(),
                labels.counts[static_cast<std::size_t>(static_cast<int>(seams::Label::Hexagonal))],
                q6.mean, n_ih, pockets.pockets);
    if (ice.count() != 128 || labels.counts[1] != 128) {
      return 1;
    }
    if (std::abs(q6.mean - 0.562448843702) > 1e-9) {
      return 1;
    }
  } catch (const seams::Error &err) {
    std::fprintf(stderr, "seams::Error %d %s\n", err.code(), err.what());
    return 1;
  }
  return 0;
}
```

The wider language ABI remains `libseams_lang` (`seams_frame_*`). A program
that only needs `seams.hpp` does not link it. `docs/examples/minimal_cpp.cpp`
is that ABI on the GenIce sI cell.

```cpp
/** Null when path is null or the file cannot be read.
 * seams_last_error() then names the reason. @ingroup lang_frame */
seams_frame *seams_frame_read(const char *path, int frame, int type_filter);
/** Read with chemfiles name filters and type assignments. options may be
 * null, which is identical to seams_frame_read. @ingroup lang_frame */
seams_frame *seams_frame_read_options(
    const char *path, int frame, int type_filter,
    const seams_chemfiles_load_options *options);
```

Its frame-mutation and order-parameter entries are
`seams_frame_set_positions`, `seams_frame_steinhardt_graph`,
`seams_frame_wl`, `seams_frame_q4_q6`, and
`seams_frame_q4_q6_chill`.

```cpp
/** Replace coordinates of a frame built by seams_frame_from_arrays.
 * n must be the atom count already stored. The box is nine LAMMPS numbers. @ingroup lang_frame */
int seams_frame_set_positions(seams_frame *frame, const double *xyz, int n,
                             const double *box9);

/** Phase times are seconds. reused is 1 when the wide neighbour list was
 * kept from an earlier call. A null times pointer skips the record. */
typedef struct seams_phase_times {
  double neighbours;
  double chill;
  double rings;
  int reused;
} seams_phase_times;

/** CHILL+ and primitive rings on one neighbour list.
 * skin > 0 builds that list at cutoff+skin and keeps it while every atom's
 * minimum-image displacement stays within skin/2. Bonds are the pairs
 * still inside cutoff, so the labels match a fresh cutoff list.
 * skin <= 0 builds a fresh cutoff list on every call. */
/** @ingroup lang_rings */
int seams_frame_chill_rings(seams_frame *frame, double cutoff, int type,
                           const char *graph, int max_ring, double skin,
                           int *labels, int counts9[9], int *ring_counts,
                           int ring_counts_len, seams_phase_times *times);

/** Build the cutoff neighbour list and discard it. graph NULL means cutoff.
 * seams_frame_neighbour_pairs then returns the directed pairs of that list. @ingroup lang_order */
int seams_frame_neighbours(seams_frame *frame, double cutoff, int type,
                           const char *graph);
/** @brief seams_frame_neighbour_pairs. @ingroup lang_order */
int seams_frame_neighbour_pairs(const seams_frame *frame);

/** H-bond edge export. Rows are sorted by (donor, acceptor, hydrogen).
 * donor and acceptor are frame indices. hydrogen is UINT32_MAX when
 * SEAMS_HBOND_EDGE_WITH_H is not requested. fields is the same bitmask as
 * seams_hbond_edges in seams.h. Any non-null output column needs cap rows.
 * n_out receives the edge count. */
int seams_frame_hbond_edges(seams_frame *frame, double cutoff, int oxygen_type,
                            int hydrogen_type, double dist_cutoff,
                            double angle_cutoff, int fields, uint32_t *donor,
                            uint32_t *acceptor, uint32_t *hydrogen,
                            double *d_DA, double *angle, int cap, int *n_out);

/** Local Steinhardt q_l on the cutoff shell. A null qlbar skips the
 * Lechner-Dellago average; the local values do not use that pass.
 * An empty shell is 0: the harmonic row stays zero and is not divided
 * by the neighbour count. freud and pyscal3 store NaN for that atom. */
/** @ingroup lang_order */
int seams_frame_steinhardt(seams_frame *frame, double cutoff, int type,
                           int order_l, double *ql, double *qlbar);
/** graph is "mutual" (four nearest, the historical default) or "cutoff"
 * (the orthoFullCsr list, the same neighbours as `seams steinhardt`). @ingroup lang_order */
int seams_frame_steinhardt_graph(seams_frame *frame, double cutoff, int type,
                                 int order_l, const char *graph, double *ql,
                                 double *qlbar);

/** Normalized Steinhardt hat{w}_l for l = 4 or 6 on the cutoff shell.
 * Returns 1 on failure. @ingroup lang_order */
int seams_frame_wl(seams_frame *frame, double cutoff, int type, int order_l,
                   double *wl);

/** Local q4 and q6 on the cutoff shell. Degrees 4 and 6 use the
 * solid-harmonic kernel. q4bar and q6bar may be null, and a null
 * pointer skips that Lechner-Dellago average. An empty shell is 0
 * for both the local value and the average. Returns 1 on failure. */
/** @ingroup lang_order */
int seams_frame_q4_q6(seams_frame *frame, double cutoff, int type, double *q4,
                      double *q6, double *q4bar, double *q6bar);

/** One cutoff list. Local q4, local q6, Lechner-Dellago q4-bar and
 * q6-bar, and published CHILL+ labels. q4bar and q6bar may be null.
 * counts9 and labels may be null. Returns 1 on failure. @ingroup lang_order */
int seams_frame_q4_q6_chill(seams_frame *frame, double cutoff, int type,
                            double *q4, double *q6, double *q4bar,
                            double *q6bar, int *labels, int counts9[9]);
```

`docs/examples/core_counts.cpp` is the same census on Ice Ih and sII.

```cpp
// Copyright (c) d-SEAMS developers
/* Same Ice Ih and sII census as core_counts.c, through seams.h. */
#include "seams.h"
#include "seams_lammps_box.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static int read_dump(const char *path, std::vector<double> &xyz, std::vector<int> &types,
                     double box9[9]) {
  FILE *in = std::fopen(path, "r");
  char line[512];
  int reading = 0;
  if (in == nullptr) {
    return 1;
  }
  while (std::fgets(line, sizeof line, in) != nullptr) {
    if (std::strncmp(line, "ITEM: BOX BOUNDS", 16) == 0) {
      double bound_lo[3], bound_hi[3], tilt[3];
      for (int axis = 0; axis < 3; ++axis) {
        char bound[128];
        if (std::fgets(bound, sizeof bound, in) == nullptr) {
          std::fclose(in);
          return 1;
        }
        tilt[axis] = 0.0;
        if (std::sscanf(bound, "%lf %lf %lf", &bound_lo[axis], &bound_hi[axis], &tilt[axis]) < 2) {
          std::fclose(in);
          return 1;
        }
      }
      seams_lammps_bounds(bound_lo, bound_hi, tilt, box9);
    } else if (std::strncmp(line, "ITEM: ATOMS", 11) == 0) {
      reading = 1;
    } else if (reading) {
      int id, kind;
      double x, y, z;
      if (std::sscanf(line, "%d %d %lf %lf %lf", &id, &kind, &x, &y, &z) != 5) {
        break;
      }
      types.push_back(kind);
      xyz.push_back(x);
      xyz.push_back(y);
      xyz.push_back(z);
    }
  }
  std::fclose(in);
  return xyz.empty() ? 1 : 0;
}

static int classify(const char *path, int &n, int counts[9], int rings[13]) {
  std::vector<double> xyz;
  std::vector<int> types;
  double box9[9] = {};
  if (read_dump(path, xyz, types, box9) != 0) {
    return 1;
  }
  n = static_cast<int>(types.size());
  std::vector<int> labels(static_cast<size_t>(n));
  seams_system *sys = seams_system_create();
  seams_rings *found = nullptr;
  if (sys == nullptr ||
      seams_system_set_frame(sys, xyz.data(), n, box9, types.data(), nullptr) != SEAMS_OK ||
      seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels.data(),
                       counts) != SEAMS_OK ||
      seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &found) != SEAMS_OK ||
      seams_rings_histogram(found, rings, 13) != SEAMS_OK) {
    if (found != nullptr) {
      seams_rings_destroy(found);
    }
    if (sys != nullptr) {
      seams_system_destroy(sys);
    }
    return 1;
  }
  seams_rings_destroy(found);
  seams_system_destroy(sys);
  return 0;
}

int main() {
  int n = 0;
  int counts[9] = {};
  int rings[13] = {};
  if (seams_set_num_threads(1) != SEAMS_OK) {
    return 1;
  }
  if (classify("input/traj/ice1h.lammpstrj", n, counts, rings) != 0 || n != 128 ||
      counts[1] != 128 || rings[6] != 256) {
    std::fprintf(stderr, "ice Ih n=%d hexagonal=%d rings6=%d\n", n, counts[1], rings[6]);
    return 1;
  }
  if (classify("input/traj/genice_sII.lammpstrj", n, counts, rings) != 0 || n != 136 ||
      counts[4] != 136 || rings[5] != 144 || rings[6] != 16) {
    std::fprintf(stderr, "sII n=%d clathrate=%d rings5=%d rings6=%d\n", n, counts[4], rings[5],
                 rings[6]);
    return 1;
  }
  std::printf("cpp hexagonal 128 rings6 256 clathrate 136 rings5 144\n");
  return 0;
}
```

```cpp
// Copyright (c) d-SEAMS developers
/* Read the GenIce sI cell through the C++ language header. */
#include "seams_lang.h"

#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char **argv) {
  const char *path = argc > 1 ? argv[1] : "input/traj/genice_sI.lammpstrj";
  seams_frame *frame = seams_frame_read(path, 1, 0);
  if (frame == nullptr || seams_frame_nop(frame) != 46) {
    std::fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  if (seams_frame_neighbours(frame, 3.5, 0, "cutoff") != 0) {
    return 1;
  }
  if (seams_frame_neighbour_pairs(frame) < 150) {
    return 1;
  }
  double q6[46];
  if (seams_frame_steinhardt(frame, 3.5, 0, 6, q6, nullptr) != 0) {
    return 1;
  }
  double mean = 0;
  for (double value : q6) {
    mean += value;
  }
  mean /= 46.0;
  int labels[46];
  int counts[9];
  if (seams_frame_chill_plus(frame, 3.5, 0, "cutoff", labels, counts) != 0) {
    return 1;
  }
  int rings[7] = {};
  if (seams_frame_rings(frame, 3.5, 0, "cutoff", 6, rings, 7) != 0) {
    return 1;
  }
  int n_five = 0;
  int per_atom[46];
  if (seams_frame_chordless(frame, 3.5, 0, "cutoff", 5, nullptr, 0, &n_five, per_atom) != 0) {
    return 1;
  }
  int five[48 * 5];
  if (n_five != 48 ||
      seams_frame_chordless(frame, 3.5, 0, "cutoff", 5, five, 48 * 5, &n_five, nullptr) != 0) {
    return 1;
  }
  int through = 0;
  for (int n : per_atom) {
    through += n;
  }
  if (through != 5 * 48 || five[0] != 0) {
    return 1;
  }
  int n_cages = 0;
  int n_sig_atoms = 0;
  int n_vertices = 0;
  if (seams_frame_signature(frame, 3.5, 0, "cutoff", "512", &n_cages, &n_sig_atoms) != 0) {
    return 1;
  }
  int n_edge = 0;
  int n_edge_atoms = 0;
  if (seams_frame_signature_edge(frame, 3.5, 0, "cutoff", "512", "none", &n_edge,
                                 &n_edge_atoms) != 0 ||
      n_edge != n_cages) {
    return 1;
  }
  int n_tight = 0;
  int n_tight_atoms = 0;
  if (seams_frame_signature_edge(frame, 3.5, 0, "cutoff", "512", "1,2", &n_tight,
                                 &n_tight_atoms) != 0 ||
      n_tight != 0) {
    return 1;
  }
  if (n_cages > 0 &&
      seams_frame_signature_vertices(frame, 3.5, 0, "cutoff", "512", 0, nullptr, 0,
                                     &n_vertices) != 0) {
    return 1;
  }
  if (n_cages > 0 &&
      seams_frame_signature_vertices_edge(frame, 3.5, 0, "cutoff", "512", "none", 0, nullptr, 0,
                                          &n_vertices) != 0) {
    return 1;
  }
  // CHILL+ and rings from one neighbour list, with the phase times.
  int fused[9] = {};
  int fused_rings[7] = {};
  seams_phase_times times{};
  if (seams_frame_chill_rings(frame, 3.5, 0, "cutoff", 6, 0.0, labels, fused, fused_rings, 7,
                              &times) != 0 ||
      fused[4] != counts[4]) {
    return 1;
  }
  int shell[46];
  if (seams_frame_symmetry(frame, 0.12, shell) != 0) {
    return 1;
  }
  if (seams_frame_rdf_nbins(6.0, 0) < 1) {
    return 1;
  }
  const int n_r2 = seams_frame_rdf2d_nbin(3.5, 0.5);
  std::vector<double> r2(static_cast<std::size_t>(n_r2));
  std::vector<double> g2(static_cast<std::size_t>(n_r2));
  if (seams_frame_rdf2d(frame, 3.5, 0.5, nullptr, r2.data(), g2.data()) != 0) {
    return 1;
  }
  int sheet[8] = {};
  if (seams_frame_monolayer(frame, 3.5, 0, 6, 1.0, nullptr, sheet, 8) != 0) {
    return 1;
  }
  // The same positions as a new frame with one molecule per atom.
  double box9[9];
  seams_frame_box(frame, box9);
  std::vector<int> types(46, 1);
  std::vector<int> mols(46);
  for (int i = 0; i < 46; ++i) {
    mols[static_cast<std::size_t>(i)] = i + 1;
  }
  seams_frame *copy = seams_frame_from_species(seams_frame_positions(frame), 46, box9,
                                               types.data(), nullptr, mols.data());
  if (copy == nullptr) {
    return 1;
  }
  double tau = 0.0;
  const int jump = seams_frames_jump_tau90(frame, copy, 1.0, 1, 2, &tau);
  seams_frame_free(copy);
  if (jump != 0) {
    return 1;
  }
  if (counts[4] != 46 || rings[5] != 48 || rings[6] != 54 || mean <= 0.5) {
    return 1;
  }
  std::printf("%s atoms 46 clathrate %d q6 %.3f rings %d\n", seams_lang_version(), counts[4],
              mean, rings[5] + rings[6]);
  const double heavy_hb[] = {0, 0, 0, 2.8, 0, 0};
  const double hyd_hb[] = {1, 0, 0};
  const int heavy_type_hb[] = {8, 8};
  const int heavy_mol_hb[] = {1, 2};
  const int hyd_type_hb[] = {2};
  const int hyd_mol_hb[] = {1};
  const double cell_hb[] = {0, 0, 0, 20, 20, 20, 0, 0, 0};
  const int donor_hb[] = {8};
  seams_frame *hb = seams_frame_from_hydrogen_cloud(
      heavy_hb, 2, heavy_type_hb, nullptr, heavy_mol_hb, hyd_hb, 1, hyd_type_hb, nullptr,
      hyd_mol_hb, cell_hb);
  int n_hb = 0;
  if (hb == nullptr ||
      seams_frame_directed_hbonds(hb, 3.5, 2.42, 30.0, donor_hb, 1, donor_hb, 1, 2, 1, nullptr,
                                  nullptr, nullptr, nullptr, nullptr, nullptr, 0, &n_hb) != 0 ||
      n_hb != 1) {
    seams_frame_free(hb);
    return 1;
  }
  seams_frame_free(hb);
  seams_frame_free(frame);
  return 0;
}
```

`docs/examples/ice.cpp` is the eight-atom C call.

```cpp
// Copyright (c) d-SEAMS developers
#include "seams.h"

#include <cmath>
#include <cstdio>

int main() {
  const double a = 2.75 * 4.0 / std::sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  double xyz[24];
  int n = 0;
  for (int f = 0; f < 4; f++) {
    for (int s = 0; s < 2; s++, n++) {
      xyz[3 * n] = (fcc[f][0] + basis[s][0]) * a;
      xyz[3 * n + 1] = (fcc[f][1] + basis[s][1]) * a;
      xyz[3 * n + 2] = (fcc[f][2] + basis[s][2]) * a;
    }
  }
  if (seams_set_num_threads(1) != SEAMS_OK || seams_get_num_threads() != 1) {
    return 1;
  }
  int labels[8];
  int counts[9];
  double q4[8];
  double q6[8];
  seams_system *sys = seams_system_create();
  if (seams_system_set_ortho(sys, xyz, 8, 0, 0, 0, a, a, a, nullptr, nullptr) != SEAMS_OK) {
    return 1;
  }
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts) !=
      SEAMS_OK) {
    return 1;
  }
  if (seams_q4q6(sys, 3.5, 0, q4, q6) != SEAMS_OK) {
    return 1;
  }
  double mean = 0;
  for (double value : q6) {
    mean += value;
  }
  mean /= 8.0;
  std::printf("%s cubic %d q6 %.3f\n", seams_version_string(), counts[0], mean);
  seams_system_destroy(sys);
  return counts[0] == 8 && mean > 0.3 ? 0 : 1;
}
```

## Concepts {#concepts}

`System`, `Neighbours`, `Rings`, `Cages`, `Fingerprint`, and `Brine` are
move-only. The destructor frees the C handle. Copying is deleted.

A status other than success throws `seams::Error`. `what()` is
`seams_handle_error`. `code()` is the status integer. `status()` is the
`Status` enumerator.

Coordinates enter as `std::span<const double>`. Paths and site-table text
enter as `std::string_view`. Results come back as values, `std::vector`,
`std::array`, or `std::string`. `Neighbours::pairs` returns the directed
pairs. `Rings` is a range: `begin` and `end` walk the rings.

`Graph`, `ChillRule`, `Label`, `CageLabel`, `SiteKind`, `RingSymmetry`,
`Axis`, and `DomainSubset` are `enum class`. Analyses with several knobs
take an options aggregate. `ChillOptions`, `SteinhardtOptions`,
`CageOptions`, `DistributionOptions`, `BondOptions`, and the rest default
to the published cutoff of 3.5 angstroms.

Methods that return a value are `[[nodiscard]]`. A method that reads a
frame takes a non-const `System` because the C call clears the handle
error. `count` is const.

A handle requires exclusive use. Two frames on two threads may run
`steinhardt`, `ql`, `q4q6`, or `q6bar` at the same time; the Steinhardt
scratch inside the library belongs to the calling thread. `set_num_threads`
is process-global.

Doxygen comments on `seams.hpp` and `seams.h` are the C++ reference at
`api/index`.

## Ice Ih {#ice-ih}

`input/traj/ice1h.lammpstrj`, frame 1, every atom, cutoff 3.5 angstrom, one
thread. These numbers are the C++ engine on that frame.

```default
ih n 128
ih chill0 0 128 0 0 0 0 0 0 0
ih q3 0.744029547233
ih q4 0.509334928671
ih q6 0.625855471508
ih q6bar 0.562448843702
ih q8 0.229702115021
ih q12 0.421988011385
ih vert 6 0
ih cages 128 0 0 0
ih knn 512
ih rdf 8 4158.664606579 0
ih rcn 26.204724409449
ih cn 4.031496062992
ih hbond 0
ih rho 4 0.030779111111
ih fp 16 256
ih site 8
ih pairs 0 128 0
ih domain 0 0 0
ih prism 0
ih sheet 256
```

`chill0` is the nine original CHILL populations in label order. q3 through
q12 are means. q6bar is the mean of the Lechner-Dellago average. `vert` is
the length of ring 0 and its first atom. `cages` is hexagonal, cubic,
water, both. `knn` is the mutual four-neighbour pair count. `rdf` is the
bin count, the volume in cubic angstroms, and g of the first bin. `rcn` is
the last running coordination. `cn` is the coordination to 3.5 angstroms.
`fp` is the fingerprint key length and the ring-census sum. `site` is atom
0 under `1=waterO`. `pairs` is contacts, cations, and anions under
`1=cationHead`. `domain` is the apolar subset. `prism` is perfect
six-stacks. `sheet` is six-rings in the monolayer histogram.

sI has 46 atoms and 2 cages of signature `512` covering 40 atoms. The XYZ
file has 256 atoms and coordination 43.325490196078 out to 6 angstroms.
The CON file has 4 atoms and 0 hydrogen bonds. Triclinic sH has 136 atoms
and mean q6 0.625601998379.

`symmetry` writes one shell label per atom. `brine` measures a planar front
from an ice flag per atom. `docs/examples/cpp_api.cpp` prints:

```default
ice Ih atoms 128 hexagonal 128 q6bar 0.562448843702 shell_ih 128 pockets 0
```

The test rejects any other hexagonal count or q6-bar mean.

## Errors and ownership {#errors-and-ownership}

```cpp
try {
  seams::System system;
  system.read("missing.lammpstrj", 1, -1);
} catch (const seams::Error &err) {
  // err.what() is seams_handle_error. err.code() is the status.
}
```

`System` and the result handles are move-only. A move leaves the source
empty. The destructor of an empty handle does nothing. The caller does not
call `seams_*_destroy`.

`bindings/cxx/test_lang.cpp` calls `libseams_lang`. That ABI is separate
from `seams.hpp`.

```cpp
ArmFn fault_arm() {
  return reinterpret_cast<ArmFn>(process_symbol("seams_test_arm"));
}
LiveFn fault_live() {
  return reinterpret_cast<LiveFn>(process_symbol("seams_test_live"));
}
LiveFn fault_peak() {
  return reinterpret_cast<LiveFn>(process_symbol("seams_test_peak"));
}
ResetFn fault_reset() {
  return reinterpret_cast<ResetFn>(process_symbol("seams_test_reset_peak"));
}
int main(int argc, char **argv) {
  const std::string root = argc > 1 ? argv[1] : ".";
  const std::string tip4p = root + "/input/traj/exampleTraj.lammpstrj";
  const std::string cubic = root + "/input/traj/cubic-real-7.fixture.lammpstrj";
  const std::string si = root + "/input/traj/genice_sI.lammpstrj";
  const std::string malformedMol2 =
      root + "/tests/fixtures/malformed/mol2_empty_counts.mol2";
  const std::string malformedMmcif =
      root + "/tests/fixtures/malformed/mmcif_no_data.mmcif";

  if (std::string(seams_lang_version()).empty()) {
    return fail("empty version");
  }

  seams_frame *tip = seams_frame_read(tip4p.c_str(), 1, 2);
  if (tip == nullptr) {
    return fail("read tip4p");
  }
  if (seams_frame_nop(tip) != 250) {
    return fail("tip4p nop");
  }
    std::cerr << "cutoff chill+ clathrate " << counts[4] << " inter "
              << counts[5] << " water " << counts[2] << "\n";
    return fail("cutoff chill+ counts");
  }
  int default_counts[9] = {};
  if (seams_frame_chill_plus(tip, 3.5, 2, nullptr, labels.data(),
                             default_counts) != 0) {
    return fail("chill+ default");
  }
  for (int i = 0; i < 9; ++i) {
    if (default_counts[i] != counts[i]) {
      return fail("chill+ default is not the cutoff list");
    }
  }
  if (seams_frame_chill_plus(tip, 3.5, 2, "knn-union", labels.data(),
                             counts) == 0) {
    return fail("chill+ accepted knn-union");
  }

  int n_hex = -1, n_cubic = -1, n_water = -1, n_both = -1;
      ++nfinite;
    }
  }
  mean = nfinite > 0 ? mean / static_cast<double>(nfinite) : 0.0;
  // Orthorhombic frames take the cutoff shell. 0.389362 is the
  // four-nearest mean of the same oxygens.
  if (std::abs(mean - 0.374897) > 1e-4) {
    std::cerr << "q3 mean " << mean << "\n";
    return fail("q3");
  }

  const int bins = 4;
  std::vector<double> r(bins), g(bins);
  std::vector<int> hist(bins);
  double volume = 0.0;
  if (seams_frame_rdf(tip, 6.0, bins, 2, 2, r.data(), g.data(), hist.data(),
                      &volume) != 0) {
    return fail("rdf");
  }
  if (hist[1] != 390 || hist[2] != 705 || hist[3] != 531) {
  if (seams_frame_chill_plus(ice, 3.5, 1, "mutual", ice_labels.data(),
                             ice_counts) != 0) {
    return fail("mW chill+");
  }
  if (ice_counts[0] != 4096) {
    std::cerr << "mW cubic " << ice_counts[0] << "\n";
    return fail("mW cubic count");
  }
  if (seams_frame_cages(ice, 3.5, 1, 4, "seeded", 0, nullptr, &n_hex, &n_cubic,
                        &n_water, &n_both) != 0) {
    return fail("mW cages");
  }
  if (n_hex != 0 || n_cubic != 4096 || n_water != 0) {
    std::cerr << "mW cages " << n_hex << " " << n_cubic << " " << n_water
              << "\n";
    return fail("mW cage counts");
  }
  if (seams_path_nframes(cubic.c_str()) != 7) {
    return fail("mW nframes");
  }

  seams_frame *guest = seams_frame_read(si.c_str(), 1, 1);
  if (guest == nullptr) {
    return fail("sI");
  }
  int n_cages = 0;
  int n_atoms = 0;
      return fail("abi disagrees with the cli");
    }
  }

  double box9[9] = {};
  seams_frame_box(tip, box9);
  if (!(box9[3] > 0.0)) {
    return fail("box");
  }
  if (seams_frame_index(tip) < 1) {
    return fail("index");
  }
  if (seams_frame_types(tip) == nullptr || seams_frame_ids(tip) == nullptr) {
    return fail("arrays");
  }
  if (seams_frame_positions(nullptr) != nullptr) {
    return fail("null positions");
  }
  if (seams_frame_chill(tip, 3.5, 2, "cutoff", labels.data(), counts) != 0) {
    return fail("plain chill");
  }
  if (seams_frame_chill_plus_rule(tip, 3.5, 2, "cutoff", "ovito-chill+", labels.data(), counts) != 0) {
    return fail("ovito rule");
  }
  if (seams_frame_chill_plus_rule(tip, 3.5, 2, "cutoff", "not-a-rule", labels.data(), counts) == 0) {
    return fail("bad rule accepted");
  }
  if (seams_frame_from_arrays(nullptr, 1, box9, nullptr, nullptr) != nullptr) {
    return fail("null xyz");
  }
  seams_frame *built = seams_frame_from_arrays(xyz, 4, box9, nullptr, nullptr);
  if (built == nullptr || seams_frame_nop(built) != 4) {
    return fail("from arrays");
  }
  seams_frame_free(built);

  const int graphOff[7] = {0, 1, 2, 3, 4, 5, 6};
  const int graphCol[6] = {1, 2, 3, 4, 5, 0};
  int graphHist[7] = {};
  int nGraphRings = 0;
  int nGraphVertices = 0;
  if (seams_graph_rings(graphOff, 7, graphCol, 6, 6, graphHist, 7,
                        nullptr, 0, nullptr, 0, &nGraphRings,
                        &nGraphVertices) != 0 ||
      graphHist[6] != 1 || nGraphRings != 1 || nGraphVertices != 6) {
    return fail("explicit graph ring query");
  }
  int packedOff[2] = {};
  int packedVertices[6] = {};
  if (seams_graph_rings(graphOff, 7, graphCol, 6, 6, graphHist, 7,
                        packedOff, 2, packedVertices, 6, &nGraphRings,
                        &nGraphVertices) != 0 ||
      packedOff[1] != 6) {
    return fail("explicit graph rings");
  }

  const double tetraXyz[18] = {
      50, 50, 50, 51, 51, 51, 51, 49, 49,
      49, 51, 49, 49, 49, 51, 50, 50, 52,
  };
  const double analyticBox[9] = {0, 0, 0, 100, 100, 100, 0, 0, 0};
  seams_frame *analytic =
      seams_frame_from_arrays(tetraXyz, 6, analyticBox, nullptr, nullptr);
  const int centre[1] = {0};
  const int tetraPartners[4] = {1, 2, 3, 4};
  double tetraQ[1] = {};
  if (analytic == nullptr ||
      seams_frame_tetrahedral(analytic, centre, 1, tetraPartners, 4, 1,
                              2.0, tetraQ) != 0 ||
      std::abs(tetraQ[0] - 1.0) > 1e-12) {
    return fail("language tetrahedral order");
  }
  seams_frame_free(analytic);

  const double radialXyz[18] = {
      50, 50, 50, 52.8, 50, 50, 50, 52.9, 50,
      50, 50, 53.1, 46.5, 50, 50, 50, 45.8, 50,
  };
  const int radialPartners[4] = {1, 2, 3, 4};
  analytic =
      seams_frame_from_arrays(radialXyz, 6, analyticBox, nullptr, nullptr);
  double lsi[1] = {};
  if (analytic == nullptr ||
      seams_frame_lsi(analytic, centre, 1, nullptr, 0, 3.7, lsi) != 0 ||
      std::abs(lsi[0] - 0.0525) > 1e-12) {
    return fail("language LSI");
  }
  if (seams_frame_lsi(analytic, centre, 1, radialPartners, 4, 3.7, lsi) !=
          0 ||
      !std::isnan(lsi[0])) {
    return fail("language LSI without a partner beyond the cutoff");
  }
  seams_frame_free(analytic);

  // A triangle at 1, 4, 7 in a cell of 9 steps +3 three times: one +L.
  const double windingBox[9] = {0, 0, 0, 9, 9, 9, 0, 0, 0};
  const int triangleOff[4] = {0, 1, 2, 3};
  const int triangleCol[3] = {1, 2, 0};
  int triangleHist[4] = {};
  int trianglePackedOff[2] = {};
  int triangleVertices[3] = {};
  int triangleRings = 0;
  int triangleVertexCount = 0;
  if (seams_graph_rings(triangleOff, 4, triangleCol, 3, 3, triangleHist, 4,
                        trianglePackedOff, 2, triangleVertices, 3,
                        &triangleRings, &triangleVertexCount) != 0 ||
      triangleRings != 1 || triangleVertices[0] != 0 ||
      triangleVertices[1] != 1 || triangleVertices[2] != 2) {
    return fail("language winding ring");
  }
  struct WindingCase {
    double xyz[9];
    int winding[3];
    int contractible;
  };
  const WindingCase windingCases[4] = {
      {{1, 4, 4, 4, 4, 4, 7, 4, 4}, {1, 0, 0}, 0},
      {{7, 4, 4, 4, 4, 4, 1, 4, 4}, {-1, 0, 0}, 0},
      {{4, 1, 4, 4, 4, 4, 4, 7, 4}, {0, 1, 0}, 0},
      {{2, 4, 4, 4, 4, 4, 3, 4, 4}, {0, 0, 0}, 1},
  };
  for (const auto &windingCase : windingCases) {
    analytic = seams_frame_from_arrays(windingCase.xyz, 3, windingBox,
                                       nullptr, nullptr);
    if (analytic == nullptr) {
      return fail("language winding frame");
    }
    for (int axis = 0; axis < 3; ++axis) {
      int axisWinding[1] = {99};
      int isContractible[1] = {99};
      if (seams_frame_ring_winding(analytic, trianglePackedOff, triangleRings,
                                   triangleVertices, triangleVertexCount,
                                   axis, axisWinding, isContractible) != 0 ||
          axisWinding[0] != windingCase.winding[axis] ||
          isContractible[0] != windingCase.contractible) {
        return fail("language ring winding");
      }
    }
    seams_frame_free(analytic);
  }

  if (seams_frame_prisms(nullptr, 3.5, 2, 6, nullptr, nullptr, nullptr) == 0) {
    return fail("null prisms");
  }
  int n_rows = 0;
  if (seams_frame_prisms(tip, 3.5, 2, 6, nullptr, nullptr, &n_rows) != 0 ||
      n_rows != 4) return fail("prism count");
  std::vector<seams_prism_row> prism_rows(static_cast<std::size_t>(n_rows));
  if (seams_frame_prisms(tip, 3.5, 2, 6, nullptr, prism_rows.data(), &n_rows) != 0) {
    return fail("prisms");
  }
  const int nbin = seams_frame_rdf2d_nbin(3.5, 0.0);
  std::vector<double> r2(static_cast<std::size_t>(nbin));
  std::vector<double> g2(static_cast<std::size_t>(nbin));
  if (seams_frame_rdf2d(nullptr, 3.5, 0.5, nullptr, r2.data(), g2.data()) == 0) {
    return fail("null rdf2d");
  }
  if (seams_frame_rdf2d(tip, 3.5, 0.5, nullptr, r2.data(), g2.data()) != 0) {
    return fail("rdf2d");
  }
  seams_frame_box(typed, filled);
  if (!(filled[6] != 0.0)) {
    return fail("shear");
  }
  seams_frame_free(typed);
  seams_frame_free(nullptr);
  if (seams_frame_nop(nullptr) != 0) {
  if (seams_frame_cages(tip, 3.5, 2, 4, "cutoff", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) != 0) {
    return fail("cutoff cages");
  }
  if (seams_frame_cages(nullptr, 3.5, 2, 4, "cutoff", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) == 0) {
    return fail("null cages");
  }
  if (seams_frame_cages(tip, 3.5, 2, 4, "banana", 0, nullptr, &cageHex, nullptr, nullptr, nullptr) == 0) {
    return fail("banana cages");
  if (seams_frame_steinhardt(nullptr, 3.5, 2, 4, ql.data(), qlbar.data()) == 0) {
    return fail("null steinhardt");
  }
    return fail("prism both null");
  }
  const std::string withSlash = outDir + "/";
  if (seams_frame_rdf2d(oxygen, 3.5, 0.5, withSlash.c_str(), rr, gg) != 0) {
        seams_frame_rdf2d(plane, 6.0, 0.5, nullptr, planeR, planeG) != 0) {
      seams_frame_free(plane);
      return fail("coplanar rdf2d");
    }
    for (double value : planeG) {
      if (!std::isfinite(value)) {
        seams_frame_free(plane);
        return fail("coplanar rdf2d values");
      }
    }
    seams_frame_free(plane);
  }
  seams_frame_free(blank2);
  seams_frame_free(oxygen);
  std::filesystem::remove_all(out);

  seams_frame_free(molten);
  seams_frame_free(tip);
  seams_frame_free(ice);
  seams_frame_free(guest);
  {
    const double brineXyz[] = {10, 10, 10, 12.8, 10, 10, 9.07, 12.42, 10,
                               9.07, 8.79, 11.62, 9.07, 8.79, 8.38, 12.0, 10, 10};
    const double brineBox[] = {0, 0, 0, 30, 30, 30, 0, 0, 0};
    const int types[] = {1, 1, 1, 1, 1, 3};
    const int mols[] = {1, 2, 3, 4, 5, 6};
    seams_frame *mixed =
        seams_frame_from_species(brineXyz, 6, brineBox, types, nullptr, mols);
    if (mixed == nullptr) {
      return fail("brine frame");
    }
    int mixedLabels[6] = {};
    if (seams_frame_chill_plus(mixed, 3.5, 0, "cutoff", mixedLabels, nullptr) !=
        0) {
      seams_frame_free(mixed);
      return fail("brine untyped chill");
    }
    if (mixedLabels[5] == 6) {
      seams_frame_free(mixed);
      return fail("brine ion was already unclassified");
    }
    if (seams_frame_select(mixed, 1, nullptr, 0) != 0) {
      seams_frame_free(mixed);
      return fail("brine select");
    }
    int selected[6] = {};
    if (seams_frame_chill_plus(mixed, 3.5, 1, "cutoff", selected, nullptr) != 0) {
      seams_frame_free(mixed);
      return fail("brine selected chill");
    }
    if (selected[5] != 6) {
      seams_frame_free(mixed);
      return fail("brine ion stayed labeled");
    }
    int ringCounts[7] = {};
    if (seams_frame_rings(mixed, 3.5, 1, "cutoff", 6, ringCounts, 7) != 0) {
      seams_frame_free(mixed);
      return fail("brine rings");
    }
    int hex = 0;
    int cubicCages = 0;
    int water = 0;
    int both = 0;
    if (seams_frame_cages(mixed, 3.5, 1, 4, "cutoff", 0, nullptr, &hex, &cubicCages,
                          &water, &both) != 0) {
      seams_frame_free(mixed);
      return fail("brine cages");
    }
    if (hex + cubicCages + water + both != 5) {
      seams_frame_free(mixed);
      return fail("brine cages counted the ion");
    }
    if (seams_frame_select(mixed, 9, nullptr, 0) == 0) {
      seams_frame_free(mixed);
      return fail("brine missing type was accepted");
    }
    if (std::string(seams_last_error()) != "atom type 9 is not in the frame") {
      seams_frame_free(mixed);
      return fail("brine missing type text");
    }
    const int absentMol = 99;
    if (seams_frame_select(mixed, 1, &absentMol, 1) == 0) {
      seams_frame_free(mixed);
      return fail("brine missing molecule was accepted");
    }
    seams_frame_free(mixed);
  }
  {
    seams_frame *fresh = seams_frame_read(tip4p.c_str(), 1, 2);
    if (fresh == nullptr) {
      return fail("reread tip4p");
    }
    const int nFresh = seams_frame_nop(fresh);
    std::vector<double> moved(static_cast<std::size_t>(nFresh) * 3);
    const double *src = seams_frame_positions(fresh);
    for (int i = 0; i < nFresh * 3; ++i) {
      moved[static_cast<std::size_t>(i)] = src[i];
    }
    moved[0] += 1.0e-4;
    double freshBox[9];
    seams_frame_box(fresh, freshBox);
    if (seams_frame_set_positions(fresh, moved.data(), nFresh, freshBox) != 0) {
      seams_frame_free(fresh);
      return fail("set positions");
    }
    if (seams_frame_neighbours(fresh, 3.5, 2, "cutoff") != 0) {
      seams_frame_free(fresh);
      return fail("neighbours");
    }
    if (seams_frame_neighbour_pairs(fresh) <= 0) {
      seams_frame_free(fresh);
      return fail("neighbour pairs");
    }
    std::vector<double> freshQl(static_cast<std::size_t>(nFresh));
    std::vector<double> bar(static_cast<std::size_t>(nFresh));
    if (seams_frame_steinhardt(fresh, 3.5, 2, 6, freshQl.data(), bar.data()) != 0) {
      seams_frame_free(fresh);
      return fail("steinhardt");
    }
    if (seams_frame_steinhardt_graph(fresh, 3.5, 2, 6, "cutoff", freshQl.data(), bar.data()) != 0) {
      seams_frame_free(fresh);
      return fail("steinhardt graph");
    }
    std::vector<double> wl(static_cast<std::size_t>(nFresh));
    if (seams_frame_wl(fresh, 3.5, 2, 4, wl.data()) != 0) {
      seams_frame_free(fresh);
      return fail("wl");
    }
    std::vector<double> q4(static_cast<std::size_t>(nFresh));
    std::vector<double> q6(static_cast<std::size_t>(nFresh));
    if (seams_frame_q4_q6(fresh, 3.5, 2, q4.data(), q6.data(), nullptr, nullptr) != 0) {
      seams_frame_free(fresh);
      return fail("q4 q6");
    }
    std::vector<int> iceLabels(static_cast<std::size_t>(nFresh));
    int fused[9] = {};
    if (seams_frame_q4_q6_chill(fresh, 3.5, 2, q4.data(), q6.data(), nullptr, nullptr,
                               iceLabels.data(), fused) != 0) {
      seams_frame_free(fresh);
      return fail("q4 q6 chill");
    }
    std::vector<int> ringCounts(9);
    seams_phase_times times{};
    if (seams_frame_chill_rings(fresh, 3.5, 2, "cutoff", 6, 0.0, iceLabels.data(), fused,
                               ringCounts.data(), 9, &times) != 0) {
      seams_frame_free(fresh);
      return fail("chill rings");
    }
```

## Symbol reference {#symbol-reference}

### seams_frame {#seams-frame}

Purpose: Opaque frame used by the language bindings.

Parameters: none

Returns: The type itself.

Ownership: No allocation. This is a type.

Errors: None.

Units: handle

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_box {#seams-frame-box}

Purpose: Copy the nine LAMMPS cell numbers.

Parameters: frame; box9

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_cages {#seams-frame-cages}

Purpose: Score hexagonal and double-diamond cage membership.

Parameters: frame; cutoff; type; k; graph; complete; per_atom; n_hex; n_cubic; n_water; n_both

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; k is a dimensionless neighbour count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill {#seams-frame-chill}

Purpose: Write one original CHILL label per atom.

Parameters: frame; cutoff; type; graph; labels; counts9

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill_plus {#seams-frame-chill-plus}

Purpose: Write one CHILL+ label per atom.

Parameters: frame; cutoff; type; graph; labels; counts9

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill_plus_rule {#seams-frame-chill-plus-rule}

Purpose: Write CHILL+ labels for a named rule.

Parameters: frame; cutoff; type; graph; rule; labels; counts9

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; rule is text

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill_rings {#seams-frame-chill-rings}

Purpose: Run CHILL+ and primitive rings on one neighbour list.

Parameters: frame; cutoff; type; graph; max_ring; skin; labels; counts9; ring_counts; ring_counts_len; times

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; skin in angstrom; phase times in seconds

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chordless {#seams-frame-chordless}

Purpose: Chordless cycles of one size on a frame.

Parameters: frame; cutoff; type; graph; size; rings; rings_len; n_rings; per_atom

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; size and indices are dimensionless, indices 0-based

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_free {#seams-frame-free}

Purpose: Free a language-binding frame.

Parameters: frame

Returns: Nothing.

Ownership: Frees the caller-owned frame. A null pointer is a no-op.

Errors: None.

Units: handle

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_from_arrays {#seams-frame-from-arrays}

Purpose: Build a frame from coordinates and a cell.

Parameters: xyz; n; box9; types; ids

Returns: A new seams_frame, or NULL when the call cannot build one.

Ownership: The caller owns the frame and frees it with seams_frame_free.

Errors: NULL. seams_last_error() names the reason.

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_from_species {#seams-frame-from-species}

Purpose: Build a frame and keep molecule ids.

Parameters: xyz; n; box9; types; ids; mols

Returns: A new seams_frame, or NULL when the call cannot build one.

Ownership: The caller owns the frame and frees it with seams_frame_free.

Errors: NULL. seams_last_error() names the reason.

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_ids {#seams-frame-ids}

Purpose: Pointer to the atom ids stored on the frame.

Parameters: frame

Returns: A pointer into the frame.

Ownership: The frame owns the memory. It stays valid until seams_frame_free or a call that replaces the coordinates.

Errors: NULL when the frame is null.

Units: dimensionless ids

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_index {#seams-frame-index}

Purpose: Return the frame index that was read.

Parameters: frame

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: dimensionless index, 1-based in the file

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_monolayer {#seams-frame-monolayer}

Purpose: Count rings in sheets of a given area.

Parameters: frame; cutoff; type; max_depth; sheet_area; output_dir; counts; counts_len

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; sheet_area in square angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_neighbour_pairs {#seams-frame-neighbour-pairs}

Purpose: Return the number of directed neighbour pairs.

Parameters: frame

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_neighbours {#seams-frame-neighbours}

Purpose: Build the neighbour list and discard the rows.

Parameters: frame; cutoff; type; graph

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_nop {#seams-frame-nop}

Purpose: Return the number of atoms in the frame.

Parameters: frame

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_positions {#seams-frame-positions}

Purpose: Pointer to the coordinates stored on the frame.

Parameters: frame

Returns: A pointer into the frame.

Ownership: The frame owns the memory. It stays valid until seams_frame_free or a call that replaces the coordinates.

Errors: NULL when the frame is null.

Units: angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_prisms {#seams-frame-prisms}

Purpose: Count perfect and imperfect prism stacks.

Parameters: frame; cutoff; type; max_depth; output_dir; rows; n_rows

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_q4_q6 {#seams-frame-q4-q6}

Purpose: Write local and averaged q4 and q6.

Parameters: frame; cutoff; type; q4; q6; q4bar; q6bar

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; order parameters are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_q4_q6_chill {#seams-frame-q4-q6-chill}

Purpose: Write q4, q6, and CHILL+ from one cutoff list.

Parameters: frame; cutoff; type; q4; q6; q4bar; q6bar; labels; counts9

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; order parameters are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf {#seams-frame-rdf}

Purpose: Radial distribution of one type pair.

Parameters: frame; rmax; nbins; type_i; type_j; r; g; count; volume

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: rmax and r in angstrom; g is dimensionless; volume in cubic angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf2d {#seams-frame-rdf2d}

Purpose: Two-dimensional radial distribution in a monolayer.

Parameters: frame; cutoff; binwidth; output_dir; r; g

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff, binwidth, and r in angstrom; g is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf2d_nbin {#seams-frame-rdf2d-nbin}

Purpose: Number of bins in the monolayer RDF.

Parameters: cutoff; binwidth

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff and binwidth in angstrom; result is a dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf_nbins {#seams-frame-rdf-nbins}

Purpose: Number of bins in the bulk RDF.

Parameters: rmax; nbins

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: rmax in angstrom; nbins is a dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_read {#seams-frame-read}

Purpose: Read one frame from a trajectory path.

Parameters: path; frame; type_filter

Returns: A new seams_frame, or NULL when the call cannot build one.

Ownership: The caller owns the frame and frees it with seams_frame_free.

Errors: NULL. seams_last_error() names the reason.

Units: frame index is dimensionless; coordinates in the file are angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rings {#seams-frame-rings}

Purpose: Histogram of primitive rings.

Parameters: frame; cutoff; type; graph; max_depth; counts; counts_len

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; counts\[k\] is a dimensionless ring count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_select {#seams-frame-select}

Purpose: Keep selected oxygens for later classifiers.

Parameters: frame; type; mols; nmol

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: type and molecule ids are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_set_positions {#seams-frame-set-positions}

Purpose: Replace coordinates on an array-built frame.

Parameters: frame; xyz; n; box9

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_signature {#seams-frame-signature}

Purpose: Count closed cages whose ring census matches a spec.

Parameters: frame; cutoff; type; graph; spec; n_cages; n_atoms

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; n_cages and n_atoms are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_signature_edge {#seams-frame-signature-edge}

Purpose: Same cage census with an explicit face-edge window.

Parameters: frame; cutoff; type; graph; spec; edge; n_cages; n_atoms

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; edge lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_signature_vertices {#seams-frame-signature-vertices}

Purpose: Atom indices of one matched cage, sorted.

Parameters: frame; cutoff; type; graph; spec; cage_index; vertices; cap; n_vertices

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner. vertices may be null when the call only asks for the length.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; cage_index, cap, and n_vertices are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_signature_vertices_edge {#seams-frame-signature-vertices-edge}

Purpose: Atom indices of one matched cage on an explicit face-edge window.

Parameters: frame; cutoff; type; graph; spec; edge; cage_index; vertices; cap; n_vertices

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner. vertices may be null when the call only asks for the length. edge is read during the call.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff and edge lengths in angstrom; cage_index, cap, and n_vertices are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_steinhardt {#seams-frame-steinhardt}

Purpose: Local Steinhardt q_l on the cutoff shell.

Parameters: frame; cutoff; type; order_l; ql; qlbar

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; ql is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_steinhardt_graph {#seams-frame-steinhardt-graph}

Purpose: Local Steinhardt q_l on a named graph.

Parameters: frame; cutoff; type; order_l; graph; ql; qlbar

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; ql is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_symmetry {#seams-frame-symmetry}

Purpose: Per-atom shell label for Ih, Ic, or stacking disorder.

Parameters: frame; shell_rmsd; labels

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: shell_rmsd is a dimensionless residual; labels are dimensionless codes

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_types {#seams-frame-types}

Purpose: Pointer to the atom types stored on the frame.

Parameters: frame

Returns: A pointer into the frame.

Ownership: The frame owns the memory. It stays valid until seams_frame_free or a call that replaces the coordinates.

Errors: NULL when the frame is null.

Units: dimensionless types

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_wl {#seams-frame-wl}

Purpose: Normalized Steinhardt w_l for degree 4 or 6.

Parameters: frame; cutoff; type; order_l; wl

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom; coordinates in angstrom; wl is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frames_jump_tau90 {#seams-frames-jump-tau90}

Purpose: Time for an H-H vector to rotate 90 degrees.

Parameters: first; second; dt; oxygen_type; hydrogen_type; tau

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: dt and tau in the caller’s time unit; types are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_lang_version {#seams-lang-version}

Purpose: Version string of the language-binding library.

Parameters: none

Returns: A pointer to library text.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/minimal_cpp.cpp`

### seams_last_error {#seams-last-error}

Purpose: Return the thread-local error detail.

Parameters: none

Returns: A pointer to library text.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/minimal_cpp.cpp`

### seams_path_nframes {#seams-path-nframes}

Purpose: Count frames in a trajectory path.

Parameters: path

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_phase_times {#seams-phase-times}

Purpose: Seconds spent in neighbours, CHILL+, and rings.

Parameters: neighbours; chill; rings; reused

Returns: The type itself.

Ownership: No allocation. This is a type.

Errors: None.

Units: times in seconds; reused is 0 or 1

Example: `docs/examples/minimal_cpp.cpp`

### seams_prism_row {#seams-prism-row}

Purpose: One row of a prism stack: size, perfect, imperfect.

Parameters: size; perfect; imperfect

Returns: The type itself.

Ownership: No allocation. This is a type.

Errors: None.

Units: dimensionless counts

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_directed_hbonds {#seams-frame-directed-hbonds}

Purpose: Directed hydrogen-bond edges sorted by donor, hydrogen, and acceptor.

Parameters: frame; pair_cutoff; dist_cutoff; angle_cutoff; donor_types; n_donor_types; acceptor_types; n_acceptor_types; hydrogen_type; exclude_same_mol; donor; hydrogen; acceptor; d_AH; d_DA; angle; cap; n_out

Returns: Zero on success. A non-zero value reports the failure.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: donor, hydrogen, and acceptor are 0-based indices; distances in angstrom; angle in degrees

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_from_hydrogen_cloud {#seams-frame-from-hydrogen-cloud}

Purpose: Build a frame from a heavy-atom cloud and a separate hydrogen cloud.

Parameters: heavy_xyz; n_heavy; heavy_types; heavy_ids; heavy_mols; hydrogen_xyz; n_hydrogen; hydrogen_types; hydrogen_ids; hydrogen_mols; box9

Returns: A new seams_frame, or NULL when the call cannot build one.

Ownership: The caller owns the frame and frees it with seams_frame_free.

Errors: NULL. seams_last_error() names the reason.

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_lsi {#seams-frame-lsi}

Purpose: Write Local Structure Index for selected centres and partners.

Parameters: frame; centres; n_centres; partners; n_partners; cutoff; lsi

Returns: Zero on success. A non-zero value reports failure.

Ownership: The caller owns every buffer. The frame stays with its owner.

Errors: A non-zero code; seams_last_error() holds the detail.

Units: cutoff in angstrom; LSI in square angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_ring_winding {#seams-frame-ring-winding}

Purpose: Write periodic winding and contractibility for packed rings.

Parameters: frame; ring_offsets; n_rings; vertices; n_vertices; axis; axis_winding; contractible

Returns: Zero on success. A non-zero value reports failure.

Ownership: The caller owns every buffer. The frame stays with its owner.

Errors: A non-zero code; seams_last_error() holds the detail.

Units: dimensionless lattice winding and flags

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_tetrahedral {#seams-frame-tetrahedral}

Purpose: Write tetrahedral order for selected centres and partners.

Parameters: frame; centres; n_centres; partners; n_partners; convention; cutoff; q

Returns: Zero on success. A non-zero value reports failure.

Ownership: The caller owns every buffer. The frame stays with its owner.

Errors: A non-zero code; seams_last_error() holds the detail.

Units: cutoff in angstrom; q is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_graph_rings {#seams-graph-rings}

Purpose: Find primitive rings in a caller-supplied CSR graph.

Parameters: offsets; n_offsets; columns; n_columns; max_depth; histogram; histogram_len; ring_offsets; ring_offsets_cap; vertices; vertices_cap; n_rings; n_vertices

Returns: Zero on success. A non-zero value reports failure.

Ownership: The caller owns every input and output buffer.

Errors: A non-zero code; seams_last_error() holds the detail.

Units: dimensionless indices and counts

Example: `docs/examples/minimal_cpp.cpp`

### SEAMS_HBOND_EDGE_WITH_H {#seams-hbond-edge-with-h}

Purpose: SEAMS_HBOND_EDGE_WITH_H

Parameters: none

Returns: The constant value.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_cpp.cpp`

### seams_chemfiles_load_options {#seams-chemfiles-load-options}

Purpose: Name filters and explicit type assignments for chemfiles-backed input. Declared in bindings/cxx/seams_lang.h.

Parameters: resnames; resname_count; excluded_resnames; excluded_resname_count; atomnames; atomname_count; excluded_atomnames; excluded_atomname_count; type_map; type_map_count

Returns: The type itself.

Ownership: No allocation. This is a type.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_hbond_edges {#seams-frame-hbond-edges}

Purpose: H-bond edge export. Declared in bindings/cxx/seams_lang.h.

Parameters: frame; cutoff; oxygen_type; hydrogen_type; dist_cutoff; angle_cutoff; fields; donor; acceptor; hydrogen; d_DA; angle; cap; n_out

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: cutoff in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_read_options {#seams-frame-read-options}

Purpose: Read with chemfiles name filters and type assignments. Declared in bindings/cxx/seams_lang.h.

Parameters: path; frame; type_filter; options

Returns: A new seams_frame, or NULL when the call cannot build one.

Ownership: The caller owns the frame and frees it with seams_frame_free.

Errors: NULL. seams_last_error() names the reason.

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_hbond_edges {#seams-hbond-edges}

Purpose: Seams_hbond_edges

Parameters: none

Returns: Zero on success. A non-zero value is the seams_status code.

Ownership: The caller owns every buffer it passed. The frame stays with its owner.

Errors: A non-zero code. seams_last_error() holds the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_cpp.cpp`

## Ice Ih {#id1}

![Ice Ih colored by published CHILL+](../_images/overlay-web15.png)

*Ice Ih colored by published CHILL+.*
