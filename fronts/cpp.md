<!-- https://docs.dseams.info/fronts/cpp.md -->

# C++ tutorial {#c-tutorial}

Classify the eight-atom cubic cell with `include/seams.hpp`, then the ice Ih dump. The program is `docs/examples/cpp_api.cpp`. Passed the repository root, the first line is `cpp api cubic 8 q6 0.629` followed by the project version. The second line is:

```default
ice Ih atoms 128 hexagonal 128 q6bar 0.562448843702 shell_ih 128 pockets 0
```

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

`seams.hpp` includes `seams.h`, which includes the generated `seams_version.h`. Meson compiles the program against `libseams`.

```bash
c++ -std=c++20 -Iinclude -Ibuild/src docs/examples/cpp_api.cpp -Lbuild/src -lseams -Wl,-rpath,build/src -o cpp_api
./cpp_api "$PWD"
```

`seams::set_num_threads(0)` restores the OpenMP width captured at the first call, `omp_get_max_threads()`. A type of 0 or less on a frame whose atoms do not all share one positive type classifies a copy in which those atoms are given one type. That copy is `unified_copy` in `src/seams_api.cpp`.

`docs/examples/core_counts.cpp` checks the same Ice Ih and sII census. It prints `cpp hexagonal 128 rings6 256 clathrate 136 rings5 144` and returns 0 only when sII also has 16 six-rings.

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

The class list is the [C++ API](../stable/cpp.md).

## Ice Ih {#ice-ih}

![Ice Ih colored by published CHILL+](../_images/overlay-web15.png)

*Ice Ih colored by published CHILL+.*
