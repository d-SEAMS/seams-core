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
