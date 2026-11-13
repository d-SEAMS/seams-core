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
