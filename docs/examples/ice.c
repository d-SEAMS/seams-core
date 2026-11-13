#include "seams.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  const double a = 2.75 * 4.0 / sqrt(3.0);
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
  int labels[8];
  int counts[9];
  seams_system *sys = seams_system_create();
  if (seams_system_set_ortho(sys, xyz, 8, 0, 0, 0, a, a, a, NULL, NULL) != SEAMS_OK) {
    return 1;
  }
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                       counts) != SEAMS_OK) {
    return 1;
  }
  printf("%s cubic %d\n", seams_version_string(), counts[0]);
  seams_system_destroy(sys);
  return counts[0] == 8 ? 0 : 1;
}
