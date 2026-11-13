#include "seams.h"
#include <math.h>
#include <stdio.h>
int main(void) {
  const double a = 2.75 * 4.0 / sqrt(3.0);
  const double s[8][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5},
                          {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  double xyz[24];
  for (int i = 0; i < 8; i++) {
    xyz[3 * i] = s[i][0] * a;
    xyz[3 * i + 1] = s[i][1] * a;
    xyz[3 * i + 2] = s[i][2] * a;
  }
  int labels[8], counts[9];
  seams_system *sys = seams_system_create();
  if (seams_system_set_ortho(sys, xyz, 8, 0, 0, 0, a, a, a, 0, 0) != SEAMS_OK) return 1;
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts) != SEAMS_OK)
    return 1;
  printf("cubic %d\n", counts[0]);
  seams_system_destroy(sys);
  return counts[0] == 8 ? 0 : 1;
}
