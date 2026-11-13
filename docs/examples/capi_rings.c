/* Ice Ih has 256 primitive six-rings at 3.5 angstroms. */
#include "seams.h"

#include <stdio.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_rings *rings = NULL;
  int hist[13] = {0};
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK ||
      seams_rings_histogram(rings, hist, 13) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  printf("rings6 %d\n", hist[6]);
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return hist[6] == 256 ? 0 : 1;
}
