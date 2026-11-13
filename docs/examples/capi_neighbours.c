/* One neighbour-list row per Ice Ih oxygen. */
#include "seams.h"

#include <stdio.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_neighbours *nb = NULL;
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  int rows = seams_neighbours_count(nb);
  printf("neighbours %d\n", rows);
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
  return rows == 128 ? 0 : 1;
}
