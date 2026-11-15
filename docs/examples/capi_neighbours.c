/* One neighbour-list row per Ice Ih oxygen. */
#include "seams.h"

#include <stdio.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_neighbours *nb = NULL;
  size_t nrows = 0;
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  int rows = seams_neighbours_count(nb);
  if (seams_neighbours_nrows(nb, &nrows) != SEAMS_OK || nrows != 128 || rows != 128) {
    fprintf(stderr, "nrows %zu count %d\n", nrows, rows);
    return 1;
  }
  printf("neighbours %d\n", rows);
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
  return 0;
}
