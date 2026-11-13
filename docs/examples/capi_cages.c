/* Every Ice Ih oxygen is in a hexagonal cage and not a cubic one. */
#include "seams.h"

#include <stdio.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_cage_set *cages = NULL;
  int hex = 0;
  int cubic = 0;
  int water = 0;
  int both = 0;
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 4, 0, &cages) != SEAMS_OK ||
      seams_cages_counts(cages, &hex, &cubic, &water, &both) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  printf("cages %d %d %d %d\n", hex, cubic, water, both);
  seams_cages_destroy(cages);
  seams_system_destroy(sys);
  return hex == 128 && cubic == 0 && water == 0 && both == 0 ? 0 : 1;
}
