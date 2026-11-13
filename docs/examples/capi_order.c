/* Published CHILL+ marks every Ice Ih oxygen hexagonal. */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>

int main(void) {
  seams_system *sys = seams_system_create();
  int counts[9] = {0};
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  int n = seams_system_count(sys);
  int *labels = calloc((size_t)n, sizeof(int));
  if (labels == NULL) {
    return 1;
  }
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                       counts) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  printf("hexagonal %d\n", counts[1]);
  free(labels);
  seams_system_destroy(sys);
  return counts[1] == 128 ? 0 : 1;
}
