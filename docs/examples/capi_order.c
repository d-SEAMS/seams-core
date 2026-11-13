/* Published CHILL+ marks every Ice Ih oxygen hexagonal. */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>

int main(void) {
  seams_system *sys = seams_system_create();
  int counts[9] = {0};
  int sized[9] = {0};
  size_t natoms = 0;
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_system_natoms(sys, &natoms) != SEAMS_OK || natoms != 128) {
    return 1;
  }
  int n = (int)natoms;
  int *labels = calloc(natoms, sizeof(int));
  double *q4 = calloc(natoms, sizeof(double));
  double *q6 = calloc(natoms, sizeof(double));
  if (labels == NULL || q4 == NULL || q6 == NULL) {
    return 1;
  }
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                       counts) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  if (seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                         natoms, sized) != SEAMS_OK ||
      sized[SEAMS_LABEL_HEXAGONAL] != 128) {
    fprintf(stderr, "chill_plus_n hexagonal %d %s\n", sized[SEAMS_LABEL_HEXAGONAL],
            seams_last_error_message());
    return 1;
  }
  if (seams_q4q6_n(sys, 3.5, 0, q4, natoms, q6, natoms) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error_message());
    return 1;
  }
  for (int i = 0; i < n; ++i) {
    if (!(q4[i] >= 0.0 && q4[i] <= 1.0 && q6[i] > 0.2 && q6[i] <= 1.0)) {
      fprintf(stderr, "q4q6 %d %g %g\n", i, q4[i], q6[i]);
      return 1;
    }
  }
  printf("hexagonal %d\n", counts[1]);
  free(labels);
  free(q4);
  free(q6);
  seams_system_destroy(sys);
  return counts[1] == 128 ? 0 : 1;
}
