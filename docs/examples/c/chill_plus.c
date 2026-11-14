/* CHILL+: per-atom ice labels and the nine class counts.
 * Usage: chill_plus <repository root> */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s <repository root>\n", argv[0]);
    return 2;
  }
  char path[4096];
  snprintf(path, sizeof path, "%s/input/traj/ice1h.lammpstrj", argv[1]);

  int rc = 1;
  int *labels = NULL;
  seams_system *sys = seams_system_create();
  size_t n = 0;
  if (sys == NULL || seams_system_read(sys, path, 1, 1) != SEAMS_OK ||
      seams_system_natoms(sys, &n) != SEAMS_OK) {
    fprintf(stderr, "read: %s\n", seams_last_error());
    goto done;
  }
  labels = malloc(n * sizeof *labels);
  int counts9[9];
  if (labels == NULL ||
      seams_chill_plus_n(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n,
                         counts9) != SEAMS_OK) {
    fprintf(stderr, "chill+: %s\n", seams_handle_error(sys));
    goto done;
  }
  /* counts9 is indexed by the SEAMS_LABEL_* codes. */
  if ((size_t)counts9[SEAMS_LABEL_HEXAGONAL] != n || labels[0] != SEAMS_LABEL_HEXAGONAL) {
    fprintf(stderr, "%d of %zu atoms are hexagonal\n", counts9[SEAMS_LABEL_HEXAGONAL], n);
    goto done;
  }
  printf("chill+ hexagonal %d of %zu\n", counts9[SEAMS_LABEL_HEXAGONAL], n);
  rc = 0;

done:
  free(labels);
  seams_system_destroy(sys);
  return rc;
}
