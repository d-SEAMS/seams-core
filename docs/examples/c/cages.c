/* Cages: the hexagonal and double-diamond cage census of hexagonal ice.
 * Usage: cages <repository root> */
#include "seams.h"

#include <stdio.h>

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s <repository root>\n", argv[0]);
    return 2;
  }
  char path[4096];
  snprintf(path, sizeof path, "%s/input/traj/ice1h.lammpstrj", argv[1]);

  int rc = 1;
  seams_system *sys = seams_system_create();
  seams_cage_set *cages = NULL;
  if (sys == NULL || seams_system_read(sys, path, 1, 1) != SEAMS_OK) {
    fprintf(stderr, "read: %s\n", seams_last_error());
    goto done;
  }
  if (seams_cages(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, 4, 0, &cages) != SEAMS_OK) {
    fprintf(stderr, "cages: %s\n", seams_handle_error(sys));
    goto done;
  }
  int n_hex = 0;
  int n_cubic = 0;
  int n_water = 0;
  int n_both = 0;
  if (seams_cages_counts(cages, &n_hex, &n_cubic, &n_water, &n_both) != SEAMS_OK) {
    fprintf(stderr, "counts: %s\n", seams_handle_error(cages));
    goto done;
  }
  /* The four counts sum to the number of atoms the search kept. */
  if (n_hex <= 0 || n_cubic != 0) {
    fprintf(stderr, "hexagonal %d cubic %d\n", n_hex, n_cubic);
    goto done;
  }
  printf("cages hexagonal %d cubic %d both %d other %d\n", n_hex, n_cubic, n_both, n_water);
  rc = 0;

done:
  seams_cages_destroy(cages);
  seams_system_destroy(sys);
  return rc;
}
