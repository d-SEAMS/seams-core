/* Rings: primitive rings of hexagonal ice, their size histogram and the
 * six-ring fraction.
 * Usage: rings <repository root> */
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
  seams_rings *rings = NULL;
  if (sys == NULL || seams_system_read(sys, path, 1, 1) != SEAMS_OK) {
    fprintf(stderr, "read: %s\n", seams_last_error());
    goto done;
  }
  if (seams_rings_find(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    fprintf(stderr, "rings: %s\n", seams_handle_error(sys));
    goto done;
  }
  size_t n = 0;
  int counts[7] = {0};
  double six = 0.0;
  if (seams_rings_nrings(rings, &n) != SEAMS_OK ||
      seams_rings_histogram_n(rings, counts, 7) != SEAMS_OK ||
      seams_rings_fraction(rings, 6, &six) != SEAMS_OK) {
    fprintf(stderr, "ring queries: %s\n", seams_handle_error(rings));
    goto done;
  }
  /* counts[k] is the number of k-membered rings; hexagonal ice has only sixes. */
  if (n == 0 || (size_t)counts[6] != n || six != 1.0) {
    fprintf(stderr, "rings %zu, six-rings %d, fraction %g\n", n, counts[6], six);
    goto done;
  }
  printf("rings %zu six-ring fraction %.2f\n", n, six);
  rc = 0;

done:
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return rc;
}
