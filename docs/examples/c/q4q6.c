/* q4 and q6: per-atom Steinhardt parameters, then the averaged q6 with its
 * mean, median and ice count.
 * Usage: q4q6 <repository root> */
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
  double *q4 = NULL;
  double *q6 = NULL;
  seams_system *sys = seams_system_create();
  size_t n = 0;
  if (sys == NULL || seams_system_read(sys, path, 1, 1) != SEAMS_OK ||
      seams_system_natoms(sys, &n) != SEAMS_OK) {
    fprintf(stderr, "read: %s\n", seams_last_error());
    goto done;
  }
  q4 = malloc(n * sizeof *q4);
  q6 = malloc(n * sizeof *q6);
  if (q4 == NULL || q6 == NULL || seams_q4q6_n(sys, 3.5, 1, q4, n, q6, n) != SEAMS_OK) {
    fprintf(stderr, "q4q6: %s\n", seams_handle_error(sys));
    goto done;
  }
  double mean = 0.0;
  double median = 0.0;
  int n_ice = 0;
  /* NULL skips the per-atom output; the summary still comes back. */
  if (seams_q6bar(sys, 3.5, 1, 0.36, NULL, 0, &mean, &median, &n_ice) != SEAMS_OK) {
    fprintf(stderr, "q6bar: %s\n", seams_handle_error(sys));
    goto done;
  }
  if (!(q6[0] > q4[0]) || (size_t)n_ice != n) {
    fprintf(stderr, "q4 %g q6 %g ice %d of %zu\n", q4[0], q6[0], n_ice, n);
    goto done;
  }
  printf("q4 %.3f q6 %.3f q6bar mean %.3f ice %d\n", q4[0], q6[0], mean, n_ice);
  rc = 0;

done:
  free(q6);
  free(q4);
  seams_system_destroy(sys);
  return rc;
}
