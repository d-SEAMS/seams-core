/* Neighbours: build the 3.5 A list of hexagonal ice and copy it out as CSR.
 * Usage: neighbours <repository root> */
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
  seams_system *sys = seams_system_create();
  seams_neighbours *nb = NULL;
  int *offsets = NULL;
  int *cols = NULL;
  if (sys == NULL || seams_system_read(sys, path, 1, 1) != SEAMS_OK) {
    fprintf(stderr, "read: %s\n", seams_last_error());
    goto done;
  }
  if (seams_neighbours_build(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_OK) {
    fprintf(stderr, "neighbours: %s\n", seams_handle_error(sys));
    goto done;
  }
  /* Query the lengths, then copy into caller buffers. */
  size_t n_offsets = 0;
  size_t nnz = 0;
  if (seams_neighbours_copy(nb, NULL, 0, &n_offsets, NULL, 0, &nnz) != SEAMS_OK) {
    goto done;
  }
  offsets = malloc(n_offsets * sizeof *offsets);
  cols = malloc((nnz > 0 ? nnz : 1) * sizeof *cols);
  if (offsets == NULL || cols == NULL ||
      seams_neighbours_copy(nb, offsets, n_offsets, &n_offsets, cols, nnz, &nnz) != SEAMS_OK) {
    goto done;
  }
  /* Every oxygen of hexagonal ice has four neighbours. */
  for (size_t i = 0; i + 1 < n_offsets; i++) {
    if (offsets[i + 1] - offsets[i] != 4) {
      fprintf(stderr, "row %zu has %d neighbours\n", i, offsets[i + 1] - offsets[i]);
      goto done;
    }
  }
  printf("neighbours rows %zu pairs %zu\n", n_offsets - 1, nnz);
  rc = 0;

done:
  free(cols);
  free(offsets);
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
  return rc;
}
