/* Neighbours, rings, and CHILL+ through the C API.
 * The caller owns every buffer and every handle, and every call
 * returns a status code. */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int die(seams_system *sys, seams_neighbours *nb, seams_rings *rings) {
  const char *err = seams_last_error();
  fprintf(stderr, "%s\n", err != NULL ? err : "seams call failed");
  seams_rings_destroy(rings);
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
  return 1;
}

int main(void) {
  const char *root = getenv("DSEAMS_ROOT");
  char path[512];
  seams_system *sys;
  seams_neighbours *nb = NULL;
  seams_rings *rings = NULL;
  seams_status status;
  int *labels = NULL;
  int counts[9];
  int histogram[7];
  size_t nnz = 0;
  size_t n_offsets = 0;
  int n;
  int i;

  if (root == NULL) {
    fprintf(stderr, "DSEAMS_ROOT is unset\n");
    return 1;
  }
  snprintf(path, sizeof path, "%s/input/traj/ice1h.lammpstrj", root);
  sys = seams_system_create();
  if (sys == NULL) {
    return die(NULL, NULL, NULL);
  }
  status = seams_system_read(sys, path, 1, -1);
  if (status != SEAMS_OK) {
    return die(sys, NULL, NULL);
  }
  n = seams_system_count(sys);
  status = seams_neighbours_build(sys, 3.5, -1, SEAMS_GRAPH_CUTOFF, &nb);
  if (status != SEAMS_OK) {
    return die(sys, nb, NULL);
  }
  status = seams_neighbours_copy(nb, NULL, 0, &n_offsets, NULL, 0, &nnz);
  if (status != SEAMS_OK) {
    return die(sys, nb, NULL);
  }
  status = seams_rings_find(sys, 3.5, -1, SEAMS_GRAPH_CUTOFF, 6, &rings);
  if (status != SEAMS_OK) {
    return die(sys, nb, rings);
  }
  memset(histogram, 0, sizeof histogram);
  status = seams_rings_histogram(rings, histogram, 7);
  if (status != SEAMS_OK) {
    return die(sys, nb, rings);
  }
  labels = calloc((size_t)n, sizeof *labels);
  if (labels == NULL) {
    return die(sys, nb, rings);
  }
  memset(counts, 0, sizeof counts);
  status = seams_chill_plus(sys, 3.5, -1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                             labels, counts);
  if (status != SEAMS_OK) {
    free(labels);
    return die(sys, nb, rings);
  }
  printf("pairs %zu rings", nnz);
  for (i = 0; i < 7; i++) {
    printf(" %d", histogram[i]);
  }
  printf(" chill");
  for (i = 0; i < 9; i++) {
    printf(" %d", counts[i]);
  }
  printf("\n");
  free(labels);
  seams_rings_destroy(rings);
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
  return 0;
}
