/* Ice Ih has 256 primitive six-rings at 3.5 angstroms. */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_rings *rings = NULL;
  int hist[13] = {0};
  int sized[13] = {0};
  size_t nrings = 0;
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK ||
      seams_rings_histogram(rings, hist, 13) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  if (seams_rings_nrings(rings, &nrings) != SEAMS_OK || nrings != 256 ||
      seams_rings_histogram_n(rings, sized, 13) != SEAMS_OK || sized[6] != 256) {
    fprintf(stderr, "nrings %zu hist6 %d\n", nrings, sized[6]);
    return 1;
  }
  int *sizes = calloc(nrings, sizeof(int));
  int *symmetry = calloc(nrings, sizeof(int));
  if (sizes == NULL || symmetry == NULL) {
    return 1;
  }
  if (seams_rings_labels_n(rings, sizes, nrings, symmetry, nrings) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error_message());
    return 1;
  }
  int basal = 0;
  int prismatic = 0;
  int hexagon = 0;
  for (size_t i = 0; i < nrings; ++i) {
    if (sizes[i] != 6) {
      fprintf(stderr, "ring %zu size %d\n", i, sizes[i]);
      return 1;
    }
    if (symmetry[i] == SEAMS_RING_HC_BASAL) {
      ++basal;
    } else if (symmetry[i] == SEAMS_RING_HC_PRISMATIC) {
      ++prismatic;
    } else if (symmetry[i] == SEAMS_RING_HEXAGON) {
      ++hexagon;
    } else {
      fprintf(stderr, "ring %zu symmetry %d\n", i, symmetry[i]);
      return 1;
    }
  }
  if (basal != 64 || prismatic != 141 || hexagon != 51) {
    fprintf(stderr, "symmetry %d %d %d\n", basal, prismatic, hexagon);
    return 1;
  }
  printf("rings6 %d\n", hist[6]);
  free(sizes);
  free(symmetry);
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return hist[6] == 256 ? 0 : 1;
}
