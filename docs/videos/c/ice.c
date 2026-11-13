/* Ice Ih census through seams.h. Input default: input/traj/ice1h.lammpstrj. */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { kMax = 256 };

static int read_dump(const char *path, double *xyz, int *n, double *lx, double *ly, double *lz) {
  FILE *in = fopen(path, "r");
  char line[512];
  int atoms = 0;
  int reading = 0;
  double lo, hi;
  if (in == NULL) {
    return 1;
  }
  while (fgets(line, sizeof line, in) != NULL) {
    if (strncmp(line, "ITEM: BOX BOUNDS", 16) == 0) {
      if (fscanf(in, "%lf %lf", &lo, &hi) != 2) {
        fclose(in);
        return 1;
      }
      *lx = hi - lo;
      if (fscanf(in, "%lf %lf", &lo, &hi) != 2) {
        fclose(in);
        return 1;
      }
      *ly = hi - lo;
      if (fscanf(in, "%lf %lf", &lo, &hi) != 2) {
        fclose(in);
        return 1;
      }
      *lz = hi - lo;
    } else if (strncmp(line, "ITEM: ATOMS", 11) == 0) {
      reading = 1;
    } else if (reading) {
      int id, type;
      double x, y, z;
      if (sscanf(line, "%d %d %lf %lf %lf", &id, &type, &x, &y, &z) != 5) {
        break;
      }
      if (atoms >= kMax) {
        fclose(in);
        return 1;
      }
      xyz[3 * atoms] = x;
      xyz[3 * atoms + 1] = y;
      xyz[3 * atoms + 2] = z;
      atoms++;
    }
  }
  fclose(in);
  *n = atoms;
  return atoms > 0 ? 0 : 1;
}

int main(int argc, char **argv) {
  double xyz[3 * kMax];
  double q4[kMax], q6[kMax], mean;
  int labels[kMax], counts[9], hist[13];
  int n = 0, i;
  double lx, ly, lz;
  seams_system *sys;
  seams_neighbours *nb = NULL;
  seams_rings *rings = NULL;
  int nI = 0;
  const char *path = argc > 1 ? argv[1] : "input/traj/ice1h.lammpstrj";
  if (read_dump(path, xyz, &n, &lx, &ly, &lz) != 0 || n != 128) {
    fprintf(stderr, "read failed\n");
    return 1;
  }
  sys = seams_system_create();
  if (seams_system_set_ortho(sys, xyz, n, 0, 0, 0, lx, ly, lz, NULL, NULL) != SEAMS_OK) {
    return 1;
  }
  if (seams_neighbours_build(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_OK) {
    return 1;
  }
  nI = seams_neighbours_count(nb);
  if (nI != 128) {
    return 1;
  }
  seams_neighbours_destroy(nb);
  nb = NULL;
  if (seams_q4q6(sys, 3.5, 1, q4, q6) != SEAMS_OK) {
    return 1;
  }
  mean = 0;
  for (i = 0; i < n; i++) {
    mean += q6[i];
  }
  mean /= n;
  if (seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts) !=
      SEAMS_OK) {
    return 1;
  }
  if (seams_rings_find(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    return 1;
  }
  if (seams_rings_histogram(rings, hist, 13) != SEAMS_OK) {
    return 1;
  }
  printf("nop %d nI %d q6 %.6f hexagonal %d rings %d\n", n, nI, mean,
         counts[SEAMS_LABEL_HEXAGONAL], hist[6]);
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return 0;
}
