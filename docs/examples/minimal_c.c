/* Read input/traj/genice_sI.lammpstrj, then neighbours, q6, CHILL+, rings. */
#include "seams.h"
#include "seams_lammps_box.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { kMax = 64 };

static int read_dump(const char *path, double *xyz, int *n, double *lx, double *ly,
                     double *lz) {
  FILE *in = fopen(path, "r");
  char line[256];
  int atoms = 0;
  int reading = 0;
  double lo, hi;
  if (in == NULL) {
    return 1;
  }
  *n = 0;
  while (fgets(line, sizeof line, in) != NULL) {
    if (strncmp(line, "ITEM: BOX BOUNDS", 16) == 0) {
      double bound_lo[3], bound_hi[3], tilt[3], box9[9];
      int axis;
      for (axis = 0; axis < 3; ++axis) {
        char bound[128];
        if (fgets(bound, sizeof bound, in) == NULL) {
          fclose(in);
          return 1;
        }
        tilt[axis] = 0.0;
        if (sscanf(bound, "%lf %lf %lf", &bound_lo[axis], &bound_hi[axis], &tilt[axis]) < 2) {
          fclose(in);
          return 1;
        }
      }
      seams_lammps_bounds(bound_lo, bound_hi, tilt, box9);
      *lx = box9[3];
      *ly = box9[4];
      *lz = box9[5];
      (void)lo;
      (void)hi;
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
  double q4[kMax], q6[kMax], mean, fraction;
  int labels[kMax], counts[9], hist[13];
  int n = 0, nnz = 0, n_offsets = 0, i;
  double lx, ly, lz;
  const int *offsets = NULL;
  const int *cols = NULL;
  seams_system *sys;
  seams_neighbours *nb = NULL;
  seams_rings *rings = NULL;
  const char *path = argc > 1 ? argv[1] : "input/traj/genice_sI.lammpstrj";
  if (read_dump(path, xyz, &n, &lx, &ly, &lz) != 0 || n != 46) {
    fprintf(stderr, "read failed\n");
    return 1;
  }
  if (seams_set_num_threads(1) != SEAMS_OK) {
    return 1;
  }
  sys = seams_system_create();
  if (seams_system_set_ortho(sys, xyz, n, 0, 0, 0, lx, ly, lz, NULL, NULL) != SEAMS_OK) {
    return 1;
  }
  if (seams_system_count(sys) != 46) {
    return 1;
  }
  if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_OK) {
    return 1;
  }
  if (seams_neighbours_count(nb) != 46) {
    return 1;
  }
  if (seams_neighbours_csr(nb, &offsets, &n_offsets, &cols, &nnz) != SEAMS_OK || nnz < 150) {
    return 1;
  }
  seams_neighbours_destroy(nb);
  if (seams_q4q6(sys, 3.5, 0, q4, q6) != SEAMS_OK) {
    return 1;
  }
  mean = 0;
  for (i = 0; i < n; i++) {
    mean += q6[i];
  }
  mean /= n;
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                       counts) != SEAMS_OK) {
    return 1;
  }
  if (counts[SEAMS_LABEL_CLATHRATE] != 46) {
    return 1;
  }
  if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    return 1;
  }
  if (seams_rings_histogram(rings, hist, 13) != SEAMS_OK) {
    return 1;
  }
  if (seams_rings_fraction(rings, 6, &fraction) != SEAMS_OK) {
    return 1;
  }
  if (hist[5] != 48 || hist[6] != 54 || seams_rings_count(rings) != 102) {
    return 1;
  }
  printf("%s atoms %d clathrate %d q6 %.3f rings %d\n", seams_version_string(), n,
         counts[SEAMS_LABEL_CLATHRATE], mean, seams_rings_count(rings));
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  (void)fraction;
  return mean > 0.5 ? 0 : 1;
}
