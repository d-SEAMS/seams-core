/* R-callable wrappers for the science tutorials. Compiled in the build tree. */
#include "seams.h"

#include <stdlib.h>
#include <string.h>

void seams_r_n(char **path, int *n_out) {
  seams_system *sys = seams_system_create();
  if (seams_system_read(sys, path[0], 1, -1) != SEAMS_OK) {
    *n_out = -1;
    seams_system_destroy(sys);
    return;
  }
  *n_out = seams_system_count(sys);
  seams_system_destroy(sys);
}

void seams_r_frame(char **path, int *n_out, int *labels, int *counts) {
  seams_system *sys = seams_system_create();
  if (seams_system_read(sys, path[0], 1, -1) != SEAMS_OK) {
    *n_out = -1;
    seams_system_destroy(sys);
    return;
  }
  int n = seams_system_count(sys);
  if (n != *n_out) {
    *n_out = -1;
    seams_system_destroy(sys);
    return;
  }
  int local[9] = {0};
  if (seams_chill_plus(sys, 3.5, 0, 0, 0, labels, local) != SEAMS_OK) {
    *n_out = -1;
  } else {
    for (int i = 0; i < 9; i++) {
      counts[i] = local[i];
    }
  }
  seams_system_destroy(sys);
}

void seams_r_cages(char **path, char **spec, int *n_cages) {
  seams_system *sys = seams_system_create();
  size_t cages = 0;
  size_t atoms = 0;
  if (seams_system_read(sys, path[0], 1, -1) != SEAMS_OK ||
      seams_cages_signature(sys, 3.5, 0, 0, spec[0], &cages, &atoms) != SEAMS_OK) {
    *n_cages = -1;
  } else {
    *n_cages = (int)cages;
  }
  seams_system_destroy(sys);
}
