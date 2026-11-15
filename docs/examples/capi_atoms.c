/* Read Ice Ih and require 128 atoms. */
#include "seams.h"

#include <stdio.h>

int main(void) {
  seams_system *sys = seams_system_create();
  size_t natoms = 0;
  if (sys == NULL) {
    return 1;
  }
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  int n = seams_system_count(sys);
  if (seams_system_natoms(sys, &natoms) != SEAMS_OK || natoms != 128 || n != 128) {
    fprintf(stderr, "natoms %zu count %d\n", natoms, n);
    return 1;
  }
  printf("atoms %d\n", n);
  seams_system_destroy(sys);
  return 0;
}
