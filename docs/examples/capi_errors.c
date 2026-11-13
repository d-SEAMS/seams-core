/* A missing path stores a thread-local error and a non-OK status. */
#include "seams.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_status status =
      seams_system_read(sys, "input/traj/no-such-frame.lammpstrj", 1, -1);
  const char *err = seams_last_error();
  const char *name = seams_status_string(status);
  printf("status %d\n", (int)status);
  if (status == SEAMS_OK || err == NULL || err[0] == '\0' || name == NULL || name[0] == '\0') {
    fprintf(stderr, "status %d err %s\n", (int)status, err ? err : "(null)");
    return 1;
  }
  seams_system_destroy(sys);
  return 0;
}
