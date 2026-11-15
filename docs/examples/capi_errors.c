/* A missing path stores a thread-local error and a non-OK status. */
#include "seams.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  seams_system *sys = seams_system_create();
  const char *null_message = NULL;
  seams_status missing =
      seams_system_natoms(NULL, NULL);
  null_message = seams_last_error_message();
  if (missing != SEAMS_ERR_NULL || null_message == NULL || null_message[0] == '\0' ||
      strcmp(null_message, seams_last_error()) != 0) {
    fprintf(stderr, "null status %d message %s\n", (int)missing,
            null_message ? null_message : "(null)");
    return 1;
  }
  seams_status status =
      seams_system_read(sys, "input/traj/no-such-frame.lammpstrj", 1, -1);
  const char *err = seams_last_error();
  const char *message = seams_last_error_message();
  const char *name = seams_status_string(status);
  printf("status %d\n", (int)status);
  if (status == SEAMS_OK || err == NULL || err[0] == '\0' || message == NULL ||
      strcmp(err, message) != 0 || name == NULL || name[0] == '\0') {
    fprintf(stderr, "status %d err %s\n", (int)status, err ? err : "(null)");
    return 1;
  }
  seams_system_destroy(sys);
  return 0;
}
