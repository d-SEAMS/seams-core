/* Errors: status codes, the thread-local message and the copy kept on the
 * handle that failed.
 * Usage: errors */
#include "seams.h"

#include <stdio.h>
#include <string.h>

int main(void) {
  int rc = 1;
  seams_system *sys = seams_system_create();
  seams_neighbours *nb = NULL;
  if (sys == NULL) {
    return 1;
  }
  /* A system with no frame cannot build neighbours. */
  seams_status status = seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb);
  if (status != SEAMS_ERR_STATE || nb != NULL) {
    fprintf(stderr, "expected SEAMS_ERR_STATE, got %d\n", (int)status);
    goto done;
  }
  /* The code's fixed text, the detail for this thread, and the handle's copy. */
  const char *text = seams_status_string(status);
  const char *thread_detail = seams_last_error();
  const char *handle_detail = seams_handle_error(sys);
  if (text[0] == '\0' || thread_detail[0] == '\0' || handle_detail[0] == '\0' ||
      strcmp(seams_status_string(SEAMS_OK), text) == 0) {
    fprintf(stderr, "messages: '%s' '%s' '%s'\n", text, thread_detail, handle_detail);
    goto done;
  }
  printf("error %s: %s\n", text, handle_detail);
  rc = 0;

done:
  seams_system_destroy(sys);
  return rc;
}
