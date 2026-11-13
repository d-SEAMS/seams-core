/* libFuzzer harness for seams_system_set_atoms / set_frame / set_ortho.
 * A rejected frame must return a status and a message, never abort.
 */
#include "seams.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int explain(seams_status status, seams_system *sys) {
  if (status == SEAMS_OK) {
    return 1;
  }
  const char *message = seams_last_error();
  if (message == NULL || message[0] == '\0') {
    abort();
  }
  (void)seams_handle_error(sys);
  return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    return 0;
  }
  int n = 0;
  if (size >= 4) {
    memcpy(&n, data, sizeof n);
    n = (n & 31) + 1;
  }
  size_t need = (size_t)n * 3 * sizeof(double) + 9 * sizeof(double);
  if (n > 0 && size >= need) {
    double xyz[32 * 3];
    double box[9];
    int types[32];
    int ids[32];
    memcpy(xyz, data, (size_t)n * 3 * sizeof(double));
    memcpy(box, data + (size_t)n * 3 * sizeof(double), sizeof box);
    for (int i = 0; i < n; ++i) {
      types[i] = (int)(data[i % size] % 4);
      ids[i] = i + 1;
    }
    seams_status status = seams_system_set_atoms(sys, xyz, (size_t)n, box, 9, types, ids);
    if (explain(status, sys)) {
      double q = 0;
      (void)seams_steinhardt(sys, 3.5, 0, 6, &q, 1, NULL, 0);
      (void)seams_system_count(sys);
    }
    status = seams_system_set_frame(sys, xyz, n, box, types, ids);
    (void)explain(status, sys);
    status = seams_system_set_ortho(sys, xyz, n, box[0], box[1], box[2], box[3], box[4],
                                    box[5], types, ids);
    (void)explain(status, sys);
  }
  if (seams_system_set_atoms(sys, NULL, 0, NULL, 0, NULL, NULL) != SEAMS_ERR_NULL) {
    abort();
  }
  seams_system_destroy(sys);
  seams_system_destroy(sys);
  return 0;
}
