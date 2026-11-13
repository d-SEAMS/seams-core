/* libFuzzer harness for seams_system_read and seams_path_nframes.
 * Each input is written once per trajectory extension the readers accept.
 */
#define _POSIX_C_SOURCE 200809L

#include "seams.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static const char *k_ext[] = {"xyz", "lammpstrj", "dump", "gro", "pdb", "dcd",
                              "xtc",  "nc",        "cif",  "mmcif", "mol2", "sdf",
                              "smi",  "con",       "data", "lmp"};

static void write_all(const char *bytes, size_t n, char paths[][64]) {
  for (size_t i = 0; i < sizeof k_ext / sizeof k_ext[0]; ++i) {
    snprintf(paths[i], 64, "capi-fuzz-%d.%s", (int)i, k_ext[i]);
    FILE *out = fopen(paths[i], "wb");
    if (out == NULL) {
      paths[i][0] = '\0';
      continue;
    }
    if (n > 0) {
      fwrite(bytes, 1, n, out);
    }
    fclose(out);
  }
}

static void remove_all(char paths[][64]) {
  for (size_t i = 0; i < sizeof k_ext / sizeof k_ext[0]; ++i) {
    if (paths[i][0] != '\0') {
      remove(paths[i]);
    }
  }
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  char paths[16][64];
  memset(paths, 0, sizeof paths);
  write_all((const char *)data, size, paths);
  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    remove_all(paths);
    return 0;
  }
  for (size_t i = 0; i < sizeof k_ext / sizeof k_ext[0]; ++i) {
    if (paths[i][0] == '\0') {
      continue;
    }
    size_t nframes = 0;
    (void)seams_path_nframes(paths[i], &nframes);
    int frame = 1;
    if (size > 0) {
      frame = 1 + (int)(data[0] % 3);
    }
    int filter = 0;
    if (size > 1) {
      filter = (int)data[1] % 5 - 1;
    }
    if (seams_system_read(sys, paths[i], frame, filter) == SEAMS_OK) {
      (void)seams_system_count(sys);
      (void)seams_handle_error(sys);
    } else if (seams_last_error()[0] == '\0') {
      /* A failed read must explain itself. Abort so the fuzzer keeps it. */
      abort();
    }
  }
  if (seams_system_read(sys, NULL, 1, 0) != SEAMS_ERR_NULL) {
    abort();
  }
  seams_system_destroy(sys);
  seams_system_destroy(sys);
  seams_system_destroy(NULL);
  remove_all(paths);
  return 0;
}
