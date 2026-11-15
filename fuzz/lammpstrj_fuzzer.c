/* libFuzzer harness for the LAMMPS dump reader behind seams_system_read. */
#define _POSIX_C_SOURCE 200809L

#include "seams.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  char path[] = "fuzz-lammpstrj-XXXXXX";
  int fd = mkstemp(path);
  if (fd < 0) {
    return 0;
  }
  if (size > 0) {
    if (write(fd, data, size) < 0) {
      close(fd);
      unlink(path);
      return 0;
    }
  }
  close(fd);
  seams_system *sys = seams_system_create();
  if (sys != NULL) {
    int frame = 1;
    int filter = 0;
    if (size > 0) {
      frame = 1 + (int)(data[0] % 4);
    }
    if (size > 1) {
      filter = (int)data[1] % 5 - 1;
    }
    size_t nframes = 0;
    (void)seams_path_nframes(path, &nframes);
    if (seams_system_read(sys, path, frame, filter) == SEAMS_OK) {
      (void)seams_system_count(sys);
    } else if (seams_last_error()[0] == '\0') {
      abort();
    }
    seams_system_destroy(sys);
  }
  unlink(path);
  return 0;
}
