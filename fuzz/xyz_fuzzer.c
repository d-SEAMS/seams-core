/* libFuzzer harness for the xyz reader behind seams_system_read. */
#define _POSIX_C_SOURCE 200809L

#include "seams.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  char path[] = "fuzz-xyz-XXXXXX";
  int fd = mkstemp(path);
  if (fd < 0) {
    return 0;
  }
  /* Keep a .xyz suffix so the reader selects the xyz path. */
  char named[80];
  snprintf(named, sizeof named, "%s.xyz", path);
  if (rename(path, named) != 0) {
    close(fd);
    unlink(path);
    return 0;
  }
  if (size > 0 && write(fd, data, size) < 0) {
    close(fd);
    unlink(named);
    return 0;
  }
  close(fd);
  seams_system *sys = seams_system_create();
  if (sys != NULL) {
    int frame = 1;
    if (size > 0) {
      frame = 1 + (int)(data[0] % 3);
    }
    size_t nframes = 0;
    (void)seams_path_nframes(named, &nframes);
    if (seams_system_read(sys, named, frame, 0) == SEAMS_OK) {
      (void)seams_system_count(sys);
    } else if (seams_last_error()[0] == '\0') {
      abort();
    }
    seams_system_destroy(sys);
  }
  unlink(named);
  return 0;
}
