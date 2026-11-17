/* Copyright (c) d-SEAMS developers
 * Installed-header and installed-library versions must agree.
 * Usage: version [expected-version-string]
 */
#include <seams.h>

#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
  if (seams_version() != SEAMS_VERSION_NUMBER || seams_abi() != SEAMS_ABI) {
    fprintf(stderr, "header %d abi %d, library %d abi %d\n", SEAMS_VERSION_NUMBER, SEAMS_ABI,
            seams_version(), seams_abi());
    return 1;
  }
  if (strcmp(seams_version_string(), SEAMS_VERSION_STRING) != 0) {
    fprintf(stderr, "header %s, library %s\n", SEAMS_VERSION_STRING, seams_version_string());
    return 1;
  }
  if (argc > 1 && strcmp(argv[1], SEAMS_VERSION_STRING) != 0) {
    fprintf(stderr, "expected %s, header %s\n", argv[1], SEAMS_VERSION_STRING);
    return 1;
  }
  printf("seams %s abi %d\n", seams_version_string(), seams_abi());
  return 0;
}
