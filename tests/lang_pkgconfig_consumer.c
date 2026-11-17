/* Copyright (c) d-SEAMS developers */
/* A downstream C program that knows seams-lang only through pkg-config.
 * tests/check_lang_pkgconfig.py builds it against meson-uninstalled/ and
 * against an install. Counting an .xyz goes through the vendored
 * chemfiles' seams_with_frame when the build has chemfiles. */
#include <seams_lang.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: %s VERSION FILE.xyz\n", argv[0]);
    return 2;
  }
  const char *version = seams_lang_version();
  if (version == NULL || strcmp(version, argv[1]) != 0) {
    fprintf(stderr, "seams_lang_version is %s, expected %s\n", version ? version : "(null)", argv[1]);
    return 1;
  }
  const int frames = seams_path_nframes(argv[2]);
  if (frames < 1) {
    fprintf(stderr, "seams_path_nframes: %s\n", seams_last_error());
    return 1;
  }
  seams_frame *frame = seams_frame_read(argv[2], 1, 0);
  if (frame == NULL) {
    fprintf(stderr, "seams_frame_read: %s\n", seams_last_error());
    return 1;
  }
  printf("seams-lang %s frames %d atoms %d\n", version, frames, seams_frame_nop(frame));
  seams_frame_free(frame);
  return 0;
}
