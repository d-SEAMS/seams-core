/* Ice Ih census through the language-front ABI.
 *
 * input/traj/ice1h.lammpstrj, every atom: 128 hexagonal, 256 six-rings.
 */
#include "seams_lang.h"

#include <stdio.h>
#include <stdlib.h>

int main(void) {
  seams_frame *frame = seams_frame_read("input/traj/ice1h.lammpstrj", 1, -1);
  if (frame == NULL || seams_frame_nop(frame) != 128) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  int n = seams_frame_nop(frame);
  int *labels = calloc((size_t)n, sizeof(int));
  int counts[9] = {0};
  int rings[13] = {0};
  if (labels == NULL ||
      seams_frame_chill_plus(frame, 3.5, 0, "cutoff", labels, counts) != 0 ||
      seams_frame_rings(frame, 3.5, 0, "cutoff", 6, rings, 13) != 0) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  printf("hexagonal %d\nrings6 %d\n", counts[1], rings[6]);
  free(labels);
  seams_frame_free(frame);
  return counts[1] == 128 && rings[6] == 256 ? 0 : 1;
}
