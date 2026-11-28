<!-- https://docs.dseams.info/fronts/c.md -->

# C tutorial {#c-tutorial}

Label the eight-atom cubic cell with `seams.h` and print the library version. The program is `docs/examples/ice.c`. It prints the project version and `cubic 8`.

```c
// Copyright (c) d-SEAMS developers
#include "seams.h"
#include <math.h>
#include <stdio.h>

int main(void) {
  const double a = 2.75 * 4.0 / sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  double xyz[24];
  int n = 0;
  for (int f = 0; f < 4; f++) {
    for (int s = 0; s < 2; s++, n++) {
      xyz[3 * n] = (fcc[f][0] + basis[s][0]) * a;
      xyz[3 * n + 1] = (fcc[f][1] + basis[s][1]) * a;
      xyz[3 * n + 2] = (fcc[f][2] + basis[s][2]) * a;
    }
  }
  int labels[8];
  int counts[9];
  seams_system *sys = seams_system_create();
  if (seams_system_set_ortho(sys, xyz, 8, 0, 0, 0, a, a, a, NULL, NULL) != SEAMS_OK) {
    return 1;
  }
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                       counts) != SEAMS_OK) {
    return 1;
  }
  printf("%s cubic %d\n", seams_version_string(), counts[0]);
  seams_system_destroy(sys);
  return counts[0] == 8 ? 0 : 1;
}
```

`seams.h` includes the generated `seams_version.h`. Meson compiles the program against `libseams` and records the library path in the binary.

```bash
cc -Iinclude -Ibuild/src docs/examples/ice.c -Lbuild/src -lseams -Wl,-rpath,build/src -lm -o ice
./ice
```

`docs/examples/core_counts.c` reads `input/traj/ice1h.lammpstrj` and `input/traj/genice_sII.lammpstrj`. It returns 0 when Ice Ih is 128 hexagonal waters and 256 six-rings, and sII is 136 clathrate waters, 144 five-rings and 16 six-rings. The lines it prints are:

```default
hexagonal 128
rings6 256
clathrate 136
rings5 144
```

```c
// Copyright (c) d-SEAMS developers
/* Ice Ih and sII counts from the C core.
 *
 * input/traj/ice1h.lammpstrj: 128 waters, all hexagonal, 256 six-rings.
 * input/traj/genice_sII.lammpstrj: 136 waters, all clathrate, 144 five-rings
 * and 16 six-rings. seams_chill_plus and seams_rings_histogram at 3.5 Å.
 */
#include "seams.h"
#include "seams_lammps_box.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int read_dump(const char *path, double **xyz, int **types, int *n, double box9[9]) {
  FILE *in = fopen(path, "r");
  char line[512];
  int reading = 0;
  int cap = 0;
  if (in == NULL) {
    return 1;
  }
  *n = 0;
  *xyz = NULL;
  *types = NULL;
  while (fgets(line, sizeof line, in) != NULL) {
    if (strncmp(line, "ITEM: BOX BOUNDS", 16) == 0) {
      double bound_lo[3], bound_hi[3], tilt[3];
      int axis;
      for (axis = 0; axis < 3; ++axis) {
        char bound[128];
        if (fgets(bound, sizeof bound, in) == NULL) {
          fclose(in);
          return 1;
        }
        tilt[axis] = 0.0;
        if (sscanf(bound, "%lf %lf %lf", &bound_lo[axis], &bound_hi[axis], &tilt[axis]) < 2) {
          fclose(in);
          return 1;
        }
      }
      seams_lammps_bounds(bound_lo, bound_hi, tilt, box9);
    } else if (strncmp(line, "ITEM: ATOMS", 11) == 0) {
      reading = 1;
    } else if (reading) {
      int id, kind;
      double x, y, z;
      if (sscanf(line, "%d %d %lf %lf %lf", &id, &kind, &x, &y, &z) != 5) {
        break;
      }
      if (*n >= cap) {
        cap = cap == 0 ? 64 : cap * 2;
        *xyz = realloc(*xyz, (size_t)cap * 3 * sizeof(double));
        *types = realloc(*types, (size_t)cap * sizeof(int));
        if (*xyz == NULL || *types == NULL) {
          fclose(in);
          return 1;
        }
      }
      (*types)[*n] = kind;
      (*xyz)[3 * *n] = x;
      (*xyz)[3 * *n + 1] = y;
      (*xyz)[3 * *n + 2] = z;
      *n += 1;
    }
  }
  fclose(in);
  return *n == 0;
}

static int classify(const char *path, int *n_out, int counts[9], int rings[13]) {
  double *xyz = NULL;
  int *types = NULL;
  int n = 0;
  double box9[9];
  seams_system *sys;
  seams_rings *found = NULL;
  int *labels;
  if (read_dump(path, &xyz, &types, &n, box9) != 0) {
    return 1;
  }
  labels = calloc((size_t)n, sizeof(int));
  sys = seams_system_create();
  if (labels == NULL || sys == NULL ||
      seams_system_set_frame(sys, xyz, n, box9, types, NULL) != SEAMS_OK ||
      seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                       counts) != SEAMS_OK ||
      seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &found) != SEAMS_OK ||
      seams_rings_histogram(found, rings, 13) != SEAMS_OK) {
    free(labels);
    free(xyz);
    free(types);
    if (found != NULL) {
      seams_rings_destroy(found);
    }
    if (sys != NULL) {
      seams_system_destroy(sys);
    }
    return 1;
  }
  *n_out = n;
  seams_rings_destroy(found);
  seams_system_destroy(sys);
  free(labels);
  free(xyz);
  free(types);
  return 0;
}

int main(void) {
  int n = 0;
  int counts[9] = {0};
  int rings[13] = {0};
  if (seams_set_num_threads(1) != SEAMS_OK) {
    return 1;
  }
  if (classify("input/traj/ice1h.lammpstrj", &n, counts, rings) != 0 || n != 128 ||
      counts[1] != 128 || rings[6] != 256) {
    fprintf(stderr, "ice Ih n=%d hexagonal=%d rings6=%d\n", n, counts[1], rings[6]);
    return 1;
  }
  printf("hexagonal 128\n");
  printf("rings6 256\n");
  if (classify("input/traj/genice_sII.lammpstrj", &n, counts, rings) != 0 || n != 136 ||
      counts[4] != 136 || rings[5] != 144 || rings[6] != 16) {
    fprintf(stderr, "sII n=%d clathrate=%d rings5=%d rings6=%d\n", n, counts[4], rings[5],
            rings[6]);
    return 1;
  }
  printf("clathrate 136\n");
  printf("rings5 144\n");
  return 0;
}
```

The function list, status codes and buffer lengths are the [C API](../stable/c.md).

## Ice Ih {#ice-ih}

![Ice Ih colored by published CHILL+](../_images/overlay-web14.png)

*Ice Ih colored by published CHILL+.*
