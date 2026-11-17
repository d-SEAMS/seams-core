/* Copyright (c) d-SEAMS developers
 * Resident memory of a trajectory loop through the C API.
 *
 * For each ice Ic supercell the program loads the frame with
 * seams_system_set_atoms, labels it with seams_chill_plus_n, and reads
 * VmRSS from /proc/self/status after each step:
 *   loaded      set_atoms and one CHILL+ call
 *   reloads     four seams_system_set_positions calls, default mode
 *   staged      four more after seams_system_stage_positions(sys, 1)
 *   unstaged    after seams_system_stage_positions(sys, 0)
 * Each line gives the change from the step before, as kB, bytes per
 * molecule and percent of the loaded RSS. Linux only.
 *
 * Usage: bench_capi_rss [cells...]   default cells 24 48
 * Cells per edge k give 8 k^3 molecules: 24 -> 110592, 48 -> 884736.
 */
#include <seams.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static long rss_kb(void) {
  FILE *status = fopen("/proc/self/status", "r");
  if (status == NULL) {
    return -1;
  }
  char line[256];
  long kb = -1;
  while (fgets(line, sizeof line, status) != NULL) {
    if (strncmp(line, "VmRSS:", 6) == 0) {
      kb = atol(line + 6);
    }
  }
  fclose(status);
  return kb;
}

static void need(seams_status status, const char *what) {
  if (status != SEAMS_OK) {
    fprintf(stderr, "%s: %s (%s)\n", what, seams_status_string(status),
            seams_last_error_message());
    exit(1);
  }
}

static void step(const char *what, long before, long after, long loaded, size_t n) {
  const long delta = after - before;
  printf("  %-9s rss %8ld kB  %+8ld kB  %+6.1f B/molecule  %+5.1f%% of loaded\n", what, after,
         delta, 1024.0 * (double)delta / (double)n, 100.0 * (double)delta / (double)loaded);
}

static void run(int k) {
  const size_t n = (size_t)8 * k * k * k;
  const double a = 4.0 * 2.75 / sqrt(3.0);
  const double basis[8][3] = {{0, 0, 0},          {0.5, 0.5, 0},      {0.5, 0, 0.5},
                              {0, 0.5, 0.5},      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  double *xyz = malloc(3 * n * sizeof *xyz);
  int *labels = calloc(n, sizeof *labels);
  if (xyz == NULL || labels == NULL) {
    exit(1);
  }
  size_t at = 0;
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) {
      for (int l = 0; l < k; ++l) {
        for (int b = 0; b < 8; ++b, ++at) {
          const double cell[3] = {i, j, l};
          for (int d = 0; d < 3; ++d) {
            xyz[3 * at + d] = (cell[d] + basis[b][d]) * a;
          }
        }
      }
    }
  }
  double box[9] = {0};
  box[3] = box[4] = box[5] = k * a;

  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    exit(1);
  }
  need(seams_system_set_atoms(sys, xyz, n, box, 9, NULL, NULL), "set_atoms");
  int counts[9];
  need(seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n,
                          counts),
       "chill_plus_n");
  const long loaded = rss_kb();
  for (int r = 0; r < 4; ++r) {
    need(seams_system_set_positions(sys, xyz, n, NULL, 0), "set_positions");
  }
  const long reloads = rss_kb();
  need(seams_system_stage_positions(sys, 1), "stage_positions");
  for (int r = 0; r < 4; ++r) {
    need(seams_system_set_positions(sys, xyz, n, NULL, 0), "set_positions");
  }
  const long staged = rss_kb();
  need(seams_system_stage_positions(sys, 0), "stage_positions");
  const long unstaged = rss_kb();

  printf("molecules %zu\n", n);
  printf("  %-9s rss %8ld kB\n", "loaded", loaded);
  step("reloads", loaded, reloads, loaded, n);
  step("staged", reloads, staged, loaded, n);
  step("unstaged", staged, unstaged, loaded, n);
  seams_system_destroy(sys);
  free(labels);
  free(xyz);
}

int main(int argc, char **argv) {
  need(seams_set_num_threads(1), "threads");
  printf("seams %s\n", seams_version_string());
  if (argc < 2) {
    run(24);
    run(48);
  }
  for (int i = 1; i < argc; ++i) {
    run(atoi(argv[i]));
  }
  return 0;
}
