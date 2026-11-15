/* Copyright (c) d-SEAMS developers
 * One process times the C API. The parent sends "threads op path".
 * k = 4 mutual neighbours, cutoff 3.5 A, q6 and q6bar, published CHILL
 * on the cutoff graph, primitive rings through size 6.
 */
#define _POSIX_C_SOURCE 200809L

#include "seams.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { kSlots = 8, kLine = 4096 };

typedef struct Slot {
  char *path;
  seams_system *sys;
  int n;
  int *labels;
  double *ql;
  double *bar;
} Slot;

static Slot g_slots[kSlots];
static double g_sink;

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void die(const char *message) {
  fprintf(stderr, "error %s\n", message);
  exit(1);
}

static void retune(int threads) {
  char text[32];
  if (threads < 1) {
    die("thread count");
  }
  snprintf(text, sizeof(text), "%d", threads);
  if (setenv("OMP_NUM_THREADS", text, 1) != 0 || setenv("PLUMED_NUM_THREADS", text, 1) != 0) {
    die("setenv");
  }
  if (seams_set_num_threads(threads) != SEAMS_OK) {
    die(seams_last_error());
  }
}

static Slot *slot_for(const char *path) {
  int i;
  for (i = 0; i < kSlots; i++) {
    if (g_slots[i].path != NULL && strcmp(g_slots[i].path, path) == 0) {
      return &g_slots[i];
    }
  }
  for (i = 0; i < kSlots; i++) {
    if (g_slots[i].path == NULL) {
      g_slots[i].path = strdup(path);
      g_slots[i].sys = seams_system_create();
      if (g_slots[i].path == NULL || g_slots[i].sys == NULL) {
        die("system");
      }
      return &g_slots[i];
    }
  }
  die("too many frames");
  return NULL;
}

static void ensure_frame(Slot *slot) {
  if (slot->n > 0) {
    return;
  }
  if (seams_system_read(slot->sys, slot->path, 1, -1) != SEAMS_OK) {
    die(seams_handle_error(slot->sys));
  }
  slot->n = seams_system_count(slot->sys);
  if (slot->n <= 0) {
    die("empty frame");
  }
}

static void buffers(Slot *slot) {
  ensure_frame(slot);
  if (slot->labels != NULL) {
    return;
  }
  slot->labels = (int *)malloc((size_t)slot->n * sizeof(int));
  slot->ql = (double *)malloc((size_t)slot->n * sizeof(double));
  slot->bar = (double *)malloc((size_t)slot->n * sizeof(double));
  if (slot->labels == NULL || slot->ql == NULL || slot->bar == NULL) {
    die("buffer");
  }
}

static double time_read(Slot *slot) {
  double t0 = now_sec();
  if (seams_system_read(slot->sys, slot->path, 1, -1) != SEAMS_OK) {
    die(seams_handle_error(slot->sys));
  }
  double dt = now_sec() - t0;
  slot->n = seams_system_count(slot->sys);
  g_sink += (double)slot->n;
  return dt;
}

static double time_knn(Slot *slot) {
  seams_neighbours *nb = NULL;
  ensure_frame(slot);
  double t0 = now_sec();
  if (seams_neighbours_knn(slot->sys, 4, 3.5, 0, 1, &nb) != SEAMS_OK) {
    die(seams_handle_error(slot->sys));
  }
  int rows = seams_neighbours_count(nb);
  seams_neighbours_destroy(nb);
  g_sink += (double)rows;
  return now_sec() - t0;
}

static double time_steinhardt(Slot *slot) {
  buffers(slot);
  double t0 = now_sec();
  if (seams_steinhardt(slot->sys, 3.5, 0, 6, slot->ql, (size_t)slot->n, slot->bar,
                       (size_t)slot->n) != SEAMS_OK) {
    die(seams_handle_error(slot->sys));
  }
  double dt = now_sec() - t0;
  g_sink += slot->ql[0] + slot->bar[0];
  return dt;
}

static double time_chill(Slot *slot) {
  int counts[9];
  buffers(slot);
  double t0 = now_sec();
  if (seams_chill(slot->sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, slot->labels,
                  (size_t)slot->n, counts) != SEAMS_OK) {
    die(seams_handle_error(slot->sys));
  }
  double dt = now_sec() - t0;
  g_sink += (double)counts[SEAMS_LABEL_HEXAGONAL];
  return dt;
}

static double time_rings(Slot *slot) {
  seams_rings *rings = NULL;
  ensure_frame(slot);
  double t0 = now_sec();
  if (seams_rings_find(slot->sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    die(seams_handle_error(slot->sys));
  }
  int found = seams_rings_count(rings);
  seams_rings_destroy(rings);
  g_sink += (double)found;
  return now_sec() - t0;
}

int main(void) {
  char line[kLine];
  setvbuf(stdout, NULL, _IOLBF, 0);
  while (fgets(line, sizeof(line), stdin) != NULL) {
    int threads = 0;
    char op[32];
    char path[kLine];
    Slot *slot;
    double seconds;
    if (strncmp(line, "quit", 4) == 0) {
      break;
    }
    if (sscanf(line, "%d %31s %4095s", &threads, op, path) != 3) {
      die("command");
    }
    retune(threads);
    slot = slot_for(path);
    if (strcmp(op, "read") == 0) {
      seconds = time_read(slot);
    } else if (strcmp(op, "knn") == 0) {
      seconds = time_knn(slot);
    } else if (strcmp(op, "steinhardt") == 0) {
      seconds = time_steinhardt(slot);
    } else if (strcmp(op, "chill") == 0) {
      seconds = time_chill(slot);
    } else if (strcmp(op, "rings") == 0) {
      seconds = time_rings(slot);
    } else {
      die(op);
    }
    printf("%.9f\n", seconds * 1e3);
    fflush(stdout);
  }
  fprintf(stderr, "sink %.6f\n", g_sink);
  return 0;
}
