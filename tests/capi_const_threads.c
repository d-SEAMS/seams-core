/* Copyright (c) d-SEAMS developers
 * Const-handle calls from many threads on shared handles.
 *
 * The header allows const calls on one handle to overlap as long as no
 * non-const call does. Eight threads start together on one frame, one
 * neighbour list and one ring set, each mixing successful reads with
 * failing ones, and the first symmetry request races to fill the ring
 * label cache. Every result must match a single-thread reference, and a
 * failed const call must leave a message on the handle. Run it under
 * -Db_sanitize=thread to check that the overlap is race-free.
 */
#define _POSIX_C_SOURCE 200809L

#include <seams.h>

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NTHREADS 8
#define ROUNDS 200

static int g_n = 0;
static double g_xyz[3 * 4096];
static int g_types[4096];
static double g_box[9];

static seams_system *g_sys = NULL;
static seams_neighbours *g_nb = NULL;
static seams_rings *g_rings = NULL;

static size_t g_nrings = 0;
static int *g_ref_sizes = NULL;
static int *g_ref_sym = NULL;
static size_t g_noff = 0, g_nnz = 0;
static int *g_ref_off = NULL;
static int *g_ref_cols = NULL;

static pthread_barrier_t g_start;

static void ice_ic(int nrep) {
  const double a = 4.0 * 2.75 / sqrt(3.0);
  const double basis[8][3] = {{0, 0, 0},          {0.5, 0.5, 0},      {0.5, 0, 0.5},
                              {0, 0.5, 0.5},      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  int at = 0;
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (int b = 0; b < 8; ++b) {
          g_xyz[3 * at] = (i + basis[b][0]) * a;
          g_xyz[3 * at + 1] = (j + basis[b][1]) * a;
          g_xyz[3 * at + 2] = (k + basis[b][2]) * a;
          g_types[at] = 1 + at % 2;
          at += 1;
        }
      }
    }
  }
  g_n = at;
  memset(g_box, 0, sizeof g_box);
  g_box[3] = g_box[4] = g_box[5] = nrep * a;
}

static void need(seams_status status, const char *what) {
  if (status != SEAMS_OK) {
    fprintf(stderr, "%s: %s (%s)\n", what, seams_status_string(status),
            seams_last_error_message());
    exit(1);
  }
}

static void *reader(void *arg) {
  const size_t n = (size_t)g_n;
  long *bad = arg;
  double *xyz = malloc(sizeof(double) * 3 * n);
  int *ints = malloc(sizeof(int) * n);
  int *sizes = malloc(sizeof(int) * g_nrings);
  int *sym = malloc(sizeof(int) * g_nrings);
  int *off = malloc(sizeof(int) * g_noff);
  int *cols = malloc(sizeof(int) * g_nnz);
  double cell[9];
  if (xyz == NULL || ints == NULL || sizes == NULL || sym == NULL || off == NULL || cols == NULL) {
    *bad += 1;
    return NULL;
  }
  pthread_barrier_wait(&g_start);
  for (int r = 0; r < ROUNDS; ++r) {
    *bad += seams_rings_labels(g_rings, sizes, sym) != SEAMS_OK;
    *bad += memcmp(sizes, g_ref_sizes, sizeof(int) * g_nrings) != 0;
    *bad += memcmp(sym, g_ref_sym, sizeof(int) * g_nrings) != 0;

    *bad += seams_system_positions(g_sys, xyz, 3 * n) != SEAMS_OK;
    *bad += memcmp(xyz, g_xyz, sizeof(double) * 3 * n) != 0;
    *bad += seams_system_types(g_sys, ints, n) != SEAMS_OK;
    *bad += memcmp(ints, g_types, sizeof(int) * n) != 0;
    *bad += seams_system_ids(g_sys, ints, n) != SEAMS_OK;
    *bad += ints[n - 1] != g_n;
    *bad += seams_system_box(g_sys, cell, 9) != SEAMS_OK;
    *bad += fabs(cell[3] - g_box[3]) > 1e-12 * g_box[3];
    *bad += seams_system_count(g_sys) != g_n;

    size_t no = 0, nz = 0;
    *bad += seams_neighbours_copy(g_nb, off, g_noff, &no, cols, g_nnz, &nz) != SEAMS_OK;
    *bad += no != g_noff || nz != g_nnz;
    *bad += memcmp(off, g_ref_off, sizeof(int) * g_noff) != 0;
    *bad += memcmp(cols, g_ref_cols, sizeof(int) * g_nnz) != 0;

    *bad += seams_system_positions(g_sys, xyz, 3 * n - 1) != SEAMS_ERR_RANGE;
    *bad += seams_last_error_message()[0] == '\0';
    *bad += seams_neighbours_copy(g_nb, off, 0, &no, cols, g_nnz, &nz) != SEAMS_ERR_RANGE;
    *bad += seams_last_error_message()[0] == '\0';
  }
  free(xyz);
  free(ints);
  free(sizes);
  free(sym);
  free(off);
  free(cols);
  return NULL;
}

int main(void) {
  ice_ic(4);
  const size_t n = (size_t)g_n;
  g_sys = seams_system_create();
  if (g_sys == NULL) {
    return 1;
  }
  need(seams_system_set_atoms(g_sys, g_xyz, n, g_box, 9, g_types, NULL), "set_atoms");
  need(seams_neighbours_build(g_sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &g_nb), "neighbours_build");
  need(seams_rings_find(g_sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &g_rings), "rings_find");
  need(seams_rings_nrings(g_rings, &g_nrings), "rings_nrings");
  need(seams_neighbours_copy(g_nb, NULL, 0, &g_noff, NULL, 0, &g_nnz), "neighbours_copy size");
  g_ref_sizes = malloc(sizeof(int) * g_nrings);
  g_ref_sym = malloc(sizeof(int) * g_nrings);
  g_ref_off = malloc(sizeof(int) * g_noff);
  g_ref_cols = malloc(sizeof(int) * g_nnz);
  if (g_nrings == 0 || g_ref_sizes == NULL || g_ref_sym == NULL || g_ref_off == NULL ||
      g_ref_cols == NULL) {
    fprintf(stderr, "no rings or no memory\n");
    return 1;
  }
  need(seams_neighbours_copy(g_nb, g_ref_off, g_noff, &g_noff, g_ref_cols, g_nnz, &g_nnz),
       "neighbours_copy");

  /* A second ring set fills its cache on one thread; the threaded run
   * fills the shared one's cache concurrently, and the two must agree. */
  seams_rings *serial = NULL;
  need(seams_rings_find(g_sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &serial), "rings_find");
  need(seams_rings_labels(serial, g_ref_sizes, g_ref_sym), "rings_labels");
  seams_rings_destroy(serial);

  pthread_t threads[NTHREADS];
  long bad[NTHREADS] = {0};
  pthread_barrier_init(&g_start, NULL, NTHREADS);
  for (int t = 0; t < NTHREADS; ++t) {
    if (pthread_create(&threads[t], NULL, reader, &bad[t]) != 0) {
      fprintf(stderr, "pthread_create failed\n");
      return 1;
    }
  }
  long total = 0;
  for (int t = 0; t < NTHREADS; ++t) {
    pthread_join(threads[t], NULL);
    total += bad[t];
  }
  pthread_barrier_destroy(&g_start);
  /* The handle copy aliases the handle, so it is read once the threads are done. */
  total += seams_handle_error(g_sys)[0] == '\0';
  total += seams_handle_error(g_nb)[0] == '\0';

  printf("capi_const_threads: %d threads x %d rounds, %zu atoms, %zu rings, %ld mismatches\n",
         NTHREADS, ROUNDS, n, g_nrings, total);
  seams_rings_destroy(g_rings);
  seams_neighbours_destroy(g_nb);
  seams_system_destroy(g_sys);
  free(g_ref_sizes);
  free(g_ref_sym);
  free(g_ref_off);
  free(g_ref_cols);
  return total == 0 ? 0 : 1;
}
