/* Per-function contract of the stable C API, compiled as C11.
 *
 * Every function in seams.h is called with bad input here: NULL pointers,
 * zero and short lengths, a frame with a wrong box, a frame with NaN
 * coordinates, and non-finite parameters. Each call must return the
 * documented status and, on failure, a non-empty error message.
 *
 * Usage: capi_contract <symbols-list> <source-root>
 * Every name in the symbols list must be exercised, or the test fails.
 */
#define _POSIX_C_SOURCE 200809L

#include "seams.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_MAX 1024
#define N_COVER 256

static int g_fails = 0;
static int g_checks = 0;
static const char *g_covered[N_COVER];
static int g_ncovered = 0;

static void cover(const char *name) {
  for (int i = 0; i < g_ncovered; ++i) {
    if (strcmp(g_covered[i], name) == 0) {
      return;
    }
  }
  if (g_ncovered < N_COVER) {
    g_covered[g_ncovered++] = name;
  }
}

static void fault(const char *fn, const char *what, const char *detail) {
  fprintf(stderr, "FAIL %s: %s (%s)\n", fn, what, detail);
  g_fails += 1;
}

static void truth(int ok, const char *fn, const char *what) {
  g_checks += 1;
  cover(fn);
  if (!ok) {
    fault(fn, what, seams_last_error_message());
  }
}

/* The status must match, and the thread-local message must be empty on
 * success and non-empty on failure. */
static void expect(seams_status got, seams_status want, const char *fn, const char *what) {
  g_checks += 1;
  cover(fn);
  const char *msg = seams_last_error_message();
  if (got != want) {
    char detail[512];
    snprintf(detail, sizeof detail, "got %s, want %s: %s", seams_status_string(got),
             seams_status_string(want), msg);
    fault(fn, what, detail);
    return;
  }
  if (want != SEAMS_OK && (msg == NULL || msg[0] == '\0')) {
    fault(fn, what, "no error message");
  }
  if (want == SEAMS_OK && msg != NULL && msg[0] != '\0') {
    fault(fn, what, "stale error message after success");
  }
}

/* ------------------------------------------------------------------ */
/* Frames                                                              */

static int g_n = 0;
static double g_xyz[3 * N_MAX];
static double g_box[9];
static int g_types[N_MAX];

static void ice_ic(void) {
  const int nrep = 2;
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
          g_types[at] = 1;
          at += 1;
        }
      }
    }
  }
  g_n = at;
  memset(g_box, 0, sizeof g_box);
  g_box[3] = g_box[4] = g_box[5] = nrep * a;
}

static seams_system *g_good = NULL;   /* holds the ice frame */
static seams_system *g_empty = NULL;  /* never filled */
static seams_system *g_nan = NULL;    /* only a NaN-coordinate load was tried */
static seams_system *g_badbox = NULL; /* only a wrong-box load was tried */

static void frames(void) {
  ice_ic();
  g_good = seams_system_create();
  g_empty = seams_system_create();
  g_nan = seams_system_create();
  g_badbox = seams_system_create();
  truth(g_good != NULL && g_empty != NULL && g_nan != NULL && g_badbox != NULL,
        "seams_system_create", "create returns handles");
  expect(seams_system_set_atoms(g_good, g_xyz, (size_t)g_n, g_box, 9, g_types, NULL), SEAMS_OK,
         "seams_system_set_atoms", "ice frame");

  double xyz[3 * N_MAX];
  memcpy(xyz, g_xyz, sizeof(double) * 3 * (size_t)g_n);
  xyz[7] = NAN;
  expect(seams_system_set_atoms(g_nan, xyz, (size_t)g_n, g_box, 9, NULL, NULL), SEAMS_ERR_ARG,
         "seams_system_set_atoms", "NaN coordinate");

  double box[9];
  memcpy(box, g_box, sizeof box);
  box[4] = -box[4];
  expect(seams_system_set_atoms(g_badbox, g_xyz, (size_t)g_n, box, 9, NULL, NULL),
         SEAMS_ERR_BOX, "seams_system_set_atoms", "negative edge");
}

static void frame_contract(void) {
  const char *fn = "seams_system_set_atoms";
  seams_system *sys = seams_system_create();
  double xyz[3 * N_MAX];
  double box[9];
  memcpy(xyz, g_xyz, sizeof(double) * 3 * (size_t)g_n);
  memcpy(box, g_box, sizeof box);
  const size_t n = (size_t)g_n;

  expect(seams_system_set_atoms(NULL, xyz, n, box, 9, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_set_atoms(sys, NULL, n, box, 9, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL xyz");
  expect(seams_system_set_atoms(sys, xyz, n, NULL, 9, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL box");
  expect(seams_system_set_atoms(sys, xyz, 0, box, 9, NULL, NULL), SEAMS_ERR_RANGE, fn, "zero atoms");
  expect(seams_system_set_atoms(sys, xyz, n, box, 0, NULL, NULL), SEAMS_ERR_RANGE, fn, "zero box length");
  expect(seams_system_set_atoms(sys, xyz, n, box, 6, NULL, NULL), SEAMS_ERR_RANGE, fn, "short box length");
  expect(seams_system_set_atoms(sys, xyz, (size_t)INT32_MAX + 1u, box, 9, NULL, NULL),
         SEAMS_ERR_RANGE, fn, "count past int");
  box[3] = 0.0;
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, NULL, NULL), SEAMS_ERR_BOX, fn, "zero edge");
  box[3] = NAN;
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, NULL, NULL), SEAMS_ERR_BOX, fn, "NaN edge");
  box[3] = g_box[3];
  box[7] = INFINITY;
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, NULL, NULL), SEAMS_ERR_BOX, fn, "infinite tilt");
  box[7] = 0.0;
  xyz[0] = INFINITY;
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, NULL, NULL), SEAMS_ERR_ARG, fn, "infinite coordinate");
  xyz[0] = g_xyz[0];
  int ids[N_MAX];
  for (int i = 0; i < g_n; ++i) {
    ids[i] = i + 1;
  }
  ids[1] = 1;
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, NULL, ids), SEAMS_ERR_ARG, fn, "repeated ids");
  truth(seams_system_count(sys) == 0, "seams_system_count", "failed loads leave no frame");
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, g_types, NULL), SEAMS_OK, fn, "valid frame");
  xyz[5] = NAN;
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, NULL, NULL), SEAMS_ERR_ARG, fn, "NaN after a frame");
  truth(seams_system_count(sys) == g_n, "seams_system_count", "failed load keeps the earlier frame");
  xyz[5] = g_xyz[5];

  fn = "seams_system_set_frame";
  expect(seams_system_set_frame(NULL, xyz, g_n, box, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_set_frame(sys, NULL, g_n, box, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL xyz");
  expect(seams_system_set_frame(sys, xyz, g_n, NULL, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL box");
  expect(seams_system_set_frame(sys, xyz, 0, box, NULL, NULL), SEAMS_ERR_RANGE, fn, "zero atoms");
  expect(seams_system_set_frame(sys, xyz, -3, box, NULL, NULL), SEAMS_ERR_RANGE, fn, "negative atoms");
  box[5] = -1.0;
  expect(seams_system_set_frame(sys, xyz, g_n, box, NULL, NULL), SEAMS_ERR_BOX, fn, "negative edge");
  box[5] = g_box[5];
  xyz[2] = NAN;
  expect(seams_system_set_frame(sys, xyz, g_n, box, NULL, NULL), SEAMS_ERR_ARG, fn, "NaN coordinate");
  xyz[2] = g_xyz[2];
  expect(seams_system_set_frame(sys, xyz, g_n, box, NULL, NULL), SEAMS_OK, fn, "valid frame");

  fn = "seams_system_set_ortho";
  const double l = g_box[3];
  expect(seams_system_set_ortho(NULL, xyz, g_n, 0, 0, 0, l, l, l, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_set_ortho(sys, NULL, g_n, 0, 0, 0, l, l, l, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL xyz");
  expect(seams_system_set_ortho(sys, xyz, 0, 0, 0, 0, l, l, l, NULL, NULL), SEAMS_ERR_RANGE, fn, "zero atoms");
  expect(seams_system_set_ortho(sys, xyz, g_n, 0, 0, 0, 0, l, l, NULL, NULL), SEAMS_ERR_BOX, fn, "zero edge");
  expect(seams_system_set_ortho(sys, xyz, g_n, NAN, 0, 0, l, l, l, NULL, NULL), SEAMS_ERR_BOX, fn, "NaN origin");
  expect(seams_system_set_ortho(sys, xyz, g_n, 0, 0, 0, l, NAN, l, NULL, NULL), SEAMS_ERR_BOX, fn, "NaN edge");
  xyz[4] = NAN;
  expect(seams_system_set_ortho(sys, xyz, g_n, 0, 0, 0, l, l, l, NULL, NULL), SEAMS_ERR_ARG, fn, "NaN coordinate");
  xyz[4] = g_xyz[4];
  expect(seams_system_set_ortho(sys, xyz, g_n, 0, 0, 0, l, l, l, NULL, NULL), SEAMS_OK, fn, "valid frame");

  fn = "seams_system_count";
  truth(seams_system_count(NULL) == 0, fn, "NULL is zero");
  truth(seams_system_count(g_empty) == 0, fn, "no frame is zero");
  truth(seams_system_count(sys) == g_n, fn, "frame count");

  fn = "seams_system_natoms";
  size_t count = 12345;
  expect(seams_system_natoms(NULL, &count), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_natoms(sys, NULL), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_system_natoms(g_empty, &count), SEAMS_ERR_STATE, fn, "no frame");
  expect(seams_system_natoms(g_nan, &count), SEAMS_ERR_STATE, fn, "NaN load left no frame");
  expect(seams_system_natoms(g_badbox, &count), SEAMS_ERR_STATE, fn, "wrong box left no frame");
  truth(count == 12345, fn, "failure leaves the output alone");
  expect(seams_system_natoms(sys, &count), SEAMS_OK, fn, "frame");
  truth(count == n, fn, "count matches");

  fn = "seams_system_destroy";
  seams_system_destroy(NULL);
  cover(fn);
  seams_system_destroy(sys);
  seams_system_destroy(sys);
  truth(seams_last_error_message()[0] != '\0', fn, "double destroy is reported");
}

static void frame_io_contract(void) {
  const size_t n = (size_t)g_n;
  seams_system *sys = seams_system_create();
  double xyz[3 * N_MAX];
  double out[3 * N_MAX];
  double box[9];
  double cell[9];
  int ids[N_MAX];
  int got[N_MAX];
  memcpy(xyz, g_xyz, sizeof(double) * 3 * n);
  memcpy(box, g_box, sizeof box);

  const char *fn = "seams_system_set_positions";
  expect(seams_system_set_positions(NULL, xyz, n, NULL, 0), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_set_positions(sys, NULL, n, NULL, 0), SEAMS_ERR_NULL, fn, "NULL xyz");
  expect(seams_system_set_positions(sys, xyz, n, NULL, 0), SEAMS_ERR_STATE, fn, "no frame");
  for (int i = 0; i < g_n; ++i) {
    ids[i] = 100 + 3 * (g_n - i);
  }
  expect(seams_system_set_atoms(sys, xyz, n, box, 9, g_types, ids), SEAMS_OK,
         "seams_system_set_atoms", "explicit ids");
  expect(seams_system_set_positions(sys, xyz, n - 1, NULL, 0), SEAMS_ERR_RANGE, fn, "count changed");
  expect(seams_system_set_positions(sys, xyz, n, box, 6), SEAMS_ERR_RANGE, fn, "short box length");
  box[4] = -1.0;
  expect(seams_system_set_positions(sys, xyz, n, box, 9), SEAMS_ERR_BOX, fn, "negative edge");
  box[4] = g_box[4];
  for (size_t k = 0; k < 3 * n; ++k) {
    xyz[k] = g_xyz[k] + 0.125;
  }
  xyz[3 * n - 1] = NAN;
  expect(seams_system_set_positions(sys, xyz, n, NULL, 0), SEAMS_ERR_ARG, fn, "NaN coordinate");
  expect(seams_system_positions(sys, out, 3 * n), SEAMS_OK, "seams_system_positions", "frame");
  truth(memcmp(out, g_xyz, sizeof(double) * 3 * n) == 0, fn, "failed reload keeps the frame");
  xyz[3 * n - 1] = g_xyz[3 * n - 1] + 0.125;
  box[3] += 1.0;
  expect(seams_system_set_positions(sys, xyz, n, box, 9), SEAMS_OK, fn, "new positions and cell");
  expect(seams_system_set_positions(sys, xyz, n, NULL, 0), SEAMS_OK, fn, "positions only");

  fn = "seams_system_positions";
  expect(seams_system_positions(NULL, out, 3 * n), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_positions(sys, NULL, 3 * n), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_system_positions(g_empty, out, 3 * n), SEAMS_ERR_STATE, fn, "no frame");
  expect(seams_system_positions(sys, out, 3 * n - 1), SEAMS_ERR_RANGE, fn, "short buffer");
  truth(seams_handle_error(sys)[0] != '\0', fn, "const failure reaches the handle");
  expect(seams_system_positions(sys, out, 3 * n), SEAMS_OK, fn, "frame");
  truth(memcmp(out, xyz, sizeof(double) * 3 * n) == 0, fn, "round trip");

  fn = "seams_system_box";
  expect(seams_system_box(NULL, cell, 9), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_box(sys, NULL, 9), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_system_box(sys, cell, 6), SEAMS_ERR_RANGE, fn, "short buffer");
  expect(seams_system_box(g_empty, cell, 9), SEAMS_ERR_STATE, fn, "no frame");
  expect(seams_system_box(sys, cell, 9), SEAMS_OK, fn, "frame");
  int same = 1;
  for (int k = 0; k < 9; ++k) {
    same &= fabs(cell[k] - box[k]) <= 1e-12 * (1.0 + fabs(box[k]));
  }
  truth(same, fn, "cell round trip");
  const double tilted[9] = {-1.0, 0.5, 2.0, 9.0, 8.0, 7.0, 1.5, -0.75, 0.25};
  expect(seams_system_set_positions(sys, xyz, n, tilted, 9), SEAMS_OK,
         "seams_system_set_positions", "triclinic cell");
  expect(seams_system_box(sys, cell, 9), SEAMS_OK, fn, "triclinic frame");
  same = 1;
  for (int k = 0; k < 9; ++k) {
    same &= fabs(cell[k] - tilted[k]) <= 1e-12 * (1.0 + fabs(tilted[k]));
  }
  truth(same, fn, "triclinic round trip");

  fn = "seams_system_types";
  expect(seams_system_types(NULL, got, n), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_types(sys, NULL, n), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_system_types(g_empty, got, n), SEAMS_ERR_STATE, fn, "no frame");
  expect(seams_system_types(sys, got, n - 1), SEAMS_ERR_RANGE, fn, "short buffer");
  expect(seams_system_types(sys, got, n), SEAMS_OK, fn, "frame");
  truth(memcmp(got, g_types, sizeof(int) * n) == 0, fn, "round trip");

  fn = "seams_system_ids";
  expect(seams_system_ids(NULL, got, n), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_ids(sys, NULL, n), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_system_ids(g_empty, got, n), SEAMS_ERR_STATE, fn, "no frame");
  expect(seams_system_ids(sys, got, n - 1), SEAMS_ERR_RANGE, fn, "short buffer");
  expect(seams_system_ids(sys, got, n), SEAMS_OK, fn, "explicit ids");
  truth(memcmp(got, ids, sizeof(int) * n) == 0, fn, "explicit ids round trip");
  expect(seams_system_set_atoms(sys, g_xyz, n, g_box, 9, NULL, NULL), SEAMS_OK,
         "seams_system_set_atoms", "default ids after explicit ids");
  expect(seams_system_ids(sys, got, n), SEAMS_OK, fn, "default ids");
  int sequential = 1;
  for (int i = 0; i < g_n; ++i) {
    sequential &= got[i] == i + 1;
  }
  truth(sequential, fn, "default ids are 1..n");
  int labels[N_MAX];
  expect(seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n,
                            NULL),
         SEAMS_OK, "seams_chill_plus_n", "after a default-id reload");
  int cubic = 0;
  for (int i = 0; i < g_n; ++i) {
    cubic += labels[i] == SEAMS_LABEL_CUBIC;
  }
  truth(cubic == g_n, "seams_chill_plus_n", "reloaded ice is all cubic");
  expect(seams_system_set_positions(sys, g_xyz, n, NULL, 0), SEAMS_OK,
         "seams_system_set_positions", "same positions");
  expect(seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n,
                            NULL),
         SEAMS_OK, "seams_chill_plus_n", "after set_positions");
  cubic = 0;
  for (int i = 0; i < g_n; ++i) {
    cubic += labels[i] == SEAMS_LABEL_CUBIC;
  }
  truth(cubic == g_n, "seams_chill_plus_n", "labels survive a position reload");
  truth(seams_handle_error(sys)[0] == '\0', "seams_handle_error", "non-const call clears the slot");

  expect(seams_system_stage_positions(NULL, 1), SEAMS_ERR_NULL, "seams_system_stage_positions",
         "NULL sys");
  expect(seams_system_stage_positions(g_empty, 1), SEAMS_OK, "seams_system_stage_positions",
         "no frame");
  expect(seams_system_stage_positions(g_empty, 0), SEAMS_OK, "seams_system_stage_positions",
         "off again");
  for (int stage = 0; stage < 2; ++stage) {
    fn = stage ? "seams_system_set_positions staged" : "seams_system_set_positions";
    expect(seams_system_stage_positions(sys, stage), SEAMS_OK, "seams_system_stage_positions",
           stage ? "on" : "off");
    expect(seams_system_set_atoms(sys, g_xyz, n, g_box, 9, g_types, NULL), SEAMS_OK,
           "seams_system_set_atoms", "before reloads");
    expect(seams_system_set_positions(sys, g_xyz, n, NULL, 0), SEAMS_OK, fn, "first reload");
    for (size_t k = 0; k < 3 * n; ++k) {
      xyz[k] = g_xyz[k] + 0.0625;
    }
    xyz[4] = INFINITY;
    expect(seams_system_set_positions(sys, xyz, n, NULL, 0), SEAMS_ERR_ARG, fn,
           "infinite coordinate after a reload");
    expect(seams_system_positions(sys, out, 3 * n), SEAMS_OK, "seams_system_positions", "frame");
    truth(memcmp(out, g_xyz, sizeof(double) * 3 * n) == 0, fn, "second failed reload keeps the frame");
    xyz[4] = g_xyz[4] + 0.0625;
    for (int r = 0; r < 3; ++r) {
      expect(seams_system_set_positions(sys, r % 2 ? g_xyz : xyz, n, NULL, 0), SEAMS_OK, fn,
             "repeated reloads");
    }
    expect(seams_system_positions(sys, out, 3 * n), SEAMS_OK, "seams_system_positions", "frame");
    truth(memcmp(out, xyz, sizeof(double) * 3 * n) == 0, fn, "last of repeated reloads is stored");
    int retyped[N_MAX];
    for (int i = 0; i < g_n; ++i) {
      retyped[i] = 1 + i % 3;
    }
    expect(seams_system_set_atoms(sys, g_xyz, n, g_box, 9, retyped, ids), SEAMS_OK,
           "seams_system_set_atoms", "new types and ids");
    expect(seams_system_set_positions(sys, xyz, n, NULL, 0), SEAMS_OK, fn, "after new types");
    expect(seams_system_types(sys, got, n), SEAMS_OK, "seams_system_types", "after new types");
    truth(memcmp(got, retyped, sizeof(int) * n) == 0, fn, "reload keeps the latest types");
    expect(seams_system_ids(sys, got, n), SEAMS_OK, "seams_system_ids", "after new ids");
    truth(memcmp(got, ids, sizeof(int) * n) == 0, fn, "reload keeps the latest ids");
    expect(seams_system_positions(sys, out, 3 * n), SEAMS_OK, "seams_system_positions", "frame");
    truth(memcmp(out, xyz, sizeof(double) * 3 * n) == 0, fn, "reload after new types is stored");
  }
  seams_system_destroy(sys);
}

static void read_contract(const char *root) {
  char path[4096];
  snprintf(path, sizeof path, "%s/input/traj/agi_water_tiny.lammpstrj", root);
  const char *fn = "seams_system_read";
  seams_system *sys = seams_system_create();
  expect(seams_system_read(NULL, path, 1, -1), SEAMS_ERR_NULL, fn, "NULL sys");
  expect(seams_system_read(sys, NULL, 1, -1), SEAMS_ERR_NULL, fn, "NULL path");
  expect(seams_system_read(sys, path, 0, -1), SEAMS_ERR_RANGE, fn, "frame 0");
  expect(seams_system_read(sys, path, -2, -1), SEAMS_ERR_RANGE, fn, "negative frame");
  expect(seams_system_read(sys, "no/such/file.lammpstrj", 1, -1), SEAMS_ERR_ARG, fn, "missing file");
  expect(seams_system_read(sys, path, 1, -1), SEAMS_OK, fn, "first frame");
  const int n = seams_system_count(sys);
  expect(seams_system_read(sys, path, 100000, -1), SEAMS_ERR_RANGE, fn, "past the last frame");
  truth(seams_system_count(sys) == n, fn, "a failed read keeps the frame");

  fn = "seams_path_nframes";
  size_t nframes = 77;
  expect(seams_path_nframes(NULL, &nframes), SEAMS_ERR_NULL, fn, "NULL path");
  expect(seams_path_nframes(path, NULL), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_path_nframes("no/such/file.lammpstrj", &nframes), SEAMS_ERR_ARG, fn, "missing file");
  truth(nframes == 77, fn, "failure leaves the output alone");
  expect(seams_path_nframes(path, &nframes), SEAMS_OK, fn, "count");
  truth(nframes >= 1, fn, "at least one frame");
  seams_system_destroy(sys);
}

/* ------------------------------------------------------------------ */
/* Metadata                                                            */

static void meta_contract(void) {
  truth(seams_abi() == SEAMS_ABI, "seams_abi", "matches header");
  truth(seams_version() == SEAMS_VERSION_NUMBER, "seams_version", "matches header");
  truth(seams_version_major() == SEAMS_VERSION_MAJOR, "seams_version_major", "matches header");
  truth(seams_version_minor() == SEAMS_VERSION_MINOR, "seams_version_minor", "matches header");
  truth(seams_version_patch() == SEAMS_VERSION_PATCH, "seams_version_patch", "matches header");
  truth(strcmp(seams_version_string(), SEAMS_VERSION_STRING) == 0, "seams_version_string",
        "matches header");
  for (int code = SEAMS_OK; code <= SEAMS_ERR_INTERNAL; ++code) {
    const char *text = seams_status_string((seams_status)code);
    truth(text != NULL && text[0] != '\0' && strcmp(text, "unknown status") != 0,
          "seams_status_string", "named code");
  }
  truth(strcmp(seams_status_string((seams_status)99), "unknown status") == 0,
        "seams_status_string", "unknown code");
  for (int label = SEAMS_LABEL_CUBIC; label <= SEAMS_LABEL_REHEX; ++label) {
    truth(strcmp(seams_label_name(label), "unknown") != 0, "seams_label_name", "named label");
  }
  truth(strcmp(seams_label_name(-1), "unknown") == 0, "seams_label_name", "negative label");
  truth(strcmp(seams_label_name(9), "unknown") == 0, "seams_label_name", "label past the end");
  for (int sym = SEAMS_RING_UNCLASSIFIED; sym <= SEAMS_RING_OTHER; ++sym) {
    truth(strcmp(seams_ring_symmetry_name(sym), "unknown") != 0, "seams_ring_symmetry_name",
          "named symmetry");
  }
  truth(strcmp(seams_ring_symmetry_name(-1), "unknown") == 0, "seams_ring_symmetry_name",
        "negative symmetry");
  truth(strcmp(seams_ring_symmetry_name(12), "unknown") == 0, "seams_ring_symmetry_name",
        "symmetry past the end");

  expect(seams_set_num_threads(-1), SEAMS_ERR_RANGE, "seams_set_num_threads", "negative width");
  truth(seams_last_error() == seams_last_error_message() ||
            strcmp(seams_last_error(), seams_last_error_message()) == 0,
        "seams_last_error", "same string as seams_last_error_message");
  truth(seams_last_error_message()[0] != '\0', "seams_last_error_message", "set after a failure");
  truth(strcmp(seams_handle_error(NULL), seams_last_error_message()) == 0, "seams_handle_error",
        "NULL reads the thread-local string");
  expect(seams_set_num_threads(1), SEAMS_OK, "seams_set_num_threads", "one thread");
  truth(seams_last_error_message()[0] == '\0', "seams_last_error_message", "cleared by success");
  truth(seams_get_num_threads() == 1, "seams_get_num_threads", "reads the stored width");
  expect(seams_set_num_threads(0), SEAMS_OK, "seams_set_num_threads", "restore the default");
  truth(seams_get_num_threads() >= 1, "seams_get_num_threads", "positive default");
  expect(seams_set_num_threads(1), SEAMS_OK, "seams_set_num_threads", "serial for the rest");

  /* A handle keeps its own message after a later call on another handle. */
  seams_neighbours *nb = NULL;
  expect(seams_neighbours_build(g_empty, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb), SEAMS_ERR_STATE,
         "seams_neighbours_build", "no frame");
  expect(seams_system_natoms(g_good, &(size_t){0}), SEAMS_OK, "seams_system_natoms", "frame");
  truth(seams_handle_error(g_empty)[0] != '\0', "seams_handle_error", "message survives on the handle");
  truth(seams_handle_error(g_good)[0] == '\0', "seams_handle_error", "clean handle has no message");
}

/* ------------------------------------------------------------------ */
/* Analysis entry points that take a frame                             */

/* Probe wrappers. Each calls one entry point with sys and cutoff and a
 * caller buffer of len entries. When g_null_out is set the primary
 * output pointer is NULL. g_graph is passed where the call takes one. */
static int g_null_out = 0;
static int g_graph = SEAMS_GRAPH_CUTOFF;
static int g_ibuf[N_MAX];
static int g_ibuf2[N_MAX];
static double g_dbuf[N_MAX];
static double g_dbuf2[N_MAX];
static double g_dbuf3[N_MAX];

typedef seams_status (*probe_fn)(seams_system *sys, double cutoff, size_t len);

enum {
  P_CUTOFF = 1u,  /* the second argument is a cutoff or radius */
  P_LEN = 2u,     /* the third argument is a checked buffer length */
  P_NULLOUT = 4u, /* NULL primary output is rejected */
  P_GRAPH = 8u    /* takes a graph code */
};

static seams_status p_neighbours_build(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  seams_neighbours *nb = NULL;
  seams_status rc = seams_neighbours_build(sys, cutoff, 0, g_graph, g_null_out ? NULL : &nb);
  seams_neighbours_destroy(nb);
  return rc;
}
static seams_status p_neighbours_knn(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  seams_neighbours *nb = NULL;
  seams_status rc = seams_neighbours_knn(sys, 4, cutoff, 0, 1, g_null_out ? NULL : &nb);
  seams_neighbours_destroy(nb);
  return rc;
}
static seams_status p_chill(seams_system *sys, double cutoff, size_t len) {
  return seams_chill(sys, cutoff, 0, g_graph, SEAMS_CHILL_PUBLISHED, g_null_out ? NULL : g_ibuf,
                     len, NULL);
}
static seams_status p_chill_plus(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  return seams_chill_plus(sys, cutoff, 0, g_graph, SEAMS_CHILL_PUBLISHED,
                          g_null_out ? NULL : g_ibuf, NULL);
}
static seams_status p_chill_plus_n(seams_system *sys, double cutoff, size_t len) {
  int counts[9];
  return seams_chill_plus_n(sys, cutoff, 0, g_graph, SEAMS_CHILL_PUBLISHED,
                            g_null_out ? NULL : g_ibuf, len, counts);
}
static seams_status p_ql(seams_system *sys, double cutoff, size_t len) {
  return seams_ql(sys, cutoff, 0, 6, g_graph, g_null_out ? NULL : g_dbuf, len,
                  g_null_out ? NULL : g_dbuf2, len);
}
static seams_status p_q4q6(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  return seams_q4q6(sys, cutoff, 0, g_null_out ? NULL : g_dbuf, g_null_out ? NULL : g_dbuf2);
}
static seams_status p_q4q6_n(seams_system *sys, double cutoff, size_t len) {
  return seams_q4q6_n(sys, cutoff, 0, g_null_out ? NULL : g_dbuf, len,
                      g_null_out ? NULL : g_dbuf2, len);
}
static seams_status p_q6bar(seams_system *sys, double cutoff, size_t len) {
  double mean = 0, median = 0;
  int n_ice = 0;
  if (g_null_out) {
    return seams_q6bar(sys, cutoff, 0, 0.36, NULL, 0, NULL, NULL, NULL);
  }
  return seams_q6bar(sys, cutoff, 0, 0.36, g_dbuf, len, &mean, &median, &n_ice);
}
static seams_status p_steinhardt(seams_system *sys, double cutoff, size_t len) {
  return seams_steinhardt(sys, cutoff, 0, 4, g_null_out ? NULL : g_dbuf, len, NULL, 0);
}
static seams_status p_symmetry(seams_system *sys, double cutoff, size_t len) {
  (void)cutoff;
  return seams_system_symmetry(sys, 0.0, g_null_out ? NULL : g_ibuf, len);
}
static seams_status p_rings_find(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  seams_rings *rings = NULL;
  seams_status rc = seams_rings_find(sys, cutoff, 0, g_graph, 6, g_null_out ? NULL : &rings);
  seams_rings_destroy(rings);
  return rc;
}
static seams_status p_cages_find(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  seams_cage_set *cages = NULL;
  seams_status rc = seams_cages_find(sys, cutoff, 0, g_graph, 0, 0, g_null_out ? NULL : &cages);
  seams_cages_destroy(cages);
  return rc;
}
static seams_status p_cages(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  seams_cage_set *cages = NULL;
  seams_status rc = seams_cages(sys, cutoff, 0, g_graph, 0, 0, g_null_out ? NULL : &cages);
  seams_cages_destroy(cages);
  return rc;
}
static seams_status p_cages_signature(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  size_t n_cages = 0, n_atoms = 0;
  return seams_cages_signature(sys, cutoff, 0, g_graph, "hc", g_null_out ? NULL : &n_cages,
                               g_null_out ? NULL : &n_atoms);
}
static seams_status p_signature_cages(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  size_t n_cages = 0, n_atoms = 0;
  return seams_signature_cages(sys, cutoff, 0, g_graph, "ddc", g_null_out ? NULL : &n_cages,
                               g_null_out ? NULL : &n_atoms);
}
static seams_status p_brine(seams_system *sys, double cutoff, size_t len) {
  (void)cutoff;
  seams_brine *front = NULL;
  seams_status rc =
      seams_brine_analyse(sys, g_ibuf2, len, 0, 0, 4, 2, g_null_out ? NULL : &front);
  seams_brine_destroy(front);
  return rc;
}
static seams_status p_rdf(seams_system *sys, double cutoff, size_t len) {
  size_t n_out = 0;
  double volume = 0;
  return seams_rdf(sys, cutoff, 20, 0, 0, g_dbuf, g_dbuf2, g_ibuf, len,
                   g_null_out ? NULL : &n_out, &volume);
}
static seams_status p_running_cn(seams_system *sys, double cutoff, size_t len) {
  size_t n_out = 0;
  return seams_running_cn(sys, cutoff, 20, 0, 0, g_dbuf, g_dbuf2, len,
                          g_null_out ? NULL : &n_out);
}
static seams_status p_coordination(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  double cn = 0;
  return seams_coordination(sys, cutoff, 0, 1, 1, g_null_out ? NULL : &cn);
}
static seams_status p_hbonds(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  int n_bonds = 0;
  return seams_hbonds(sys, cutoff, 1, 2, 2.42, 30.0, g_null_out ? NULL : &n_bonds);
}
static seams_status p_density(seams_system *sys, double cutoff, size_t len) {
  (void)cutoff;
  size_t n_out = 0;
  return seams_density_axis(sys, 0, 10, 2, g_dbuf, g_dbuf3, len, g_null_out ? NULL : &n_out);
}
static seams_status p_fingerprint(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  seams_topo *fp = NULL;
  seams_status rc = seams_fingerprint(sys, cutoff, 0, 4, g_graph, 2, 0, g_null_out ? NULL : &fp);
  seams_fingerprint_destroy(fp);
  return rc;
}
static seams_status p_contacts(seams_system *sys, double cutoff, size_t len) {
  (void)cutoff;
  (void)len;
  size_t n_pairs = 0;
  int n_cation = 0, n_anion = 0;
  if (g_null_out) {
    return seams_contact_pairs(sys, "1=cationHead", NULL, NULL, NULL);
  }
  return seams_contact_pairs(sys, "1=cationHead", &n_pairs, &n_cation, &n_anion);
}
static seams_status p_domains(seams_system *sys, double cutoff, size_t len) {
  (void)len;
  int n_subset = 0, largest = 0;
  double percolation = 0;
  if (g_null_out) {
    return seams_domains(sys, cutoff, "1=waterO", 0, NULL, NULL, NULL);
  }
  return seams_domains(sys, cutoff, "1=waterO", 0, &n_subset, &largest, &percolation);
}
static seams_status p_prisms(seams_system *sys, double cutoff, size_t len) {
  return seams_prisms(sys, cutoff, 0, g_null_out ? NULL : g_ibuf, g_ibuf2, len);
}
static seams_status p_monolayer(seams_system *sys, double cutoff, size_t len) {
  return seams_monolayer(sys, cutoff, 0, 6, 100.0, g_null_out ? NULL : g_ibuf, len);
}
static seams_status p_site_table(seams_system *sys, double cutoff, size_t len) {
  (void)cutoff;
  return seams_site_table(sys, "1=waterO", g_null_out ? NULL : g_ibuf, len);
}

/* need is the shortest accepted buffer length for P_LEN probes. */
static void probe(const char *fn, probe_fn f, unsigned flags, size_t need) {
  const size_t len = need > 0 ? need : (size_t)g_n;
  expect(f(NULL, 3.5, len), SEAMS_ERR_NULL, fn, "NULL frame");
  expect(f(g_empty, 3.5, len), SEAMS_ERR_STATE, fn, "frame never filled");
  expect(f(g_nan, 3.5, len), SEAMS_ERR_STATE, fn, "frame load rejected a NaN coordinate");
  expect(f(g_badbox, 3.5, len), SEAMS_ERR_STATE, fn, "frame load rejected a wrong box");
  if (flags & P_CUTOFF) {
    expect(f(g_good, NAN, len), SEAMS_ERR_RANGE, fn, "NaN cutoff");
    expect(f(g_good, INFINITY, len), SEAMS_ERR_RANGE, fn, "infinite cutoff");
    expect(f(g_good, 0.0, len), SEAMS_ERR_RANGE, fn, "zero cutoff");
    expect(f(g_good, -3.5, len), SEAMS_ERR_RANGE, fn, "negative cutoff");
  }
  if (flags & P_LEN) {
    expect(f(g_good, 3.5, 0), SEAMS_ERR_RANGE, fn, "zero length");
    expect(f(g_good, 3.5, len - 1), SEAMS_ERR_RANGE, fn, "short length");
  }
  if (flags & P_NULLOUT) {
    g_null_out = 1;
    expect(f(g_good, 3.5, len), SEAMS_ERR_NULL, fn, "NULL output");
    g_null_out = 0;
  }
  if (flags & P_GRAPH) {
    g_graph = 99;
    expect(f(g_good, 3.5, len), SEAMS_ERR_ARG, fn, "unknown graph");
    g_graph = SEAMS_GRAPH_CUTOFF;
  }
  expect(f(g_good, 3.5, len), SEAMS_OK, fn, "valid call");
}

static void analysis_contract(void) {
  const unsigned C = P_CUTOFF, L = P_LEN, O = P_NULLOUT, G = P_GRAPH;
  for (int i = 0; i < g_n; ++i) {
    g_ibuf2[i] = i % 2;
  }
  probe("seams_neighbours_build", p_neighbours_build, C | O | G, 0);
  probe("seams_neighbours_knn", p_neighbours_knn, C | O, 0);
  probe("seams_chill", p_chill, C | L | O | G, 0);
  probe("seams_chill_plus", p_chill_plus, C | O | G, 0);
  probe("seams_chill_plus_n", p_chill_plus_n, C | L | O | G, 0);
  probe("seams_ql", p_ql, C | L | O | G, 0);
  probe("seams_q4q6", p_q4q6, C | O, 0);
  probe("seams_q4q6_n", p_q4q6_n, C | L | O, 0);
  probe("seams_q6bar", p_q6bar, C | L | O, 0);
  probe("seams_steinhardt", p_steinhardt, C | L | O, 0);
  probe("seams_system_symmetry", p_symmetry, L | O, 0);
  probe("seams_rings_find", p_rings_find, C | O | G, 0);
  probe("seams_cages_find", p_cages_find, C | O | G, 0);
  probe("seams_cages", p_cages, C | O | G, 0);
  probe("seams_cages_signature", p_cages_signature, C | O | G, 0);
  probe("seams_signature_cages", p_signature_cages, C | O | G, 0);
  probe("seams_brine_analyse", p_brine, L | O, 0);
  probe("seams_rdf", p_rdf, C | L | O, 20);
  probe("seams_running_cn", p_running_cn, C | L | O, 20);
  probe("seams_coordination", p_coordination, C | O, 0);
  probe("seams_hbonds", p_hbonds, C | O, 0);
  probe("seams_density_axis", p_density, L | O, 10);
  probe("seams_fingerprint", p_fingerprint, C | O | G, 0);
  probe("seams_contact_pairs", p_contacts, O, 0);
  probe("seams_domains", p_domains, C | O, 0);
  probe("seams_prisms", p_prisms, C | L | O, 7);
  probe("seams_monolayer", p_monolayer, C | L | O, 7);
  probe("seams_site_table", p_site_table, L | O, 0);

  /* Arguments the generic probe does not reach. */
  size_t a = 0, b = 0;
  expect(seams_cages_signature(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, NULL, &a, &b), SEAMS_ERR_NULL,
         "seams_cages_signature", "NULL spec");
  expect(seams_contact_pairs(g_good, NULL, &a, NULL, NULL), SEAMS_ERR_NULL, "seams_contact_pairs",
         "NULL spec");
  expect(seams_domains(g_good, 3.5, NULL, 0, NULL, &(int){0}, NULL), SEAMS_ERR_NULL,
         "seams_domains", "NULL spec");
  expect(seams_domains(g_good, 3.5, "1=waterO", 7, NULL, &(int){0}, NULL), SEAMS_ERR_ARG,
         "seams_domains", "unknown subset");
  expect(seams_site_table(g_good, NULL, g_ibuf, (size_t)g_n), SEAMS_ERR_NULL, "seams_site_table",
         "NULL spec");
  expect(seams_chill(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, 9, g_ibuf, (size_t)g_n, NULL),
         SEAMS_ERR_ARG, "seams_chill", "unknown rule");
  expect(seams_chill_plus_n(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, 9, g_ibuf, (size_t)g_n, NULL),
         SEAMS_ERR_ARG, "seams_chill_plus_n", "unknown rule");
  expect(seams_ql(g_good, 3.5, 0, 5, SEAMS_GRAPH_CUTOFF, g_dbuf, (size_t)g_n, NULL, 0),
         SEAMS_ERR_RANGE, "seams_ql", "unsupported order");
  expect(seams_q6bar(g_good, 3.5, 0, NAN, NULL, 0, &(double){0}, NULL, NULL), SEAMS_ERR_RANGE,
         "seams_q6bar", "NaN ice cut");
  expect(seams_rings_find(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, 2, &(seams_rings *){NULL}),
         SEAMS_ERR_RANGE, "seams_rings_find", "ring size below 3");
  expect(seams_rings_find(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, 13, &(seams_rings *){NULL}),
         SEAMS_ERR_RANGE, "seams_rings_find", "ring size above 12");
  expect(seams_neighbours_build(g_good, 3.5, 0, SEAMS_GRAPH_SEEDED, &(seams_neighbours *){NULL}),
         SEAMS_ERR_ARG, "seams_neighbours_build", "seeded graph is cages only");
  expect(seams_neighbours_knn(g_good, 0, 3.5, 0, 1, &(seams_neighbours *){NULL}), SEAMS_ERR_RANGE,
         "seams_neighbours_knn", "zero k");
  expect(seams_brine_analyse(g_good, NULL, (size_t)g_n, 0, 0, 4, 2, &(seams_brine *){NULL}),
         SEAMS_ERR_NULL, "seams_brine_analyse", "NULL ice flags");
  expect(seams_brine_analyse(g_good, g_ibuf2, (size_t)g_n + 1, 0, 0, 4, 2, &(seams_brine *){NULL}),
         SEAMS_ERR_RANGE, "seams_brine_analyse", "long ice flags");
  expect(seams_brine_analyse(g_good, g_ibuf2, (size_t)g_n, 0, 0, 4, 3, &(seams_brine *){NULL}),
         SEAMS_ERR_RANGE, "seams_brine_analyse", "axis 3");
  expect(seams_density_axis(g_good, 0, 10, -1, g_dbuf, g_dbuf3, 10, &a), SEAMS_ERR_RANGE,
         "seams_density_axis", "negative axis");
  expect(seams_hbonds(g_good, 3.5, 1, 0, 2.42, 30.0, &(int){0}), SEAMS_ERR_ARG, "seams_hbonds",
         "hydrogen type 0");
  expect(seams_hbonds(g_good, 3.5, 1, 2, NAN, 30.0, &(int){0}), SEAMS_ERR_RANGE, "seams_hbonds",
         "NaN distance cutoff");
  expect(seams_monolayer(g_good, 3.5, 0, 6, NAN, g_ibuf, 7), SEAMS_ERR_RANGE, "seams_monolayer",
         "NaN sheet area");
  expect(seams_monolayer(g_good, 3.5, 0, 6, 0.0, g_ibuf, 7), SEAMS_ERR_RANGE, "seams_monolayer",
         "zero sheet area");
  expect(seams_prisms(g_good, 3.5, 0, g_ibuf, NULL, 7), SEAMS_ERR_NULL, "seams_prisms",
         "NULL imperfect");
  expect(seams_q4q6_n(g_good, 3.5, 0, g_dbuf, (size_t)g_n, g_dbuf2, 1), SEAMS_ERR_RANGE,
         "seams_q4q6_n", "short q6 with a full q4");
  expect(seams_q4q6_n(g_good, 3.5, 0, NULL, 0, g_dbuf2, (size_t)g_n), SEAMS_OK, "seams_q4q6_n",
         "q6 alone");

  /* The length-checked spellings agree with the ABI 1 ones. */
  int a_labels[N_MAX], b_labels[N_MAX];
  expect(seams_chill_plus(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, a_labels, NULL),
         SEAMS_OK, "seams_chill_plus", "reference");
  expect(seams_chill_plus_n(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, b_labels,
                            (size_t)g_n, NULL),
         SEAMS_OK, "seams_chill_plus_n", "same frame");
  truth(memcmp(a_labels, b_labels, sizeof(int) * (size_t)g_n) == 0, "seams_chill_plus_n",
        "labels match seams_chill_plus");
}

/* ------------------------------------------------------------------ */
/* Result handles                                                      */

/* Cell edges far shorter or longer than the cutoff once overflowed the
 * cell count and stencil reach of the cutoff search (found by
 * capi_reload_fuzzer under UBSan). Any status but an internal error passes;
 * the sanitizer builds catch the overflow itself. */
static void extreme_cells_contract(void) {
  const double xyz[12] = {0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0, 3.1, 3.2};
  const double edges[2] = {1e-20, 1e250};
  seams_system *sys = seams_system_create();
  for (int e = 0; e < 2; ++e) {
    for (int tilted = 0; tilted < 2; ++tilted) {
      const double edge = edges[e];
      const double tilt = tilted ? 0.25 * edge : 0.0;
      const double box[9] = {0, 0, 0, edge, edge, edge, tilt, 0, 0};
      const char *what = e == 0 ? (tilted ? "thin tilted cell" : "thin cell")
                                : (tilted ? "long tilted cell" : "long cell");
      if (seams_system_set_atoms(sys, xyz, 4, box, 9, NULL, NULL) != SEAMS_OK) {
        continue;
      }
      int labels[4];
      truth(seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                               4, NULL) != SEAMS_ERR_INTERNAL,
            "seams_chill_plus_n", what);
      seams_neighbours *nb = NULL;
      truth(seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_ERR_INTERNAL,
            "seams_neighbours_build", what);
      seams_neighbours_destroy(nb);
    }
  }
  seams_system_destroy(sys);
}

static void neighbours_contract(void) {
  seams_neighbours *nb = NULL;
  expect(seams_neighbours_build(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb), SEAMS_OK,
         "seams_neighbours_build", "list");
  if (nb == NULL) {
    return;
  }
  const char *fn = "seams_neighbours_count";
  truth(seams_neighbours_count(NULL) == 0, fn, "NULL is zero");
  truth(seams_neighbours_count(nb) == g_n, fn, "one row per atom");

  fn = "seams_neighbours_nrows";
  size_t rows = 0;
  expect(seams_neighbours_nrows(NULL, &rows), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_neighbours_nrows(nb, NULL), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_neighbours_nrows(nb, &rows), SEAMS_OK, fn, "rows");
  truth(rows == (size_t)g_n, fn, "matches the count");

  fn = "seams_neighbours_copy";
  size_t n_off = 0, nnz = 0;
  expect(seams_neighbours_copy(NULL, NULL, 0, &n_off, NULL, 0, &nnz), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_neighbours_copy(nb, NULL, 0, NULL, NULL, 0, &nnz), SEAMS_ERR_NULL, fn, "NULL n_offsets");
  expect(seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, NULL), SEAMS_ERR_NULL, fn, "NULL nnz");
  expect(seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz), SEAMS_OK, fn, "query");
  truth(n_off == (size_t)g_n + 1 && nnz == (size_t)g_n * 4, fn, "ice Ic has four neighbours each");
  expect(seams_neighbours_copy(nb, g_ibuf, 0, &n_off, NULL, 0, &nnz), SEAMS_ERR_RANGE, fn,
         "zero offset capacity");
  expect(seams_neighbours_copy(nb, g_ibuf, n_off - 1, &n_off, NULL, 0, &nnz), SEAMS_ERR_RANGE, fn,
         "short offsets");
  int *cols = malloc(sizeof(int) * (nnz > 0 ? nnz : 1));
  expect(seams_neighbours_copy(nb, g_ibuf, n_off, &n_off, cols, 0, &nnz), SEAMS_ERR_RANGE, fn,
         "zero column capacity");
  expect(seams_neighbours_copy(nb, g_ibuf, n_off, &n_off, cols, nnz, &nnz), SEAMS_OK, fn, "copy");
  truth(g_ibuf[0] == 0 && (size_t)g_ibuf[g_n] == nnz, fn, "offsets are a prefix sum");

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
  fn = "seams_neighbours_csr";
  const int *off = NULL, *col = NULL;
  int ioff = 0, innz = 0;
  expect(seams_neighbours_csr(NULL, &off, &ioff, &col, &innz), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_neighbours_csr(nb, NULL, &ioff, &col, &innz), SEAMS_ERR_NULL, fn, "NULL offsets");
  expect(seams_neighbours_csr(nb, &off, &ioff, &col, NULL), SEAMS_ERR_NULL, fn, "NULL nnz");
  expect(seams_neighbours_csr(nb, &off, &ioff, &col, &innz), SEAMS_OK, fn, "borrow");
  truth(off != NULL && ioff == g_n + 1 && (size_t)innz == nnz &&
            memcmp(col, cols, sizeof(int) * nnz) == 0,
        fn, "borrowed view matches the copy");
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
  free(cols);

  fn = "seams_neighbours_destroy";
  seams_neighbours_destroy(NULL);
  cover(fn);
  seams_neighbours_destroy(nb);
  seams_neighbours_destroy(nb);
  truth(seams_last_error_message()[0] != '\0', fn, "double destroy is reported");
}

static void rings_contract(void) {
  seams_rings *rings = NULL;
  expect(seams_rings_find(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings), SEAMS_OK,
         "seams_rings_find", "rings");
  if (rings == NULL) {
    return;
  }
  const char *fn = "seams_rings_count";
  truth(seams_rings_count(NULL) == 0, fn, "NULL is zero");
  const int count = seams_rings_count(rings);
  truth(count > 0, fn, "ice Ic has six-rings");

  fn = "seams_rings_nrings";
  size_t nr = 0;
  expect(seams_rings_nrings(NULL, &nr), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_rings_nrings(rings, NULL), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_rings_nrings(rings, &nr), SEAMS_OK, fn, "count");
  truth(nr == (size_t)count, fn, "matches the int count");

  fn = "seams_rings_histogram";
  int hist[13];
  expect(seams_rings_histogram(NULL, hist, 13), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_rings_histogram(rings, NULL, 13), SEAMS_ERR_NULL, fn, "NULL counts");
  expect(seams_rings_histogram(rings, hist, 0), SEAMS_ERR_RANGE, fn, "zero length");
  expect(seams_rings_histogram(rings, hist, -1), SEAMS_ERR_RANGE, fn, "negative length");
  expect(seams_rings_histogram(rings, hist, 13), SEAMS_OK, fn, "histogram");
  truth(hist[6] == count, fn, "every ring is a six-ring");

  fn = "seams_rings_histogram_n";
  int hist_n[13];
  expect(seams_rings_histogram_n(NULL, hist_n, 13), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_rings_histogram_n(rings, NULL, 13), SEAMS_ERR_NULL, fn, "NULL counts");
  expect(seams_rings_histogram_n(rings, hist_n, 0), SEAMS_ERR_RANGE, fn, "zero length");
  expect(seams_rings_histogram_n(rings, hist_n, 13), SEAMS_OK, fn, "histogram");
  truth(memcmp(hist, hist_n, sizeof hist) == 0, fn, "matches the int-length histogram");
  expect(seams_rings_histogram_n(rings, hist_n, 4), SEAMS_OK, fn, "sizes past the buffer");
  truth(hist_n[0] == 0 && hist_n[3] == 0, fn, "sizes past the buffer are dropped");

  int *sizes = malloc(sizeof(int) * nr);
  int *sym = malloc(sizeof(int) * nr);
  fn = "seams_rings_labels";
  expect(seams_rings_labels(NULL, sizes, sym), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_rings_labels(rings, NULL, NULL), SEAMS_ERR_NULL, fn, "both outputs NULL");
  expect(seams_rings_labels(rings, sizes, sym), SEAMS_OK, fn, "labels");

  fn = "seams_rings_labels_n";
  int *sizes_n = malloc(sizeof(int) * nr);
  int *sym_n = malloc(sizeof(int) * nr);
  expect(seams_rings_labels_n(NULL, sizes_n, nr, sym_n, nr), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_rings_labels_n(rings, NULL, 0, NULL, 0), SEAMS_ERR_NULL, fn, "both outputs NULL");
  expect(seams_rings_labels_n(rings, sizes_n, 0, NULL, 0), SEAMS_ERR_RANGE, fn, "zero sizes length");
  expect(seams_rings_labels_n(rings, sizes_n, nr, sym_n, nr - 1), SEAMS_ERR_RANGE, fn,
         "short symmetry");
  expect(seams_rings_labels_n(rings, sizes_n, nr, sym_n, nr), SEAMS_OK, fn, "labels");
  truth(memcmp(sizes, sizes_n, sizeof(int) * nr) == 0 && memcmp(sym, sym_n, sizeof(int) * nr) == 0,
        fn, "matches seams_rings_labels");
  free(sizes);
  free(sym);
  free(sizes_n);
  free(sym_n);

  fn = "seams_rings_vertices";
  int atoms[12];
  size_t nv = 0;
  expect(seams_rings_vertices(NULL, 0, atoms, 12, &nv), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_rings_vertices(rings, 0, atoms, 12, NULL), SEAMS_ERR_NULL, fn, "NULL n_out");
  expect(seams_rings_vertices(rings, -1, atoms, 12, &nv), SEAMS_ERR_RANGE, fn, "negative index");
  expect(seams_rings_vertices(rings, count, atoms, 12, &nv), SEAMS_ERR_RANGE, fn, "index past the end");
  expect(seams_rings_vertices(rings, 0, NULL, 0, &nv), SEAMS_OK, fn, "query");
  truth(nv == 6, fn, "six vertices");
  expect(seams_rings_vertices(rings, 0, atoms, 0, &nv), SEAMS_ERR_RANGE, fn, "zero capacity");
  expect(seams_rings_vertices(rings, 0, atoms, 5, &nv), SEAMS_ERR_RANGE, fn, "short capacity");
  expect(seams_rings_vertices(rings, 0, atoms, 12, &nv), SEAMS_OK, fn, "vertices");

  fn = "seams_rings_fraction";
  double frac = -1;
  expect(seams_rings_fraction(NULL, 6, &frac), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_rings_fraction(rings, 6, NULL), SEAMS_ERR_NULL, fn, "NULL out");
  expect(seams_rings_fraction(rings, 6, &frac), SEAMS_OK, fn, "fraction");
  truth(frac > 0.99 && frac <= 1.0, fn, "every atom sits on a six-ring");

  fn = "seams_rings_destroy";
  seams_rings_destroy(NULL);
  cover(fn);
  seams_rings_destroy(rings);
  seams_rings_destroy(rings);
  truth(seams_last_error_message()[0] != '\0', fn, "double destroy is reported");
}

static void cages_contract(void) {
  seams_cage_set *cages = NULL;
  expect(seams_cages_find(g_good, 3.5, 0, SEAMS_GRAPH_CUTOFF, 0, 0, &cages), SEAMS_OK,
         "seams_cages_find", "cages");
  if (cages == NULL) {
    return;
  }
  const char *fn = "seams_cages_counts";
  int hex = -1, cubic = -1, water = -1, both = -1;
  expect(seams_cages_counts(NULL, &hex, &cubic, &water, &both), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_cages_counts(cages, NULL, NULL, NULL, NULL), SEAMS_ERR_NULL, fn, "every output NULL");
  expect(seams_cages_counts(cages, &hex, &cubic, &water, &both), SEAMS_OK, fn, "counts");
  truth(hex + cubic + water + both == g_n, fn, "counts partition the frame");

  fn = "seams_cages_labels";
  expect(seams_cages_labels(NULL, g_ibuf, (size_t)g_n), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_cages_labels(cages, NULL, (size_t)g_n), SEAMS_ERR_NULL, fn, "NULL buffer");
  expect(seams_cages_labels(cages, g_ibuf, 0), SEAMS_ERR_RANGE, fn, "zero length");
  expect(seams_cages_labels(cages, g_ibuf, (size_t)g_n - 1), SEAMS_ERR_RANGE, fn, "short length");
  expect(seams_cages_labels(cages, g_ibuf, (size_t)g_n), SEAMS_OK, fn, "labels");

  fn = "seams_cages_destroy";
  seams_cages_destroy(NULL);
  cover(fn);
  seams_cages_destroy(cages);
  seams_cages_destroy(cages);
  truth(seams_last_error_message()[0] != '\0', fn, "double destroy is reported");
}

static void brine_contract(void) {
  seams_brine *front = NULL;
  int ice[N_MAX];
  for (int i = 0; i < g_n; ++i) {
    ice[i] = g_xyz[3 * i + 2] < 0.5 * g_box[5];
  }
  expect(seams_brine_analyse(g_good, ice, (size_t)g_n, 0, 0, 4, 2, &front), SEAMS_OK,
         "seams_brine_analyse", "front");
  if (front == NULL) {
    return;
  }
  double x = 0;
  int flag = 0;
  const char *fn = "seams_brine_fronts";
  expect(seams_brine_fronts(NULL, &x, NULL, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_brine_fronts(front, NULL, NULL, NULL, NULL), SEAMS_ERR_NULL, fn, "every output NULL");
  expect(seams_brine_fronts(front, &x, &x, &x, &flag), SEAMS_OK, fn, "fronts");

  fn = "seams_brine_rejection";
  expect(seams_brine_rejection(NULL, &x, NULL, NULL, NULL), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_brine_rejection(front, NULL, NULL, NULL, NULL), SEAMS_ERR_NULL, fn, "every output NULL");
  expect(seams_brine_rejection(front, &x, &x, &x, &flag), SEAMS_OK, fn, "rejection");

  fn = "seams_brine_pockets";
  expect(seams_brine_pockets(NULL, &flag, NULL), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_brine_pockets(front, NULL, NULL), SEAMS_ERR_NULL, fn, "both outputs NULL");
  expect(seams_brine_pockets(front, &flag, &flag), SEAMS_OK, fn, "pockets");

  fn = "seams_brine_profile";
  size_t bins = 0;
  expect(seams_brine_profile(NULL, NULL, NULL, NULL, 0, &bins), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_brine_profile(front, NULL, NULL, NULL, 0, NULL), SEAMS_ERR_NULL, fn, "NULL n_out");
  expect(seams_brine_profile(front, NULL, NULL, NULL, 0, &bins), SEAMS_OK, fn, "query");
  truth(bins == 4, fn, "four bins");
  expect(seams_brine_profile(front, g_dbuf, NULL, NULL, 0, &bins), SEAMS_ERR_RANGE, fn, "zero capacity");
  expect(seams_brine_profile(front, g_dbuf, g_dbuf2, g_dbuf3, bins - 1, &bins), SEAMS_ERR_RANGE, fn,
         "short capacity");
  expect(seams_brine_profile(front, g_dbuf, g_dbuf2, g_dbuf3, bins, &bins), SEAMS_OK, fn, "profile");

  fn = "seams_brine_destroy";
  seams_brine_destroy(NULL);
  cover(fn);
  seams_brine_destroy(front);
  seams_brine_destroy(front);
  truth(seams_last_error_message()[0] != '\0', fn, "double destroy is reported");
}

static void fingerprint_contract(void) {
  seams_topo *fp = NULL;
  expect(seams_fingerprint(g_good, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 0, &fp), SEAMS_OK,
         "seams_fingerprint", "fingerprint");
  if (fp == NULL) {
    return;
  }
  char key[4096];
  size_t len = 0;
  const char *fn = "seams_fingerprint_key";
  expect(seams_fingerprint_key(NULL, key, sizeof key, &len), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_fingerprint_key(fp, key, sizeof key, NULL), SEAMS_ERR_NULL, fn, "NULL n_out");
  expect(seams_fingerprint_key(fp, NULL, 0, &len), SEAMS_OK, fn, "query");
  truth(len > 0 && len < sizeof key, fn, "non-empty key");
  expect(seams_fingerprint_key(fp, key, 0, &len), SEAMS_ERR_RANGE, fn, "zero capacity");
  expect(seams_fingerprint_key(fp, key, len, &len), SEAMS_ERR_RANGE, fn, "no room for the NUL");
  expect(seams_fingerprint_key(fp, key, len + 1, &len), SEAMS_OK, fn, "key");
  truth(strlen(key) == len, fn, "NUL-terminated");

  fn = "seams_fingerprint_atom";
  expect(seams_fingerprint_atom(NULL, 0, key, sizeof key, &len), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_fingerprint_atom(fp, 0, key, sizeof key, NULL), SEAMS_ERR_NULL, fn, "NULL n_out");
  expect(seams_fingerprint_atom(fp, -1, key, sizeof key, &len), SEAMS_ERR_RANGE, fn, "negative index");
  expect(seams_fingerprint_atom(fp, g_n, key, sizeof key, &len), SEAMS_ERR_RANGE, fn,
         "index past the end");
  expect(seams_fingerprint_atom(fp, 0, key, 0, &len), SEAMS_ERR_RANGE, fn, "zero capacity");
  expect(seams_fingerprint_atom(fp, 0, key, sizeof key, &len), SEAMS_OK, fn, "atom key");

  fn = "seams_fingerprint_rings";
  int census[13];
  expect(seams_fingerprint_rings(NULL, census, 13), SEAMS_ERR_NULL, fn, "NULL handle");
  expect(seams_fingerprint_rings(fp, NULL, 13), SEAMS_ERR_NULL, fn, "NULL counts");
  expect(seams_fingerprint_rings(fp, census, 0), SEAMS_ERR_RANGE, fn, "zero length");
  expect(seams_fingerprint_rings(fp, census, 13), SEAMS_OK, fn, "census");

  fn = "seams_fingerprint_destroy";
  seams_fingerprint_destroy(NULL);
  cover(fn);
  seams_fingerprint_destroy(fp);
  seams_fingerprint_destroy(fp);
  truth(seams_last_error_message()[0] != '\0', fn, "double destroy is reported");
}

/* ------------------------------------------------------------------ */

static int coverage(const char *list) {
  FILE *in = fopen(list, "r");
  if (in == NULL) {
    fprintf(stderr, "FAIL cannot open %s\n", list);
    return 1;
  }
  char line[256];
  int missing = 0;
  int listed = 0;
  while (fgets(line, sizeof line, in) != NULL) {
    line[strcspn(line, "\r\n")] = '\0';
    if (line[0] == '\0' || line[0] == '#') {
      continue;
    }
    listed += 1;
    int seen = 0;
    for (int i = 0; i < g_ncovered; ++i) {
      if (strcmp(g_covered[i], line) == 0) {
        seen = 1;
        break;
      }
    }
    if (!seen) {
      fprintf(stderr, "FAIL %s has no contract test\n", line);
      missing += 1;
    }
  }
  fclose(in);
  printf("capi_contract: %d listed functions, %d covered names\n", listed, g_ncovered);
  return missing;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s <symbols-list> <source-root>\n", argv[0]);
    return 2;
  }
  frames();
  meta_contract();
  frame_contract();
  frame_io_contract();
  read_contract(argv[2]);
  analysis_contract();
  extreme_cells_contract();
  neighbours_contract();
  rings_contract();
  cages_contract();
  brine_contract();
  fingerprint_contract();
  seams_system_destroy(g_good);
  seams_system_destroy(g_empty);
  seams_system_destroy(g_nan);
  seams_system_destroy(g_badbox);
  g_fails += coverage(argv[1]);
  printf("capi_contract: %d checks, %d failures\n", g_checks, g_fails);
  return g_fails == 0 ? 0 : 1;
}
