/* Copyright (c) d-SEAMS developers
 * Steinhardt and CHILL+ calls on separate handles from many threads.
 *
 * Each of eight threads owns a handle holding its own jittered ice frame,
 * so every thread's q values differ. All threads start together and run
 * seams_q4q6_n, seams_ql (cutoff and mutual graphs, with the average),
 * seams_q6bar, seams_steinhardt and seams_chill_plus_n in a loop. Every result must match the
 * values that thread's frame gave on one thread, so scratch shared
 * between calls shows up as a mismatch, and under -Db_sanitize=thread as
 * a race.
 */
#define _POSIX_C_SOURCE 200809L

#include <seams.h>

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NTHREADS 8
#define ROUNDS 20
#define NATOMS 512

struct job {
  seams_system *sys;
  double q4[NATOMS], q6[NATOMS];
  double ql[NATOMS], qlbar[NATOMS];
  double mql[NATOMS], mqlbar[NATOMS];
  double q6bar[NATOMS];
  double s8[NATOMS], s8bar[NATOMS];
  double mean;
  int n_ice;
  int labels[NATOMS];
  long bad;
};

static struct job g_jobs[NTHREADS];
static pthread_barrier_t g_start;

static void ice_ic(int seed, double *xyz, double *box) {
  const int nrep = 4;
  const double a = 4.0 * 2.75 / sqrt(3.0);
  const double basis[8][3] = {{0, 0, 0},          {0.5, 0.5, 0},      {0.5, 0, 0.5},
                              {0, 0.5, 0.5},      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  unsigned state = 2463534242u + 977u * (unsigned)seed;
  int at = 0;
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (int b = 0; b < 8; ++b) {
          const double cell[3] = {i, j, k};
          for (int d = 0; d < 3; ++d) {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            const double jitter = 0.3 * ((double)(state % 2001u) / 1000.0 - 1.0);
            xyz[3 * at + d] = (cell[d] + basis[b][d]) * a + jitter;
          }
          at += 1;
        }
      }
    }
  }
  memset(box, 0, 9 * sizeof *box);
  box[3] = box[4] = box[5] = nrep * a;
}

static int run(struct job *out) {
  const size_t n = NATOMS;
  double median = 0;
  int ok = 1;
  ok &= seams_q4q6_n(out->sys, 3.5, 0, out->q4, n, out->q6, n) == SEAMS_OK;
  ok &= seams_ql(out->sys, 3.5, 0, 6, SEAMS_GRAPH_CUTOFF, out->ql, n, out->qlbar, n) == SEAMS_OK;
  ok &= seams_ql(out->sys, 3.5, 0, 4, SEAMS_GRAPH_MUTUAL, out->mql, n, out->mqlbar, n) == SEAMS_OK;
  ok &= seams_q6bar(out->sys, 3.5, 0, 0.4, out->q6bar, n, &out->mean, &median, &out->n_ice) ==
        SEAMS_OK;
  ok &= seams_steinhardt(out->sys, 3.5, 0, 8, out->s8, n, out->s8bar, n) == SEAMS_OK;
  ok &= seams_chill_plus_n(out->sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                           out->labels, n, NULL) == SEAMS_OK;
  return ok;
}

static long differ(const double *a, const double *b) {
  long bad = 0;
  for (int i = 0; i < NATOMS; ++i) {
    bad += !(fabs(a[i] - b[i]) <= 1e-12);
  }
  return bad;
}

static void *worker(void *arg) {
  struct job *ref = arg;
  struct job *mine = malloc(sizeof *mine);
  if (mine == NULL) {
    ref->bad += 1;
    return NULL;
  }
  mine->sys = ref->sys;
  pthread_barrier_wait(&g_start);
  for (int r = 0; r < ROUNDS; ++r) {
    if (!run(mine)) {
      fprintf(stderr, "call failed: %s\n", seams_last_error_message());
      ref->bad += 1;
      continue;
    }
    ref->bad += differ(mine->q4, ref->q4) + differ(mine->q6, ref->q6);
    ref->bad += differ(mine->ql, ref->ql) + differ(mine->qlbar, ref->qlbar);
    ref->bad += differ(mine->mql, ref->mql) + differ(mine->mqlbar, ref->mqlbar);
    ref->bad += differ(mine->q6bar, ref->q6bar);
    ref->bad += differ(mine->s8, ref->s8) + differ(mine->s8bar, ref->s8bar);
    ref->bad += memcmp(mine->labels, ref->labels, sizeof mine->labels) != 0;
    ref->bad += mine->n_ice != ref->n_ice || !(fabs(mine->mean - ref->mean) <= 1e-12);
  }
  free(mine);
  return NULL;
}

int main(void) {
  static double xyz[3 * NATOMS];
  double box[9];
  for (int t = 0; t < NTHREADS; ++t) {
    ice_ic(t, xyz, box);
    g_jobs[t].sys = seams_system_create();
    if (g_jobs[t].sys == NULL ||
        seams_system_set_atoms(g_jobs[t].sys, xyz, NATOMS, box, 9, NULL, NULL) != SEAMS_OK) {
      fprintf(stderr, "load failed: %s\n", seams_last_error_message());
      return 1;
    }
  }
  for (int t = 0; t < NTHREADS; ++t) {
    if (!run(&g_jobs[t])) {
      fprintf(stderr, "reference failed: %s\n", seams_last_error_message());
      return 1;
    }
  }
  if (memcmp(g_jobs[0].q6, g_jobs[1].q6, sizeof g_jobs[0].q6) == 0) {
    fprintf(stderr, "frames are not distinct\n");
    return 1;
  }

  pthread_t threads[NTHREADS];
  pthread_barrier_init(&g_start, NULL, NTHREADS);
  for (int t = 0; t < NTHREADS; ++t) {
    if (pthread_create(&threads[t], NULL, worker, &g_jobs[t]) != 0) {
      fprintf(stderr, "pthread_create failed\n");
      return 1;
    }
  }
  long total = 0;
  for (int t = 0; t < NTHREADS; ++t) {
    pthread_join(threads[t], NULL);
    total += g_jobs[t].bad;
    seams_system_destroy(g_jobs[t].sys);
  }
  pthread_barrier_destroy(&g_start);
  printf("capi_q_threads: %d threads x %d rounds, %d atoms each, %ld mismatches\n", NTHREADS,
         ROUNDS, NATOMS, total);
  return total == 0 ? 0 : 1;
}
