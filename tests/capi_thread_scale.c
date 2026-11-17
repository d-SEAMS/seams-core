/* seams_cutoff_csr and seams_ql_cutoff at 1, 2, 4 and 8 OpenMP threads on
 * periodic ice Ic. Every width must return the bytes the one-thread call
 * returned. The per-stage wall times from seams_c_last_stages go to a text
 * log; bench/capi_thread_scale.txt is one such log from a release build.
 *
 *   capi_thread_scale LOG [CELLS [REPEATS]]
 *
 * CELLS conventional diamond cells per edge (default 24, 110592 oxygens),
 * REPEATS timed calls per width after one untimed call (default 5). The
 * medians are logged. Exit 77 when the build has no OpenMP. */
#define _POSIX_C_SOURCE 199309L

#include <seams_c_api.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define N_WIDTHS 4
#define N_STAGES 4

static const int kWidths[N_WIDTHS] = {1, 2, 4, 8};
static const char *const kStageNames[N_STAGES] = {"cloud", "build", "ql", "copy"};
static const double kLattice = 6.37;
static const double kCutoff = 3.5;
static const int kOrder = 6;

struct timing {
  double total;
  double stage[N_STAGES];
};

static double now_seconds(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int by_value(const void *a, const void *b) {
  const double x = *(const double *)a;
  const double y = *(const double *)b;
  return (x > y) - (x < y);
}

static double median(double *values, int n) {
  qsort(values, (size_t)n, sizeof(double), by_value);
  return n % 2 ? values[n / 2] : 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

static void *must_alloc(size_t bytes) {
  void *p = malloc(bytes > 0 ? bytes : 1);
  if (p == NULL) {
    fprintf(stderr, "out of memory for %zu bytes\n", bytes);
    exit(1);
  }
  return p;
}

/* Oxygens of cubic ice: the diamond sites of an fcc lattice. */
static double *ice_ic(int cells, int *n_out) {
  static const double fcc[4][3] = {{0, 0, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0}};
  const int n = 8 * cells * cells * cells;
  double *xyz = must_alloc(sizeof(double) * 3 * (size_t)n);
  int at = 0;
  for (int i = 0; i < cells; i++) {
    for (int j = 0; j < cells; j++) {
      for (int k = 0; k < cells; k++) {
        for (int b = 0; b < 4; b++) {
          for (int shift = 0; shift < 2; shift++) {
            const double s = 0.25 * shift;
            xyz[3 * at] = kLattice * (i + fcc[b][0] + s);
            xyz[3 * at + 1] = kLattice * (j + fcc[b][1] + s);
            xyz[3 * at + 2] = kLattice * (k + fcc[b][2] + s);
            at++;
          }
        }
      }
    }
  }
  *n_out = n;
  return xyz;
}

static void set_width(int width) {
#ifdef _OPENMP
  omp_set_num_threads(width);
#else
  (void)width;
#endif
}

static void take_stages(struct timing *t, double total) {
  t->total = total;
  seams_c_last_stages(t->stage);
}

static void reduce(struct timing *runs, int repeats, struct timing *out) {
  double *column = must_alloc(sizeof(double) * (size_t)repeats);
  for (int r = 0; r < repeats; r++) {
    column[r] = runs[r].total;
  }
  out->total = median(column, repeats);
  for (int s = 0; s < N_STAGES; s++) {
    for (int r = 0; r < repeats; r++) {
      column[r] = runs[r].stage[s];
    }
    out->stage[s] = median(column, repeats);
  }
  free(column);
}

static void log_rows(FILE *out, const char *name, const struct timing rows[N_WIDTHS]) {
  fprintf(out, "%-16s %7s %9s %9s %9s %9s %9s %9s %8s\n", name, "threads", "total_ms",
          "cloud_ms", "build_ms", "ql_ms", "copy_ms", "other_ms", "speedup");
  for (int w = 0; w < N_WIDTHS; w++) {
    const struct timing *t = &rows[w];
    double inside = 0.0;
    for (int s = 0; s < N_STAGES; s++) {
      inside += t->stage[s];
    }
    fprintf(out, "%-16s %7d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %8.2f\n", name, kWidths[w],
            1e3 * t->total, 1e3 * t->stage[0], 1e3 * t->stage[1], 1e3 * t->stage[2],
            1e3 * t->stage[3], 1e3 * (t->total - inside), rows[0].total / t->total);
  }
}

static void log_split(FILE *out, const char *name, const struct timing *t) {
  fprintf(out, "%s at %d threads:", name, kWidths[N_WIDTHS - 1]);
  for (int s = 0; s < N_STAGES; s++) {
    fprintf(out, " %s %.1f%%", kStageNames[s], 100.0 * t->stage[s] / t->total);
  }
  fprintf(out, "\n");
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 4) {
    fprintf(stderr, "usage: %s LOG [CELLS [REPEATS]]\n", argv[0]);
    return 2;
  }
#ifndef _OPENMP
  printf("built without OpenMP; one width only\n");
  return 77;
#else
  const int cells = argc > 2 ? atoi(argv[2]) : 24;
  const int repeats = argc > 3 ? atoi(argv[3]) : 5;
  if (cells < 2 || repeats < 1) {
    fprintf(stderr, "CELLS must be at least 2 and REPEATS at least 1\n");
    return 2;
  }
  int n = 0;
  double *xyz = ice_ic(cells, &n);
  const double edge = kLattice * cells;
  const double box9[9] = {0, 0, 0, edge, edge, edge, 0, 0, 0};

  int *offsets[N_WIDTHS];
  int *cols[N_WIDTHS];
  double *ql[N_WIDTHS];
  int n_cols = 0;
  struct timing csr_rows[N_WIDTHS];
  struct timing ql_rows[N_WIDTHS];
  struct timing *runs = must_alloc(sizeof(struct timing) * (size_t)repeats);
  int failed = 0;

  for (int w = 0; w < N_WIDTHS; w++) {
    set_width(kWidths[w]);
    offsets[w] = must_alloc(sizeof(int) * ((size_t)n + 1));
    int count = 0;
    if (seams_cutoff_csr(xyz, n, box9, kCutoff, offsets[w], NULL, 0, &count) != 0) {
      fprintf(stderr, "seams_cutoff_csr count failed: %s\n", seams_c_error());
      return 1;
    }
    if (w == 0) {
      n_cols = count;
    } else if (count != n_cols) {
      fprintf(stderr, "%d threads count %d pairs, one thread %d\n", kWidths[w], count, n_cols);
      failed = 1;
    }
    cols[w] = must_alloc(sizeof(int) * (size_t)count);
    for (int r = 0; r < repeats; r++) {
      const double t0 = now_seconds();
      const int rc = seams_cutoff_csr(xyz, n, box9, kCutoff, offsets[w], cols[w], count, &count);
      const double t1 = now_seconds();
      if (rc != 0) {
        fprintf(stderr, "seams_cutoff_csr failed at %d threads\n", kWidths[w]);
        return 1;
      }
      take_stages(&runs[r], t1 - t0);
    }
    reduce(runs, repeats, &csr_rows[w]);

    ql[w] = must_alloc(sizeof(double) * (size_t)n);
    if (seams_ql_cutoff(xyz, n, box9, kCutoff, kOrder, ql[w]) != 0) {
      fprintf(stderr, "seams_ql_cutoff failed: %s\n", seams_c_error());
      return 1;
    }
    for (int r = 0; r < repeats; r++) {
      const double t0 = now_seconds();
      const int rc = seams_ql_cutoff(xyz, n, box9, kCutoff, kOrder, ql[w]);
      const double t1 = now_seconds();
      if (rc != 0) {
        fprintf(stderr, "seams_ql_cutoff failed at %d threads\n", kWidths[w]);
        return 1;
      }
      take_stages(&runs[r], t1 - t0);
    }
    reduce(runs, repeats, &ql_rows[w]);

    if (w > 0 && !failed) {
      if (memcmp(offsets[w], offsets[0], sizeof(int) * ((size_t)n + 1)) != 0 ||
          memcmp(cols[w], cols[0], sizeof(int) * (size_t)n_cols) != 0) {
        fprintf(stderr, "seams_cutoff_csr at %d threads differs from one thread\n", kWidths[w]);
        failed = 1;
      }
      if (memcmp(ql[w], ql[0], sizeof(double) * (size_t)n) != 0) {
        fprintf(stderr, "seams_ql_cutoff at %d threads differs from one thread\n", kWidths[w]);
        failed = 1;
      }
    }
  }

  /* Ideal ice Ic: four neighbours inside 3.5 A and one q6 for every site. */
  if (n_cols != 4 * n) {
    fprintf(stderr, "%d pairs, expected %d\n", n_cols, 4 * n);
    failed = 1;
  }
  for (int i = 0; i < n; i++) {
    if (!(fabs(ql[0][i] - ql[0][0]) < 1e-9) || !(ql[0][0] > 0.6)) {
      fprintf(stderr, "q6[%d] is %.12f, q6[0] is %.12f\n", i, ql[0][i], ql[0][0]);
      failed = 1;
      break;
    }
  }

  FILE *log = fopen(argv[1], "w");
  if (log == NULL) {
    perror(argv[1]);
    return 1;
  }
  FILE *sinks[2] = {stdout, log};
  for (int k = 0; k < 2; k++) {
    FILE *out = sinks[k];
    fprintf(out, "ice Ic %d^3 cells, %d oxygens, %d pairs, cutoff %.2f A, q%d\n", cells, n,
            n_cols, kCutoff, kOrder);
    fprintf(out, "omp_get_num_procs %d, median of %d calls after one untimed call\n",
            omp_get_num_procs(), repeats);
    fprintf(out, "outputs at 2, 4 and 8 threads %s the one-thread bytes\n",
            failed ? "DIFFER FROM" : "equal");
    log_rows(out, "seams_cutoff_csr", csr_rows);
    log_rows(out, "seams_ql_cutoff", ql_rows);
    log_split(out, "seams_cutoff_csr", &csr_rows[N_WIDTHS - 1]);
    log_split(out, "seams_ql_cutoff", &ql_rows[N_WIDTHS - 1]);
  }
  fclose(log);

  for (int w = 0; w < N_WIDTHS; w++) {
    free(offsets[w]);
    free(cols[w]);
    free(ql[w]);
  }
  free(runs);
  free(xyz);
  return failed;
#endif
}
