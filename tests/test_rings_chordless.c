/* Copyright (c) d-SEAMS developers */
/* seams_rings_chordless against a plain depth-first enumeration. */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fails = 0;

static void expect(int cond, const char *message) {
  if (!cond) {
    fprintf(stderr, "fail: %s (%s)\n", message, seams_last_error());
    g_fails += 1;
  }
}

typedef struct {
  const int *off;
  const int *cols;
  int size;
  int *path;
  char *on;
  int *out;
  size_t count;
  size_t cap;
} brute_t;

static int adjacent(const brute_t *b, int u, int v) {
  for (int a = b->off[u]; a < b->off[u + 1]; ++a) {
    if (b->cols[a] == v) {
      return 1;
    }
  }
  return 0;
}

/* Every simple path from s through higher vertices; a cycle is kept when
** no vertex touches a path vertex other than its two cycle neighbours. */
static void brute_extend(brute_t *b, int k) {
  const int s = b->path[0];
  const int last = b->path[k - 1];
  for (int a = b->off[last]; a < b->off[last + 1]; ++a) {
    const int v = b->cols[a];
    if (v <= s || b->on[v]) {
      continue;
    }
    int chord = 0;
    for (int i = 1; i + 1 < k && !chord; ++i) {
      chord = adjacent(b, v, b->path[i]);
    }
    if (chord) {
      continue;
    }
    const int closes = adjacent(b, v, s);
    if (k == b->size - 1) {
      if (closes && b->path[1] < v) {
        if (b->count < b->cap) {
          memcpy(b->out + (size_t)b->size * b->count, b->path, (size_t)k * sizeof(int));
          b->out[(size_t)b->size * b->count + k] = v;
        }
        b->count += 1;
      }
      continue;
    }
    if (closes && k > 1) {
      continue;
    }
    b->on[v] = 1;
    b->path[k] = v;
    brute_extend(b, k + 1);
    b->on[v] = 0;
  }
}

static int g_size = 0;

static int row_lex(const void *x, const void *y) {
  const int *a = x;
  const int *b = y;
  for (int i = 0; i < g_size; ++i) {
    if (a[i] != b[i]) {
      return a[i] < b[i] ? -1 : 1;
    }
  }
  return 0;
}

static void check_graph(seams_neighbours *nb, const char *tag) {
  size_t n_off = 0;
  size_t nnz = 0;
  expect(seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz) == SEAMS_OK, "csr size");
  int *off = malloc(n_off * sizeof(int));
  int *cols = malloc((nnz > 0 ? nnz : 1) * sizeof(int));
  expect(seams_neighbours_copy(nb, off, n_off, &n_off, cols, nnz, &nnz) == SEAMS_OK, "csr");
  const int n = (int)n_off - 1;
  int *per_atom = malloc((size_t)n * sizeof(int));
  int *ref_atom = malloc((size_t)n * sizeof(int));
  int path[12];
  char *on = calloc((size_t)n, 1);
  for (int size = 3; size <= 8; ++size) {
    size_t count = 0;
    expect(seams_rings_chordless(nb, size, NULL, 0, &count, per_atom, (size_t)n) == SEAMS_OK,
           "chordless query");
    int *rings = malloc((count > 0 ? count : 1) * (size_t)size * sizeof(int));
    size_t again = 0;
    expect(seams_rings_chordless(nb, size, rings, count * (size_t)size, &again, NULL, 0) ==
               SEAMS_OK,
           "chordless rings");
    int *ref = malloc((count + 1) * (size_t)size * sizeof(int));
    brute_t b = {off, cols, size, path, on, ref, 0, count + 1};
    for (int s = 0; s < n; ++s) {
      path[0] = s;
      on[s] = 1;
      brute_extend(&b, 1);
      on[s] = 0;
    }
    g_size = size;
    if (b.count <= b.cap) {
      qsort(ref, b.count, (size_t)size * sizeof(int), row_lex);
    }
    memset(ref_atom, 0, (size_t)n * sizeof(int));
    for (size_t r = 0; r < b.count && r < b.cap; ++r) {
      for (int k = 0; k < size; ++k) {
        ref_atom[ref[r * (size_t)size + (size_t)k]] += 1;
      }
    }
    char what[96];
    snprintf(what, sizeof what, "%s size %d: %zu cycles, brute force %zu", tag, size, count,
             b.count);
    printf("%s\n", what);
    expect(again == count && b.count == count, what);
    expect(count == 0 || memcmp(rings, ref, count * (size_t)size * sizeof(int)) == 0, what);
    expect(memcmp(per_atom, ref_atom, (size_t)n * sizeof(int)) == 0, what);

    if (count > 1) {
      size_t short_count = 0;
      int *rows = malloc((count - 1) * (size_t)size * sizeof(int));
      expect(seams_rings_chordless(nb, size, rows, (count - 1) * (size_t)size, &short_count,
                                   per_atom, (size_t)n) == SEAMS_ERR_RANGE,
             "short rings buffer");
      expect(short_count == count, "short buffer still counts");
      expect(memcmp(per_atom, ref_atom, (size_t)n * sizeof(int)) == 0,
             "short buffer still fills per_atom");
      free(rows);
    }
    free(rings);
    free(ref);
  }
  free(on);
  free(per_atom);
  free(ref_atom);
  free(off);
  free(cols);
}

int main(void) {
  const char *path = getenv("SEAMS_TRAJ");
  if (path == NULL) {
    fprintf(stderr, "SEAMS_TRAJ is not set\n");
    return 1;
  }
  seams_system *sys = seams_system_create();
  expect(seams_system_read(sys, path, 1, -1) == SEAMS_OK, "read");
  seams_neighbours *nb = NULL;
  expect(seams_neighbours_build(sys, 3.5, 2, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK, "graph");
  check_graph(nb, "oxygen cutoff");
  seams_neighbours_destroy(nb);
  nb = NULL;
  expect(seams_neighbours_build(sys, 3.2, 0, SEAMS_GRAPH_MUTUAL, &nb) == SEAMS_OK,
         "mutual graph");
  check_graph(nb, "all atoms mutual");

  size_t n_off = 0;
  size_t nnz = 0;
  seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz);
  const size_t n = n_off - 1;
  size_t one = 0;
  size_t four = 0;
  expect(seams_set_num_threads(1) == SEAMS_OK, "one thread");
  seams_rings_chordless(nb, 6, NULL, 0, &one, NULL, 0);
  int *a = malloc((one > 0 ? one : 1) * 6 * sizeof(int));
  int *b = malloc((one > 0 ? one : 1) * 6 * sizeof(int));
  expect(seams_rings_chordless(nb, 6, a, one * 6, &one, NULL, 0) == SEAMS_OK, "one thread rings");
  expect(seams_set_num_threads(4) == SEAMS_OK, "four threads");
  expect(seams_rings_chordless(nb, 6, b, one * 6, &four, NULL, 0) == SEAMS_OK,
         "four thread rings");
  expect(four == one && memcmp(a, b, one * 6 * sizeof(int)) == 0, "thread count independent");
  seams_set_num_threads(0);
  free(a);
  free(b);

  int *per = malloc(n * sizeof(int));
  size_t got = 0;
  expect(seams_rings_chordless(NULL, 6, NULL, 0, &got, NULL, 0) == SEAMS_ERR_NULL, "null graph");
  expect(seams_rings_chordless(nb, 6, NULL, 0, NULL, NULL, 0) == SEAMS_ERR_NULL, "no outputs");
  expect(seams_rings_chordless(nb, 6, per, 6, NULL, NULL, 0) == SEAMS_ERR_NULL,
         "rings without n_rings");
  expect(seams_rings_chordless(nb, 2, NULL, 0, &got, NULL, 0) == SEAMS_ERR_RANGE, "size 2");
  expect(seams_rings_chordless(nb, 13, NULL, 0, &got, NULL, 0) == SEAMS_ERR_RANGE, "size 13");
  expect(seams_rings_chordless(nb, 6, NULL, 0, NULL, per, n - 1) == SEAMS_ERR_RANGE,
         "short per_atom");
  free(per);
  seams_neighbours_destroy(nb);
  seams_system_destroy(sys);
  if (g_fails != 0) {
    fprintf(stderr, "%d checks failed\n", g_fails);
    return 1;
  }
  printf("seams_rings_chordless ok\n");
  return 0;
}
