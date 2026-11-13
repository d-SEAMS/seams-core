/* Pure C11 proof of the stable C API. Compiled as C and linked to libseams. */
#define _POSIX_C_SOURCE 200809L

#include "seams.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int g_fails = 0;

static void fault(const char *message) {
  fprintf(stderr, "fail: %s (%s)\n", message, seams_last_error());
  g_fails += 1;
}

static void expect_status(seams_status got, seams_status want, const char *what) {
  if (got != want || seams_last_error()[0] == '\0') {
    fprintf(stderr, "fail: %s status %d (%s)\n", what, (int)got, seams_last_error());
    g_fails += 1;
  }
}

static int message_ok(const void *handle) {
  const char *text = seams_handle_error(handle);
  return text != NULL && text[0] != '\0';
}

typedef struct {
  int n;
  double *xyz;
  int *types;
  double box[9];
} Frame;

static void frame_clear(Frame *frame) {
  free(frame->xyz);
  free(frame->types);
  frame->xyz = NULL;
  frame->types = NULL;
  frame->n = 0;
}

static int frame_alloc(Frame *frame, int n) {
  frame->n = n;
  frame->xyz = calloc((size_t)n * 3, sizeof(double));
  frame->types = calloc((size_t)n, sizeof(int));
  if (frame->xyz == NULL || frame->types == NULL) {
    return -1;
  }
  for (int i = 0; i < n; ++i) {
    frame->types[i] = 1;
  }
  return 0;
}

static void ice_ic(Frame *frame) {
  const int nrep = 2;
  const double bond = 2.75;
  const double a = 4.0 * bond / sqrt(3.0);
  const double basis[8][3] = {{0, 0, 0},          {0.5, 0.5, 0},      {0.5, 0, 0.5},
                              {0, 0.5, 0.5},      {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25},
                              {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
  frame_alloc(frame, nrep * nrep * nrep * 8);
  int at = 0;
  for (int i = 0; i < nrep; ++i) {
    for (int j = 0; j < nrep; ++j) {
      for (int k = 0; k < nrep; ++k) {
        for (int b = 0; b < 8; ++b) {
          frame->xyz[3 * at] = (i + basis[b][0]) * a;
          frame->xyz[3 * at + 1] = (j + basis[b][1]) * a;
          frame->xyz[3 * at + 2] = (k + basis[b][2]) * a;
          at += 1;
        }
      }
    }
  }
  frame->box[3] = frame->box[4] = frame->box[5] = nrep * a;
}

static int load_dump(const char *path, Frame *frame) {
  FILE *in = fopen(path, "r");
  if (in == NULL) {
    return -1;
  }
  char line[1024];
  int n = -1;
  int got = 0;
  double lo[3] = {0}, hi[3] = {0}, tilt[3] = {0};
  Frame tmp;
  memset(&tmp, 0, sizeof tmp);
  while (fgets(line, sizeof line, in) != NULL) {
    if (strncmp(line, "ITEM: NUMBER OF ATOMS", 21) == 0) {
      if (fgets(line, sizeof line, in) == NULL) {
        break;
      }
      n = atoi(line);
    } else if (strncmp(line, "ITEM: BOX BOUNDS", 16) == 0) {
      for (int r = 0; r < 3; ++r) {
        if (fgets(line, sizeof line, in) == NULL) {
          break;
        }
        double a = 0, b = 0, c = 0;
        int nread = sscanf(line, "%lf %lf %lf", &a, &b, &c);
        lo[r] = a;
        hi[r] = nread >= 2 ? b : a;
        tilt[r] = nread >= 3 ? c : 0;
      }
    } else if (strncmp(line, "ITEM: ATOMS", 11) == 0) {
      if (n < 1 || frame_alloc(&tmp, n) != 0) {
        fclose(in);
        return -1;
      }
      while (got < n && fgets(line, sizeof line, in) != NULL) {
        int id = 0, type = 1;
        double x = 0, y = 0, z = 0;
        if (sscanf(line, "%d %d %lf %lf %lf", &id, &type, &x, &y, &z) != 5) {
          continue;
        }
        tmp.xyz[3 * got] = x;
        tmp.xyz[3 * got + 1] = y;
        tmp.xyz[3 * got + 2] = z;
        tmp.types[got] = type;
        got += 1;
      }
      break;
    }
  }
  fclose(in);
  if (got != n || n < 1) {
    frame_clear(&tmp);
    return -1;
  }
  tmp.box[0] = lo[0];
  tmp.box[1] = lo[1];
  tmp.box[2] = lo[2];
  tmp.box[3] = hi[0] - lo[0];
  tmp.box[4] = hi[1] - lo[1];
  tmp.box[5] = hi[2] - lo[2];
  tmp.box[6] = tilt[0];
  tmp.box[7] = tilt[1];
  tmp.box[8] = tilt[2];
  *frame = tmp;
  return 0;
}

static int load_oxygen(const char *path, Frame *frame) {
  FILE *in = fopen(path, "r");
  if (in == NULL) {
    return -1;
  }
  int n = 0;
  double lx = 0, ly = 0, lz = 0;
  if (fscanf(in, "%d %lf %lf %lf", &n, &lx, &ly, &lz) != 4 || n < 1 ||
      frame_alloc(frame, n) != 0) {
    fclose(in);
    return -1;
  }
  for (int i = 0; i < n; ++i) {
    if (fscanf(in, "%lf %lf %lf", &frame->xyz[3 * i], &frame->xyz[3 * i + 1],
               &frame->xyz[3 * i + 2]) != 3) {
      frame_clear(frame);
      fclose(in);
      return -1;
    }
  }
  fclose(in);
  frame->box[3] = lx;
  frame->box[4] = ly;
  frame->box[5] = lz;
  return 0;
}

static int write_frame(const char *path, const Frame *frame) {
  FILE *out = fopen(path, "w");
  if (out == NULL) {
    return -1;
  }
  fprintf(out, "%d\n", frame->n);
  for (int i = 0; i < 9; ++i) {
    fprintf(out, "%.17g%c", frame->box[i], i == 8 ? '\n' : ' ');
  }
  for (int i = 0; i < frame->n; ++i) {
    fprintf(out, "%.17g %.17g %.17g %d\n", frame->xyz[3 * i], frame->xyz[3 * i + 1],
            frame->xyz[3 * i + 2], frame->types[i]);
  }
  fclose(out);
  return 0;
}

static int store(seams_system *sys, const Frame *frame) {
  return seams_system_set_frame(sys, frame->xyz, frame->n, frame->box, frame->types, NULL) ==
         SEAMS_OK;
}

static void names_and_versions(void) {
  if (seams_abi() != SEAMS_ABI) {
    fault("abi");
  }
  if (seams_version() != SEAMS_VERSION_NUMBER) {
    fault("packed version");
  }
  if (seams_version_major() != SEAMS_VERSION_MAJOR ||
      seams_version_minor() != SEAMS_VERSION_MINOR ||
      seams_version_patch() != SEAMS_VERSION_PATCH) {
    fault("version parts");
  }
  if (strcmp(seams_version_string(), SEAMS_VERSION_STRING) != 0) {
    fault("version string");
  }
  for (int code = 0; code <= 7; ++code) {
    if (seams_status_string((seams_status)code)[0] == '\0') {
      fault("status string");
    }
  }
  for (int label = 0; label <= 8; ++label) {
    if (seams_label_name(label)[0] == '\0') {
      fault("label name");
    }
  }
  for (int sym = 0; sym <= 11; ++sym) {
    if (seams_ring_symmetry_name(sym)[0] == '\0') {
      fault("ring name");
    }
  }
  if (seams_last_error() == NULL || seams_handle_error(NULL) == NULL) {
    fault("error pointers");
  }
}

static void error_paths(const Frame *frame) {
  double bad[9] = {0, 0, 0, 0, 10, 10, 0, 0, 0};
  double nanbox[9];
  memcpy(nanbox, frame->box, sizeof nanbox);
  nanbox[4] = nan("");
  expect_status(seams_system_set_frame(NULL, frame->xyz, frame->n, frame->box, NULL, NULL),
                SEAMS_ERR_NULL, "null system");
  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    fault("create");
    return;
  }
  expect_status(seams_system_set_atoms(sys, frame->xyz, 0, frame->box, 9, NULL, NULL),
                SEAMS_ERR_RANGE, "zero atoms");
  expect_status(seams_system_set_atoms(sys, NULL, (size_t)frame->n, frame->box, 9, NULL, NULL),
                SEAMS_ERR_NULL, "null coordinates");
  expect_status(seams_system_set_atoms(sys, frame->xyz, (size_t)frame->n, frame->box, 3, NULL, NULL),
                SEAMS_ERR_RANGE, "short box buffer");
  expect_status(seams_system_set_frame(sys, frame->xyz, frame->n, bad, NULL, NULL), SEAMS_ERR_BOX,
                "flat box");
  if (!message_ok(sys)) {
    fault("handle error after a bad box");
  }
  expect_status(seams_system_set_frame(sys, frame->xyz, frame->n, nanbox, NULL, NULL),
                SEAMS_ERR_BOX, "non-finite box");
  expect_status(seams_system_set_ortho(sys, frame->xyz, 0, 0, 0, 0, 1, 1, 1, NULL, NULL),
                SEAMS_ERR_RANGE, "zero-length ortho");
  expect_status(seams_system_read(NULL, "missing", 1, 0), SEAMS_ERR_NULL, "null read");
  expect_status(seams_system_read(sys, NULL, 1, 0), SEAMS_ERR_NULL, "null path");
  expect_status(seams_system_read(sys, "input/traj/ice1h.lammpstrj", 0, 0), SEAMS_ERR_RANGE,
                "frame index");
  size_t nframes = 99;
  expect_status(seams_path_nframes(NULL, &nframes), SEAMS_ERR_NULL, "null nframes");
  expect_status(seams_path_nframes("input/traj/ice1h.lammpstrj", NULL), SEAMS_ERR_NULL,
                "null nframes out");
  seams_neighbours *nb = NULL;
  expect_status(seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, NULL), SEAMS_ERR_NULL,
                "null neighbour out");
  expect_status(seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb), SEAMS_ERR_STATE,
                "neighbours before a frame");
  expect_status(seams_set_num_threads(-1), SEAMS_ERR_RANGE, "negative threads");
  if (seams_system_set_frame(sys, frame->xyz, frame->n, frame->box, NULL, NULL) != SEAMS_OK) {
    fault("frame for short buffers");
  } else {
    int one = 0;
    expect_status(seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, &one, 0,
                              NULL),
                  SEAMS_ERR_RANGE, "zero-length chill buffer");
    double q = 0;
    expect_status(seams_ql(sys, 3.5, 0, 3, SEAMS_GRAPH_CUTOFF, &q, 0, NULL, 0), SEAMS_ERR_RANGE,
                  "zero-length q buffer");
    expect_status(seams_ql(sys, 3.5, 0, 5, SEAMS_GRAPH_CUTOFF, &q, 1, NULL, 0), SEAMS_ERR_RANGE,
                  "bad order");
    if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK) {
      size_t n_off = 0, nnz = 0;
      int tiny = 0;
      expect_status(seams_neighbours_copy(nb, &tiny, 1, &n_off, NULL, 0, &nnz), SEAMS_ERR_RANGE,
                    "short neighbour buffer");
      seams_neighbours_destroy(nb);
      seams_neighbours_destroy(nb);
      if (seams_last_error()[0] == '\0') {
        fault("double destroy of neighbours");
      }
    }
    char key[4];
    size_t n_out = 0;
    seams_topo *fp = NULL;
    if (seams_fingerprint(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 1, 0, &fp) == SEAMS_OK) {
      expect_status(seams_fingerprint_key(fp, key, 1, &n_out), SEAMS_ERR_RANGE,
                    "short fingerprint buffer");
      seams_fingerprint_destroy(fp);
      seams_fingerprint_destroy(fp);
      if (seams_last_error()[0] == '\0') {
        fault("double destroy of fingerprint");
      }
    }
  }
  seams_system_destroy(sys);
  seams_system_destroy(sys);
  if (seams_last_error()[0] == '\0') {
    fault("double destroy of system");
  }
  seams_system_destroy(NULL);
  seams_neighbours_destroy(NULL);
  seams_rings_destroy(NULL);
  seams_cages_destroy(NULL);
  seams_brine_destroy(NULL);
  seams_fingerprint_destroy(NULL);
}

static double fingerprint(const double *values, int n) {
  long double sum = 0;
  for (int i = 0; i < n; ++i) {
    double x = isfinite(values[i]) ? values[i] : 0;
    sum += (long double)x * (i + 1);
  }
  return (double)sum;
}

static double mean_of(const double *values, int n) {
  long double sum = 0;
  int count = 0;
  for (int i = 0; i < n; ++i) {
    if (!isfinite(values[i])) {
      continue;
    }
    sum += values[i];
    count += 1;
  }
  return count > 0 ? (double)(sum / count) : 0;
}

typedef struct {
  char key[32];
  int is_int;
  int ival;
  double dval;
} Rec;

static Rec g_rec[160];
static int g_nrec;

static void put_i(const char *key, int value) {
  Rec *rec = &g_rec[g_nrec++];
  snprintf(rec->key, sizeof rec->key, "%s", key);
  rec->is_int = 1;
  rec->ival = value;
}

static void put_d(const char *key, double value) {
  Rec *rec = &g_rec[g_nrec++];
  snprintf(rec->key, sizeof rec->key, "%s", key);
  rec->is_int = 0;
  rec->dval = value;
}

static void hist9(const char *prefix, const int *labels, int n) {
  int counts[9] = {0};
  for (int i = 0; i < n; ++i) {
    int label = labels[i];
    if (label < 0 || label > 8) {
      label = 6;
    }
    counts[label] += 1;
  }
  for (int i = 0; i < 9; ++i) {
    char key[32];
    snprintf(key, sizeof key, "%s_%d", prefix, i);
    put_i(key, counts[i]);
  }
}

static const Rec *find_rec(const char *key) {
  for (int i = 0; i < g_nrec; ++i) {
    if (strcmp(g_rec[i].key, key) == 0) {
      return &g_rec[i];
    }
  }
  return NULL;
}

static void call_every(seams_system *sys, const Frame *frame, const char *kind) {
  const int n = frame->n;
  g_nrec = 0;
  put_i("nop", seams_system_count(sys));
  seams_neighbours *nb = NULL;
  if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_OK) {
    fault("neighbours");
  } else {
    const int *offsets = NULL;
    const int *cols = NULL;
    int n_offsets = 0, nnz = 0;
    seams_status csr_rc;
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    csr_rc = seams_neighbours_csr(nb, &offsets, &n_offsets, &cols, &nnz);
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    if (csr_rc != SEAMS_OK || seams_neighbours_count(nb) != n) {
      fault("csr");
    } else {
      int lo = 1000, hi = 0;
      for (int i = 0; i < n; ++i) {
        int degree = offsets[i + 1] - offsets[i];
        if (degree < lo) {
          lo = degree;
        }
        if (degree > hi) {
          hi = degree;
        }
      }
      put_i("deg_min", lo);
      put_i("deg_max", hi);
      put_i("nnz", nnz);
      size_t need_off = 0, need_col = 0;
      if (seams_neighbours_copy(nb, NULL, 0, &need_off, NULL, 0, &need_col) != SEAMS_OK) {
        fault("csr query");
      }
      put_i("csr_nnz", (int)need_col);
    }
    seams_neighbours_destroy(nb);
  }
  seams_neighbours *knn = NULL;
  if (seams_neighbours_knn(sys, 4, 3.5, 0, 1, &knn) != SEAMS_OK) {
    fault("knn");
  } else {
    size_t need_off = 0, need_col = 0;
    if (seams_neighbours_copy(knn, NULL, 0, &need_off, NULL, 0, &need_col) != SEAMS_OK) {
      fault("knn copy");
    }
    put_i("knn_nnz", (int)need_col);
    seams_neighbours_destroy(knn);
  }
  seams_rings *rings = NULL;
  if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    fault("rings");
  } else {
    int hist[7] = {0};
    if (seams_rings_histogram(rings, hist, 7) != SEAMS_OK) {
      fault("histogram");
    }
    for (int i = 0; i < 7; ++i) {
      char key[32];
      snprintf(key, sizeof key, "ring_%d", i);
      put_i(key, hist[i]);
    }
    int nr = seams_rings_count(rings);
    if (nr > 0) {
      int *sizes = calloc((size_t)nr, sizeof(int));
      int *sym = calloc((size_t)nr, sizeof(int));
      if (seams_rings_labels(rings, sizes, sym) != SEAMS_OK) {
        fault("ring labels");
      }
      size_t n_vert = 0;
      if (seams_rings_vertices(rings, 0, NULL, 0, &n_vert) != SEAMS_OK) {
        fault("vertex query");
      } else if (n_vert > 0) {
        int *atoms = calloc(n_vert, sizeof(int));
        if (seams_rings_vertices(rings, 0, atoms, n_vert, &n_vert) != SEAMS_OK) {
          fault("vertices");
        }
        free(atoms);
      }
      free(sizes);
      free(sym);
    }
    double fraction = 0;
    if (seams_rings_fraction(rings, 6, &fraction) != SEAMS_OK) {
      fault("fraction");
    }
    seams_rings_destroy(rings);
  }
  for (int order = 3; order <= 12; order += order == 3 ? 1 : order == 4 ? 2 : order == 6 ? 2 : 4) {
    double *ql = calloc((size_t)n, sizeof(double));
    double *bar = calloc((size_t)n, sizeof(double));
    seams_status rc =
        seams_ql(sys, 3.5, 0, order, SEAMS_GRAPH_CUTOFF, ql, (size_t)n, bar, (size_t)n);
    if (rc != SEAMS_OK) {
      fault("ql");
    } else if (order == 3 || order == 4 || order == 6 || order == 8 || order == 12) {
      char key[32];
      snprintf(key, sizeof key, "q%d_mean", order);
      put_d(key, mean_of(ql, n));
      snprintf(key, sizeof key, "q%d_fp", order);
      put_d(key, fingerprint(ql, n));
    }
    if (seams_steinhardt(sys, 3.5, 0, order, ql, (size_t)n, NULL, 0) != SEAMS_OK) {
      fault("steinhardt");
    }
    free(ql);
    free(bar);
  }
  double *q4 = calloc((size_t)n, sizeof(double));
  double *q6 = calloc((size_t)n, sizeof(double));
  double *q6bar = calloc((size_t)n, sizeof(double));
  if (seams_q4q6(sys, 3.5, 0, q4, q6) != SEAMS_OK) {
    fault("q4q6");
  }
  double mean = 0, median = 0;
  int n_ice = 0;
  if (seams_q6bar(sys, 3.5, 0, 0.36, q6bar, (size_t)n, &mean, &median, &n_ice) != SEAMS_OK) {
    fault("q6bar");
  } else {
    put_d("q6bar_mean", mean);
    put_d("q6bar_p50", median);
    put_d("q6bar_fp", fingerprint(q6bar, n));
    put_i("q6bar_ice", n_ice);
    if (strcmp(kind, "nucleation") == 0) {
      put_d("nuc_mean", mean);
      put_d("nuc_p50", median);
      put_i("nuc_ice", n_ice);
    }
  }
  int *labels = calloc((size_t)n, sizeof(int));
  int counts[9];
  if (seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, (size_t)n,
                  counts) != SEAMS_OK) {
    fault("chill");
  } else {
    hist9("chill", labels, n);
  }
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts) !=
      SEAMS_OK) {
    fault("chill+");
  } else {
    hist9("plus", labels, n);
  }
  if (seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_OVITO, labels, counts) !=
      SEAMS_OK) {
    fault("ovito chill+");
  } else {
    hist9("ovito", labels, n);
  }
  if (seams_system_symmetry(sys, 0.12, labels, (size_t)n) != SEAMS_OK) {
    fault("symmetry");
  } else {
    int ih = 0, ic = 0, isd = 0, other = 0;
    for (int i = 0; i < n; ++i) {
      if (labels[i] == 0) {
        ih += 1;
      } else if (labels[i] == 1) {
        ic += 1;
      } else if (labels[i] == 2) {
        isd += 1;
      } else {
        other += 1;
      }
    }
    put_i("sym_ih", ih);
    put_i("sym_ic", ic);
    put_i("sym_isd", isd);
    put_i("sym_other", other);
  }
  seams_cage_set *cages = NULL;
  if (seams_cages_find(sys, 3.5, 0, SEAMS_GRAPH_SEEDED, 4, 0, &cages) != SEAMS_OK) {
    fault("cages");
  } else {
    int n_hex = 0, n_cubic = 0, n_water = 0, n_both = 0;
    if (seams_cages_counts(cages, &n_hex, &n_cubic, &n_water, &n_both) != SEAMS_OK) {
      fault("cage counts");
    }
    put_i("cage_hex", n_hex);
    put_i("cage_cubic", n_cubic);
    put_i("cage_water", n_water);
    put_i("cage_both", n_both);
    if (seams_cages_labels(cages, labels, (size_t)n) != SEAMS_OK) {
      fault("cage labels");
    }
    seams_cages_destroy(cages);
    seams_cages_destroy(cages);
  }
  seams_cage_set *again = NULL;
  if (seams_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 4, 0, &again) == SEAMS_OK) {
    seams_cages_destroy(again);
  } else {
    fault("cages alias");
  }
  size_t n_cages = 0, n_atoms = 0;
  if (seams_cages_signature(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", &n_cages, &n_atoms) !=
      SEAMS_OK) {
    fault("signature");
  } else {
    put_i("sig_cages", (int)n_cages);
    put_i("sig_atoms", (int)n_atoms);
  }
  if (seams_signature_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", &n_cages, &n_atoms) !=
      SEAMS_OK) {
    fault("signature alias");
  }
  size_t n_out = 0;
  double volume = 0;
  if (seams_rdf(sys, 10, 20, 0, 0, NULL, NULL, NULL, 0, &n_out, &volume) != SEAMS_OK ||
      n_out == 0) {
    fault("rdf query");
  } else {
    double *r = calloc(n_out, sizeof(double));
    double *g = calloc(n_out, sizeof(double));
    int *count = calloc(n_out, sizeof(int));
    double *cn = calloc(n_out, sizeof(double));
    if (seams_rdf(sys, 10, 20, 0, 0, r, g, count, n_out, &n_out, &volume) != SEAMS_OK) {
      fault("rdf");
    } else {
      double gsum = 0;
      int pairs = 0;
      for (size_t i = 0; i < n_out; ++i) {
        gsum += g[i];
        pairs += count[i];
      }
      put_d("rdf_vol", volume);
      put_d("rdf_gsum", gsum);
      put_i("rdf_count", pairs);
    }
    if (seams_running_cn(sys, 10, 20, 0, 0, r, cn, n_out, &n_out) != SEAMS_OK) {
      fault("running cn");
    } else if (n_out > 0) {
      put_d("cn_last", cn[n_out - 1]);
    }
    double coordination = 0;
    if (seams_coordination(sys, 3.5, 20, 0, 0, &coordination) != SEAMS_OK) {
      fault("coordination");
    } else {
      put_d("cn", coordination);
    }
    free(r);
    free(g);
    free(count);
    free(cn);
  }
  if (seams_density_axis(sys, 0, 10, 2, NULL, NULL, 0, &n_out) != SEAMS_OK) {
    fault("density query");
  } else {
    double *z = calloc(n_out, sizeof(double));
    double *rho = calloc(n_out, sizeof(double));
    if (seams_density_axis(sys, 0, 10, 2, z, rho, n_out, &n_out) != SEAMS_OK) {
      fault("density");
    } else {
      double sum = 0;
      for (size_t i = 0; i < n_out; ++i) {
        sum += rho[i];
      }
      put_d("dens_sum", sum);
    }
    free(z);
    free(rho);
  }
  int n_bonds = -1;
  if (seams_hbonds(sys, 3.5, 1, 2, 2.5, 30, &n_bonds) != SEAMS_OK) {
    fault("hbonds");
  } else {
    put_i("hbond_n", n_bonds);
  }
  int perfect[7] = {0}, imperfect[7] = {0};
  if (seams_prisms(sys, 3.5, 0, perfect, imperfect, 7) != SEAMS_OK) {
    fault("prisms");
  } else {
    int psum = 0, isum = 0;
    for (int i = 3; i <= 6; ++i) {
      psum += perfect[i];
      isum += imperfect[i];
    }
    put_i("prism_p", psum);
    put_i("prism_i", isum);
  }
  int mono[7] = {0};
  if (seams_monolayer(sys, 3.5, 0, 6, frame->box[3] * frame->box[4], mono, 7) != SEAMS_OK) {
    fault("monolayer");
  } else {
    int total = 0;
    for (int i = 0; i < 7; ++i) {
      total += mono[i];
    }
    put_i("mono_n", total);
  }
  if (seams_site_table(sys, "1=waterO", labels, (size_t)n) != SEAMS_OK) {
    fault("site table");
  } else {
    int water = 0;
    for (int i = 0; i < n; ++i) {
      if (labels[i] == 8) {
        water += 1;
      }
    }
    put_i("site_water", water);
  }
  size_t n_pairs = 0;
  int n_cation = 0, n_anion = 0;
  if (seams_contact_pairs(sys, "1=cationHead,2=anion", &n_pairs, &n_cation, &n_anion) !=
      SEAMS_OK) {
    fault("contacts");
  } else {
    put_i("contact_n", (int)n_pairs);
    put_i("contact_cat", n_cation);
    put_i("contact_an", n_anion);
  }
  int n_subset = 0, largest = 0;
  double percolation = 0;
  if (seams_domains(sys, 3.5, "1=apolar", 0, &n_subset, &largest, &percolation) != SEAMS_OK) {
    fault("domains");
  } else {
    put_i("domain_n", n_subset);
    put_i("domain_largest", largest);
    put_d("domain_p", percolation);
  }
  seams_topo *fp = NULL;
  if (seams_fingerprint(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 1, &fp) != SEAMS_OK) {
    fault("fingerprint");
  } else {
    size_t key_len = 0;
    size_t frame_key = 0;
    if (seams_fingerprint_key(fp, NULL, 0, &key_len) != SEAMS_OK) {
      fault("fingerprint query");
    } else {
      frame_key = key_len;
      char *key = calloc(key_len + 1, 1);
      if (seams_fingerprint_key(fp, key, key_len + 1, &key_len) != SEAMS_OK) {
        fault("fingerprint key");
      }
      free(key);
    }
    if (seams_fingerprint_atom(fp, 0, NULL, 0, &key_len) == SEAMS_OK && key_len > 0) {
      char *atom = calloc(key_len + 1, 1);
      (void)seams_fingerprint_atom(fp, 0, atom, key_len + 1, &key_len);
      free(atom);
    }
    int census[8] = {0};
    if (seams_fingerprint_rings(fp, census, 8) != SEAMS_OK) {
      fault("fingerprint rings");
    } else {
      int total = 0;
      for (int i = 0; i < 8; ++i) {
        total += census[i];
      }
      put_i("fp_len", (int)frame_key);
      put_i("fp_rings", total);
    }
    seams_fingerprint_destroy(fp);
  }
  if (strcmp(kind, "brine") == 0) {
    int *ice = calloc((size_t)n, sizeof(int));
    if (seams_chill_plus(sys, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                         counts) == SEAMS_OK) {
      for (int i = 0; i < n; ++i) {
        ice[i] = labels[i] == 0 || labels[i] == 1 || labels[i] == 7 || labels[i] == 8;
      }
      seams_brine *front = NULL;
      if (seams_brine_analyse(sys, ice, (size_t)n, 1, 3, 20, 2, &front) != SEAMS_OK) {
        fault("brine");
      } else {
        double enter = 0, exit_z = 0, thick = 0;
        int found = 0;
        if (seams_brine_fronts(front, &enter, &exit_z, &thick, &found) != SEAMS_OK) {
          fault("brine fronts");
        }
        put_i("brine_found", found);
        put_d("brine_enter", enter);
        put_d("brine_exit", exit_z);
        put_d("brine_thick", thick);
        double c_ice = 0, c_liquid = 0, k = 0;
        int defined = 0;
        (void)seams_brine_rejection(front, &c_ice, &c_liquid, &k, &defined);
        size_t bins = 0;
        (void)seams_brine_profile(front, NULL, NULL, NULL, 0, &bins);
        int pockets = 0, reservoir = 0;
        (void)seams_brine_pockets(front, &pockets, &reservoir);
        seams_brine_destroy(front);
        seams_brine_destroy(front);
      }
    }
    free(ice);
  }
  free(labels);
  free(q4);
  free(q6);
  free(q6bar);
}

static void compare_line(const char *who, const char *key, const char *value, int n, double tol) {
  const Rec *rec = find_rec(key);
  if (rec == NULL) {
    return;
  }
  if (rec->is_int) {
    if (atoi(value) != rec->ival) {
      fprintf(stderr, "fail: %s %s C %d other %s\n", who, key, rec->ival, value);
      g_fails += 1;
    }
  } else {
    double got = strtod(value, NULL);
    double limit = tol;
    if (strstr(key, "_fp") != NULL) {
      limit = 1e-6 * (double)n * ((double)n + 1.0) / 2.0;
    }
    if (!isfinite(got) || fabs(got - rec->dval) > limit) {
      fprintf(stderr, "fail: %s %s C %.17g other %.17g\n", who, key, rec->dval, got);
      g_fails += 1;
    }
  }
}

static void compare_cmd(const char *cmd, const char *who, int n, double tol, const char *only) {
  FILE *pipe = popen(cmd, "r");
  if (pipe == NULL) {
    fault(who);
    return;
  }
  char line[256];
  while (fgets(line, sizeof line, pipe) != NULL) {
    char key[32], value[128];
    if (sscanf(line, "%31s %127s", key, value) == 2) {
      if (only != NULL && strncmp(key, only, strlen(only)) != 0 &&
          strncmp(key, "knn_", 4) != 0 && strncmp(key, "rdf_", 4) != 0 &&
          strncmp(key, "cn", 2) != 0 && strncmp(key, "hbond_", 6) != 0 &&
          strncmp(key, "dens_", 5) != 0 && strncmp(key, "prism_", 6) != 0 &&
          strncmp(key, "mono_", 5) != 0 && strncmp(key, "site_", 5) != 0 &&
          strncmp(key, "contact_", 8) != 0 && strncmp(key, "domain_", 7) != 0 &&
          strncmp(key, "fp_", 3) != 0 && strcmp(key, "ovito_q4") != 0 &&
          strcmp(key, "ovito_q6") != 0 && strcmp(key, "freud_q4") != 0 &&
          strcmp(key, "freud_q6") != 0) {
        continue;
      }
      if (strcmp(key, "ovito_q4") == 0) {
        compare_line(who, "q4_mean", value, n, 1e-6);
      } else if (strcmp(key, "ovito_q6") == 0) {
        compare_line(who, "q6_mean", value, n, 1e-6);
      } else if (strcmp(key, "freud_q4") == 0) {
        compare_line(who, "q4_mean", value, n, 1e-4);
      } else if (strcmp(key, "freud_q6") == 0) {
        compare_line(who, "q6_mean", value, n, 1e-4);
      } else {
        compare_line(who, key, value, n, tol);
      }
    }
  }
  int status = pclose(pipe);
  if (status != 0 && strcmp(who, "c++") == 0) {
    fprintf(stderr, "fail: %s exited %d\n", who, status);
    g_fails += 1;
  }
}

static void threads(seams_system *sys, int n) {
  double *base = calloc((size_t)n, sizeof(double));
  if (seams_set_num_threads(1) != SEAMS_OK ||
      seams_q4q6(sys, 3.5, 0, NULL, base) != SEAMS_OK) {
    fault("serial q6");
    free(base);
    return;
  }
  int widths[4] = {1, 2, 4, 8};
  for (int i = 0; i < 4; ++i) {
    double *q6 = calloc((size_t)n, sizeof(double));
    if (seams_set_num_threads(widths[i]) != SEAMS_OK || seams_get_num_threads() != widths[i] ||
        seams_q4q6(sys, 3.5, 0, NULL, q6) != SEAMS_OK) {
      fault("threaded q6");
    } else if (fabs(fingerprint(q6, n) - fingerprint(base, n)) >
               1e-6 * (double)n * ((double)n + 1) / 2.0) {
      fprintf(stderr, "fail: q6 drifted at %d threads\n", widths[i]);
      g_fails += 1;
    }
    free(q6);
  }
  seams_set_num_threads(1);
  free(base);
}

static void one(const char *src, const char *cpp, const char *py, const char *work,
                const char *name, const Frame *frame, const char *read_path) {
  char path[1024];
  snprintf(path, sizeof path, "%s/%s.frame", work, name);
  if (write_frame(path, frame) != 0) {
    fault("write frame");
    return;
  }
  fprintf(stderr, "capi fixture %s n=%d\n", name, frame->n);
  seams_system *sys = seams_system_create();
  if (!store(sys, frame)) {
    fault("store");
    seams_system_destroy(sys);
    return;
  }
  if (read_path != NULL) {
    size_t nframes = 0;
    if (seams_path_nframes(read_path, &nframes) != SEAMS_OK || nframes < 1) {
      fault("nframes");
    }
    if (seams_system_read(sys, read_path, 1, 0) != SEAMS_OK) {
      fault("read");
    }
    if (!store(sys, frame)) {
      fault("restore after read");
    }
  }
  call_every(sys, frame, name);
  threads(sys, frame->n);
  seams_system_destroy(sys);
  char cmd[4096];
  const char *only = strcmp(name, "brine") == 0 ? "brine_" : NULL;
  snprintf(cmd, sizeof cmd, "\"%s\" \"%s\" %s", cpp, path, name);
  compare_cmd(cmd, "c++", frame->n, 1e-6, only);
  /* Neighbour lists, Steinhardt parameters, CHILL+, template labels,
   * the RDF, and coordination are compared with freud and OVITO by
   * tests/capi_ovito_freud.py. That program calls this library through
   * ctypes. A missing import used to print nothing and exit 0. */
  (void)py;
  (void)src;
}

int main(int argc, char **argv) {
  if (argc != 5) {
    fprintf(stderr, "usage: capi_suite SRC CPP OVITO_PY WORK\n");
    return 2;
  }
  mkdir(argv[4], 0755);
  /* The reference process stays serial. The thread loop below is the
   * OpenMP check, at 1, 2, 4 and 8 threads. */
  setenv("OMP_NUM_THREADS", "1", 1);
  setenv("RAYON_NUM_THREADS", "1", 1);
  (void)seams_set_num_threads(1);
  names_and_versions();
  Frame ic;
  memset(&ic, 0, sizeof ic);
  ice_ic(&ic);
  error_paths(&ic);
  char ice_path[1024], si_path[1024], liquid_path[1024], brine_path[1024], nuc_path[1024],
      xyz_path[1024], con_path[1024];
  snprintf(ice_path, sizeof ice_path, "%s/input/traj/ice1h.lammpstrj", argv[1]);
  snprintf(si_path, sizeof si_path, "%s/input/traj/genice_sI.lammpstrj", argv[1]);
  snprintf(liquid_path, sizeof liquid_path, "%s/benchmarks/paper-cells/tip4pice-liquid-4096.txt",
           argv[1]);
  snprintf(brine_path, sizeof brine_path, "%s/input/traj/nacl_front.lammpstrj", argv[1]);
  snprintf(nuc_path, sizeof nuc_path, "%s/input/traj/nucleation_f1.lammpstrj", argv[1]);
  snprintf(xyz_path, sizeof xyz_path, "%s/input/xyz/configuration.xyz", argv[1]);
  snprintf(con_path, sizeof con_path, "%s/input/con/tiny_multi_cuh2.con", argv[1]);
  char mmcif[1024];
  snprintf(mmcif, sizeof mmcif, "%s/empty.mmcif", argv[4]);
  FILE *empty = fopen(mmcif, "w");
  if (empty != NULL) {
    fputs(" \n", empty);
    fclose(empty);
  }
  seams_system *mmcif_sys = seams_system_create();
  size_t nframes = 0;
  if (seams_path_nframes(mmcif, &nframes) == SEAMS_OK || seams_last_error()[0] == '\0' ||
      seams_system_read(mmcif_sys, mmcif, 1, 0) == SEAMS_OK || seams_last_error()[0] == '\0') {
    fault("empty mmcif");
  }
  seams_system_destroy(mmcif_sys);
  seams_system *reader = seams_system_create();
  if (seams_system_read(reader, ice_path, 1, 0) != SEAMS_OK ||
      seams_system_read(reader, xyz_path, 1, -1) != SEAMS_OK) {
    fault("fixture read");
  }
  if (seams_system_read(reader, con_path, 1, -1) != SEAMS_OK && seams_last_error()[0] == '\0') {
    fault("con read message");
  }
  (void)seams_path_nframes(ice_path, &nframes);
  seams_system_destroy(reader);
  Frame ih, si, liquid, brine, nuc;
  memset(&ih, 0, sizeof ih);
  memset(&si, 0, sizeof si);
  memset(&liquid, 0, sizeof liquid);
  memset(&brine, 0, sizeof brine);
  memset(&nuc, 0, sizeof nuc);
  if (load_dump(ice_path, &ih) != 0 || load_dump(si_path, &si) != 0 ||
      load_oxygen(liquid_path, &liquid) != 0 || load_dump(brine_path, &brine) != 0 ||
      load_dump(nuc_path, &nuc) != 0) {
    fault("load fixtures");
  } else {
    one(argv[1], argv[2], argv[3], argv[4], "ice_ic", &ic, ice_path);
    one(argv[1], argv[2], argv[3], argv[4], "ice_ih", &ih, NULL);
    one(argv[1], argv[2], argv[3], argv[4], "si", &si, NULL);
    one(argv[1], argv[2], argv[3], argv[4], "liquid", &liquid, NULL);
    one(argv[1], argv[2], argv[3], argv[4], "brine", &brine, NULL);
    one(argv[1], argv[2], argv[3], argv[4], "nucleation", &nuc, NULL);
  }
  frame_clear(&ic);
  frame_clear(&ih);
  frame_clear(&si);
  frame_clear(&liquid);
  frame_clear(&brine);
  frame_clear(&nuc);
  if (g_fails != 0) {
    fprintf(stderr, "%d capi checks failed\n", g_fails);
    return 1;
  }
  printf("capi suite ok\n");
  return 0;
}
