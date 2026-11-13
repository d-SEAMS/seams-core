/* Copyright (c) d-SEAMS developers */
#include "seams.h"

#include <math.h>
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

static int diamond(int nrep, double *xyz, double box[9]) {
  const double a = 2.75 * 4.0 / sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  int n = 0;
  for (int i = 0; i < nrep; i++) {
    for (int j = 0; j < nrep; j++) {
      for (int k = 0; k < nrep; k++) {
        for (int f = 0; f < 4; f++) {
          for (int s = 0; s < 2; s++) {
            xyz[3 * n] = (i + fcc[f][0] + basis[s][0]) * a;
            xyz[3 * n + 1] = (j + fcc[f][1] + basis[s][1]) * a;
            xyz[3 * n + 2] = (k + fcc[f][2] + basis[s][2]) * a;
            n += 1;
          }
        }
      }
    }
  }
  box[0] = box[1] = box[2] = 0.0;
  box[3] = box[4] = box[5] = nrep * a;
  box[6] = box[7] = box[8] = 0.0;
  return n;
}

int main(void) {
  expect(seams_abi() == SEAMS_ABI, "abi");
  expect(seams_version_major() == SEAMS_VERSION_MAJOR, "major");
  expect(seams_version_minor() == SEAMS_VERSION_MINOR, "minor");
  expect(seams_version_patch() == SEAMS_VERSION_PATCH, "patch");
  expect(strcmp(seams_version_string(), SEAMS_VERSION_STRING) == 0, "version string");
  expect(strcmp(seams_status_string(SEAMS_OK), "ok") == 0, "status ok");
  expect(strcmp(seams_status_string(SEAMS_ERR_NULL), "null pointer") == 0, "status null");
  expect(strcmp(seams_status_string(SEAMS_ERR_RANGE), "value out of range") == 0, "status range");
  expect(strcmp(seams_status_string(SEAMS_ERR_BOX), "box is not a positive cell") == 0, "status box");
  expect(strcmp(seams_status_string(SEAMS_ERR_ARG), "bad argument") == 0, "status arg");
  expect(strcmp(seams_status_string(SEAMS_ERR_STATE), "system has no frame") == 0, "status state");
  expect(strcmp(seams_status_string(SEAMS_ERR_NOMEM), "out of memory") == 0, "status nomem");
  expect(strcmp(seams_status_string(SEAMS_ERR_INTERNAL), "internal error") == 0, "status internal");
  expect(strcmp(seams_label_name(SEAMS_LABEL_CUBIC), "cubic") == 0, "label cubic");
  expect(strcmp(seams_label_name(SEAMS_LABEL_HEXAGONAL), "hexagonal") == 0, "label hexagonal");
  expect(strcmp(seams_label_name(SEAMS_LABEL_WATER), "water") == 0, "label water");
  expect(strcmp(seams_label_name(SEAMS_LABEL_INTERFACIAL), "interfacial") == 0, "label interfacial");
  expect(strcmp(seams_label_name(SEAMS_LABEL_CLATHRATE), "clathrate") == 0, "label clathrate");
  expect(strcmp(seams_label_name(SEAMS_LABEL_INTERCLATHRATE), "interClathrate") == 0,
         "label interClathrate");
  expect(strcmp(seams_label_name(SEAMS_LABEL_UNCLASSIFIED), "unclassified") == 0,
         "label unclassified");
  expect(strcmp(seams_label_name(SEAMS_LABEL_RECUBIC), "reCubic") == 0, "label reCubic");
  expect(strcmp(seams_label_name(SEAMS_LABEL_REHEX), "reHex") == 0, "label reHex");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_UNCLASSIFIED), "unclassified") == 0,
         "ring unclassified");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_DDC), "ddc") == 0, "ring ddc");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_HC_BASAL), "hc-basal") == 0, "ring hc basal");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_HC_PRISMATIC), "hc-prismatic") == 0,
         "ring hc prismatic");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_BOTH_BASAL), "both-basal") == 0,
         "ring both basal");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_BOTH_PRISMATIC), "both-prismatic") == 0,
         "ring both prismatic");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_PRISM), "prism") == 0, "ring prism");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_DEFORMED_PRISM), "deformed-prism") == 0,
         "ring deformed");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_MIXED_PRISM), "mixed-prism") == 0,
         "ring mixed");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_PENTAGON), "pentagon") == 0, "ring pentagon");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_HEXAGON), "hexagon") == 0, "ring hexagon");
  expect(strcmp(seams_ring_symmetry_name(SEAMS_RING_OTHER), "other") == 0, "ring other");
  expect(seams_last_error() != NULL, "last error pointer");
  const int original_threads = seams_get_num_threads();
  expect(seams_set_num_threads(2) == SEAMS_OK, "set threads");
  expect(seams_get_num_threads() == 2, "get threads");
  expect(seams_set_num_threads(-1) == SEAMS_ERR_RANGE, "negative threads");
  expect(seams_set_num_threads(1) == SEAMS_OK, "one thread");
  expect(seams_set_num_threads(0) == SEAMS_OK, "restore threads");
  expect(seams_get_num_threads() == original_threads, "restored omp max");

  const int nrep = 2;
  const int cap = nrep * nrep * nrep * 8;
  double *xyz = calloc((size_t)cap * 3, sizeof(double));
  double box[9];
  const int n = diamond(nrep, xyz, box);
  expect(n == cap, "diamond count");

  seams_system *sys = seams_system_create();
  expect(sys != NULL, "create");
  expect(seams_system_set_frame(NULL, xyz, n, box, NULL, NULL) == SEAMS_ERR_NULL,
         "null system");
  box[3] = 0.0;
  expect(seams_system_set_frame(sys, xyz, n, box, NULL, NULL) == SEAMS_ERR_BOX,
         "flat box");
  box[3] = box[4];
  expect(seams_system_set_ortho(sys, xyz, n, 0, 0, 0, box[3], box[4], box[5], NULL,
                                NULL) == SEAMS_OK,
         "ortho frame");
  expect(seams_system_count(sys) == n, "count");

  seams_neighbours *nb = NULL;
  expect(seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK,
         "neighbours");
  const int *offsets = NULL;
  const int *cols = NULL;
  int n_offsets = 0;
  int nnz = 0;
  expect(seams_neighbours_csr(nb, &offsets, &n_offsets, &cols, &nnz) == SEAMS_OK, "csr");
  expect(n_offsets == n + 1, "csr rows");
  expect(seams_neighbours_count(nb) == n, "neighbour rows");
  int degree_ok = 1;
  for (int i = 0; i < n; i++) {
    if (offsets[i + 1] - offsets[i] != 4) {
      degree_ok = 0;
    }
  }
  expect(degree_ok, "four neighbours");
  seams_neighbours_destroy(nb);

  int *labels = calloc((size_t)n, sizeof(int));
  int counts[9];
  expect(seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                          counts) == SEAMS_OK,
         "published chill+");
  int cubic = 1;
  for (int i = 0; i < n; i++) {
    if (labels[i] != SEAMS_LABEL_CUBIC) {
      cubic = 0;
    }
  }
  expect(cubic && counts[SEAMS_LABEL_CUBIC] == n, "published labels are cubic");
  expect(seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_MUTUAL, SEAMS_CHILL_OVITO, labels,
                          counts) == SEAMS_OK,
         "ovito chill+");
  cubic = 1;
  for (int i = 0; i < n; i++) {
    if (labels[i] != SEAMS_LABEL_CUBIC) {
      cubic = 0;
    }
  }
  expect(cubic, "ovito labels are cubic");

  double *q4 = calloc((size_t)n, sizeof(double));
  double *q6 = calloc((size_t)n, sizeof(double));
  expect(seams_q4q6(sys, 3.5, 0, q4, q6) == SEAMS_OK, "q4 q6");
  double mean6 = 0.0;
  int finite = 1;
  for (int i = 0; i < n; i++) {
    if (!isfinite(q4[i]) || !isfinite(q6[i])) {
      finite = 0;
    }
    mean6 += q6[i];
  }
  mean6 /= n;
  expect(finite && mean6 > 0.3, "q6 of cubic ice");

  box[6] = 0.4;
  expect(seams_system_set_frame(sys, xyz, n, box, NULL, NULL) == SEAMS_OK, "triclinic");
  expect(seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                          counts) == SEAMS_OK,
         "triclinic chill+");
  box[6] = 0.0;
  expect(seams_system_set_frame(sys, xyz, n, box, NULL, NULL) == SEAMS_OK, "restore ortho");

  seams_rings *rings = NULL;
  expect(seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK, "rings");
  int hist[7];
  expect(seams_rings_histogram(rings, hist, 7) == SEAMS_OK, "histogram");
  expect(hist[6] > 0, "six-rings");
  const int nr = seams_rings_count(rings);
  int *sizes = calloc((size_t)nr, sizeof(int));
  int *sym = calloc((size_t)nr, sizeof(int));
  expect(seams_rings_labels(rings, sizes, sym) == SEAMS_OK, "symmetry");
  int labelled = 0;
  for (int i = 0; i < nr; i++) {
    if (sizes[i] == 6 && sym[i] != SEAMS_RING_OTHER) {
      labelled += 1;
    }
  }
  expect(labelled == hist[6], "every six-ring has a symmetry label");
  double fraction = 0.0;
  expect(seams_rings_fraction(rings, 6, &fraction) == SEAMS_OK, "fraction");
  expect(fraction > 0.5, "most atoms sit on a six-ring");

  seams_rings_destroy(rings);

  {
    double pair[9] = {0, 0, 0, 2.0, 0, 0, 12.0, 0, 0};
    int types[3] = {2, 2, 1};
    double wide[9] = {0, 0, 0, 30, 30, 30, 0, 0, 0};
    expect(seams_system_set_frame(sys, pair, 3, wide, types, NULL) == SEAMS_OK,
           "two-type frame");
    double qa[3], qb[3];
    expect(seams_q4q6(sys, 3.5, 0, qa, qb) == SEAMS_OK, "type<=0 uses a copy");
    double q2[3] = {0, 0, 0};
    double q6b[3] = {0, 0, 0};
    expect(seams_q4q6(sys, 3.5, 2, q2, q6b) == SEAMS_OK, "type 2 after unify");
    expect(q2[0] > 0.2 && q2[1] > 0.2, "type 2 pair was not rewritten");
    expect(!(q2[2] > 0.2), "the other type stays out of the type 2 shell");
  }

  seams_system_destroy(sys);
  free(xyz);
  free(labels);
  free(q4);
  free(q6);
  free(sizes);
  free(sym);
  if (g_fails != 0) {
    fprintf(stderr, "%d checks failed\n", g_fails);
    return 1;
  }
  printf("seams_c_api ok %s\n", seams_version_string());
  return 0;
}
