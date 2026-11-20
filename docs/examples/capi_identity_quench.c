// Copyright (c) d-SEAMS developers
/* Export hydrogen-bond edges for a real fixture frame, then run the
 * quenched-analysis hook with an identity callback. An identity callback
 * returns its input unchanged, so the quenched frame is bit-identical to
 * the raw frame and every raw-vs-quenched delta must be exactly zero.
 */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static seams_status identity_quench(const double *xyz, size_t n, const double *box9,
                                    size_t box_len, double *xyz_out, size_t xyz_out_len,
                                    double *box9_out, size_t box9_out_len, void *user_data) {
  (void)user_data;
  if (xyz == NULL || box9 == NULL || xyz_out == NULL || xyz_out_len < 3 * n || box_len != 9) {
    return SEAMS_ERR_NULL;
  }
  memcpy(xyz_out, xyz, sizeof(double) * 3 * n);
  if (box9_out != NULL && box9_out_len == 9) {
    memcpy(box9_out, box9, sizeof(double) * 9);
  }
  return SEAMS_OK;
}

int main(int argc, char **argv) {
  const char *root = argc > 1 ? argv[1] : ".";
  char sample[4096];
  snprintf(sample, sizeof sample, "%s/input/traj/exampleTraj.lammpstrj", root);

  seams_system *sys = seams_system_create();
  if (sys == NULL || seams_system_read(sys, sample, 1, -1) != SEAMS_OK) {
    fprintf(stderr, "read: %s\n", seams_last_error());
    return 1;
  }

  const double cutoff = 3.5;
  const int oxygen_type = 2;
  const int hydrogen_type = 1;
  const double dist_cutoff = 2.42;
  const double angle_cutoff = 30.0;

  size_t n_edges = 0;
  if (seams_hbond_edges(sys, cutoff, oxygen_type, hydrogen_type, dist_cutoff, angle_cutoff,
                        SEAMS_HBOND_EDGE_WITH_H | SEAMS_HBOND_EDGE_WITH_GEOMETRY, NULL, 0,
                        &n_edges) != SEAMS_OK) {
    fprintf(stderr, "edge query: %s\n", seams_last_error());
    return 1;
  }
  seams_hbond_edge *edges = calloc(n_edges > 0 ? n_edges : 1, sizeof(seams_hbond_edge));
  if (edges == NULL) {
    return 1;
  }
  for (size_t i = 0; i < (n_edges > 0 ? n_edges : 1); ++i) {
    edges[i].size = sizeof(seams_hbond_edge);
  }
  if (seams_hbond_edges(sys, cutoff, oxygen_type, hydrogen_type, dist_cutoff, angle_cutoff,
                        SEAMS_HBOND_EDGE_WITH_H | SEAMS_HBOND_EDGE_WITH_GEOMETRY, edges, n_edges,
                        &n_edges) != SEAMS_OK) {
    fprintf(stderr, "edges: %s\n", seams_last_error());
    return 1;
  }
  printf("capi_identity_quench edges %zu\n", n_edges);
  if (n_edges == 0) {
    fprintf(stderr, "expected a real fixture frame with hydrogen bonds\n");
    return 1;
  }

  const int n = seams_system_count(sys);
  int *labels_raw = calloc((size_t)n, sizeof(int));
  int *labels_q = calloc((size_t)n, sizeof(int));
  int counts_raw[9] = {0};
  int counts_q[9] = {0};
  int rings_raw[8] = {0};
  int rings_q[8] = {0};
  int hbond_raw = -1;
  int hbond_q = -1;
  int edge_delta = -1;
  int ring_delta = -1;
  int class_delta = -1;
  if (labels_raw == NULL || labels_q == NULL) {
    return 1;
  }
  if (seams_quenched_analysis(sys, identity_quench, NULL, cutoff, oxygen_type, SEAMS_GRAPH_CUTOFF,
                              6, oxygen_type, hydrogen_type, dist_cutoff, angle_cutoff, labels_raw,
                              (size_t)n, labels_q, (size_t)n, counts_raw, counts_q, rings_raw,
                              rings_q, 8, &hbond_raw, &hbond_q, &edge_delta, &ring_delta,
                              &class_delta) != SEAMS_OK) {
    fprintf(stderr, "quenched analysis: %s\n", seams_last_error());
    return 1;
  }
  printf("capi_identity_quench raw_hbonds %d quenched_hbonds %d\n", hbond_raw, hbond_q);
  printf("capi_identity_quench edge_delta %d ring_delta %d class_delta %d\n", edge_delta,
         ring_delta, class_delta);

  free(labels_raw);
  free(labels_q);
  free(edges);
  seams_system_destroy(sys);

  if (hbond_raw != (int)n_edges || hbond_raw != hbond_q || edge_delta != 0 || ring_delta != 0 ||
      class_delta != 0) {
    fprintf(stderr, "identity callback must give zero deltas\n");
    return 1;
  }
  return 0;
}
