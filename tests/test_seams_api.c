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
  expect(seams_system_count(sys) == n, "clear the error before destroy");
  seams_neighbours_destroy(nb);
  expect(strcmp(seams_last_error(), "double destroy") != 0, "first destroy frees the list");
  seams_neighbours_destroy(nb);
  expect(strcmp(seams_last_error(), "double destroy") == 0, "second destroy is refused");

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

  expect(seams_version() == SEAMS_VERSION_NUMBER, "packed version");
  expect(seams_system_set_frame(sys, xyz, n, box, NULL, NULL) == SEAMS_OK, "restore diamond");
  expect(seams_system_set_atoms(sys, xyz, (size_t)n, box, 9, NULL, NULL) == SEAMS_OK,
         "set atoms");
  expect(seams_system_set_atoms(sys, xyz, (size_t)n, box, 8, NULL, NULL) == SEAMS_ERR_RANGE,
         "box length");
  expect(seams_handle_error(sys)[0] != '\0', "handle keeps the box error");
  expect(seams_system_set_atoms(sys, xyz, (size_t)n, box, 9, NULL, NULL) == SEAMS_OK,
         "set atoms again");

  seams_neighbours *copy_nb = NULL;
  expect(seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &copy_nb) == SEAMS_OK,
         "neighbours again");
  size_t n_off = 0;
  size_t n_col = 0;
  expect(seams_neighbours_copy(copy_nb, NULL, 0, &n_off, NULL, 0, &n_col) == SEAMS_OK,
         "csr query");
  expect(n_off == (size_t)n + 1, "query rows");
  int *off_copy = calloc(n_off, sizeof(int));
  int *col_copy = calloc(n_col > 0 ? n_col : 1, sizeof(int));
  expect(seams_neighbours_copy(copy_nb, off_copy, n_off, &n_off, col_copy, n_col, &n_col) ==
             SEAMS_OK,
         "csr copy");
  expect(seams_neighbours_copy(copy_nb, off_copy, 1, &n_off, col_copy, n_col, &n_col) ==
             SEAMS_ERR_RANGE,
         "short csr buffer");
  seams_neighbours_destroy(copy_nb);

  double *q3 = calloc((size_t)n, sizeof(double));
  double *q3bar = calloc((size_t)n, sizeof(double));
  double *ql4 = calloc((size_t)n, sizeof(double));
  double *ql6 = calloc((size_t)n, sizeof(double));
  expect(seams_ql(sys, 3.5, 0, 3, SEAMS_GRAPH_CUTOFF, q3, (size_t)n, q3bar, (size_t)n) ==
             SEAMS_OK,
         "q3");
  expect(seams_ql(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, ql4, (size_t)n, NULL, 0) == SEAMS_OK,
         "q4");
  expect(seams_ql(sys, 3.5, 0, 6, SEAMS_GRAPH_MUTUAL, ql6, (size_t)n, NULL, 0) == SEAMS_OK,
         "q6 mutual");
  int q_ok = 1;
  for (int i = 0; i < n; i++) {
    if (!isfinite(q3[i]) || !isfinite(ql4[i]) || !isfinite(ql6[i])) {
      q_ok = 0;
    }
  }
  expect(q_ok, "q3 q4 q6 are finite");

  double *q6bar = calloc((size_t)n, sizeof(double));
  double q6_mean = 0.0;
  double q6_median = 0.0;
  int q6_ice = 0;
  expect(seams_q6bar(sys, 3.5, 0, 0.36, q6bar, (size_t)n, &q6_mean, &q6_median, &q6_ice) ==
             SEAMS_OK,
         "q6bar");
  expect(isfinite(q6_mean) && isfinite(q6_median) && q6_ice >= 0, "q6bar summary");

  expect(seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, (size_t)n,
                     counts) == SEAMS_OK,
         "published chill");
  expect(counts[SEAMS_LABEL_CUBIC] > 0, "published chill finds cubic ice");
  expect(seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, labels, (size_t)n,
                     counts) == SEAMS_OK,
         "ovito chill");

  seams_rings *again = NULL;
  expect(seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &again) == SEAMS_OK,
         "rings again");
  size_t n_vert = 0;
  expect(seams_rings_vertices(again, 0, NULL, 0, &n_vert) == SEAMS_OK, "ring query");
  int *verts = calloc(n_vert > 0 ? n_vert : 1, sizeof(int));
  expect(seams_rings_vertices(again, 0, verts, n_vert, &n_vert) == SEAMS_OK, "ring vertices");
  seams_rings_destroy(again);

  seams_cage_set *cages = NULL;
  expect(seams_cages(sys, 3.5, 0, SEAMS_GRAPH_SEEDED, 4, 0, &cages) == SEAMS_OK,
         "seeded cages");
  int n_hex = 0;
  int n_cubic = 0;
  int n_water = 0;
  int n_both = 0;
  expect(seams_cages_counts(cages, &n_hex, &n_cubic, &n_water, &n_both) == SEAMS_OK,
         "cage counts");
  expect(n_hex + n_cubic + n_water + n_both == n, "cage counts cover the frame");
  int *per_atom = calloc((size_t)n, sizeof(int));
  expect(seams_cages_labels(cages, per_atom, (size_t)n) == SEAMS_OK, "cage labels");
  seams_cage_label tag = (seams_cage_label)per_atom[0];
  expect(tag == SEAMS_CAGE_WATER || tag == SEAMS_CAGE_HEX || tag == SEAMS_CAGE_CUBIC ||
             tag == SEAMS_CAGE_BOTH,
         "cage code");
  seams_cages_destroy(cages);
  size_t n_sig = 0;
  size_t n_sig_atoms = 0;
  expect(seams_cages_signature(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", &n_sig, &n_sig_atoms) ==
             SEAMS_OK,
         "signature alias");
  expect(seams_signature_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", &n_sig, &n_sig_atoms) ==
             SEAMS_OK,
         "signature cages");

  int *shell = calloc((size_t)n, sizeof(int));
  expect(seams_system_symmetry(sys, 0.0, shell, (size_t)n) == SEAMS_OK, "symmetry");
  int shell_ok = 1;
  for (int i = 0; i < n; i++) {
    if (shell[i] < -1 || shell[i] > 2) {
      shell_ok = 0;
    }
  }
  expect(shell_ok, "shell labels");

  int *ice = calloc((size_t)n, sizeof(int));
  for (int i = 0; i < n; i++) {
    ice[i] = 1;
  }
  seams_brine *front = NULL;
  expect(seams_brine_analyse(sys, ice, (size_t)n, 0, 0, 4, 2, &front) == SEAMS_OK, "brine");
  double z_enter = 0.0;
  double z_exit = 0.0;
  double thickness = 0.0;
  int found = 0;
  expect(seams_brine_fronts(front, &z_enter, &z_exit, &thickness, &found) == SEAMS_OK,
         "fronts");
  double c_ice = 0.0;
  double c_liquid = 0.0;
  double rejection = 0.0;
  int defined = 0;
  expect(seams_brine_rejection(front, &c_ice, &c_liquid, &rejection, &defined) == SEAMS_OK,
         "rejection");
  size_t n_bins = 0;
  expect(seams_brine_profile(front, NULL, NULL, NULL, 0, &n_bins) == SEAMS_OK, "profile query");
  double *z_mid = calloc(n_bins > 0 ? n_bins : 1, sizeof(double));
  double *molality = calloc(n_bins > 0 ? n_bins : 1, sizeof(double));
  double *ice_fraction = calloc(n_bins > 0 ? n_bins : 1, sizeof(double));
  expect(seams_brine_profile(front, z_mid, molality, ice_fraction, n_bins, &n_bins) == SEAMS_OK,
         "profile");
  int n_pockets = 0;
  int reservoir = 0;
  expect(seams_brine_pockets(front, &n_pockets, &reservoir) == SEAMS_OK, "pockets");
  seams_brine_destroy(front);

  size_t nframes = 0;
  expect(seams_path_nframes(NULL, &nframes) == SEAMS_ERR_NULL, "null path");
  const char *sample = getenv("SEAMS_SAMPLE");
  if (sample != NULL) {
    expect(seams_path_nframes(sample, &nframes) == SEAMS_OK, "nframes");
    expect(nframes >= 1, "at least one frame");
    seams_system *file = seams_system_create();
    expect(seams_system_read(file, sample, 1, -1) == SEAMS_OK, "read frame");
    expect(seams_system_count(file) > 0, "read count");
    seams_system_destroy(file);
  }
  const char *xyz_path = getenv("SEAMS_XYZ");
  if (xyz_path != NULL) {
    seams_system *xyz_frame = seams_system_create();
    expect(seams_system_read(xyz_frame, xyz_path, 1, -1) == SEAMS_OK, "read xyz");
    expect(seams_system_count(xyz_frame) > 0, "xyz count");
    seams_system_destroy(xyz_frame);
  }
  const char *con_path = getenv("SEAMS_CON");
  if (con_path != NULL) {
    size_t con_frames = 0;
    expect(seams_path_nframes(con_path, &con_frames) == SEAMS_OK, "con frames");
    expect(con_frames >= 1, "con has a frame");
    seams_system *con_frame = seams_system_create();
    expect(seams_system_read(con_frame, con_path, 1, -1) == SEAMS_OK, "read con");
    expect(seams_system_count(con_frame) > 0, "con count");
    seams_system_destroy(con_frame);
  }
  {
    FILE *dump = fopen("capi_triclinic.lammpstrj", "w");
    expect(dump != NULL, "write triclinic dump");
    if (dump != NULL) {
      fputs("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n"
            "ITEM: BOX BOUNDS xy xz yz pp pp pp\n0.0 10.0 0.5\n0.0 10.0 0.2\n"
            "0.0 10.0 0.1\nITEM: ATOMS id type x y z\n1 1 1.0 1.0 1.0\n"
            "2 1 2.0 1.0 1.0\n",
            dump);
      fclose(dump);
      seams_system *tri = seams_system_create();
      expect(seams_system_read(tri, "capi_triclinic.lammpstrj", 1, -1) == SEAMS_OK,
             "read triclinic dump");
      expect(seams_system_count(tri) == 2, "triclinic count");
      seams_system_destroy(tri);
      remove("capi_triclinic.lammpstrj");
    }
  }

  expect(seams_steinhardt(sys, 3.5, 0, 6, q6bar, (size_t)n, NULL, 0) == SEAMS_OK, "steinhardt");
  seams_neighbours *knn = NULL;
  expect(seams_neighbours_knn(sys, 4, 3.5, 0, 1, &knn) == SEAMS_OK, "knn mutual");
  expect(seams_neighbours_count(knn) == n, "knn rows");
  seams_neighbours_destroy(knn);
  expect(seams_neighbours_knn(sys, 4, 3.5, 0, 0, &knn) == SEAMS_OK, "knn one way");
  seams_neighbours_destroy(knn);
  size_t n_rdf = 0;
  double volume = 0.0;
  expect(seams_rdf(sys, 6.0, 20, 0, 0, NULL, NULL, NULL, 0, &n_rdf, &volume) == SEAMS_OK,
         "rdf query");
  expect(n_rdf > 0 && volume > 0.0, "rdf bins");
  double *rdf_r = calloc(n_rdf, sizeof(double));
  double *rdf_g = calloc(n_rdf, sizeof(double));
  expect(seams_rdf(sys, 6.0, 20, 0, 0, rdf_r, rdf_g, NULL, n_rdf, &n_rdf, NULL) == SEAMS_OK,
         "rdf");
  size_t n_cn = 0;
  expect(seams_running_cn(sys, 6.0, 20, 0, 0, NULL, NULL, 0, &n_cn) == SEAMS_OK, "running cn");
  double cn = 0.0;
  expect(seams_coordination(sys, 3.5, 20, 0, 0, &cn) == SEAMS_OK, "coordination");
  expect(isfinite(cn) && cn > 0.0, "coordination value");
  int n_hbonds = -1;
  expect(seams_hbonds(sys, 3.5, 1, 2, 2.42, 30.0, &n_hbonds) == SEAMS_OK, "hbonds");
  expect(n_hbonds == 0, "no hydrogens");
  size_t n_rho = 0;
  expect(seams_density_axis(sys, 0, 8, 2, NULL, NULL, 0, &n_rho) == SEAMS_OK, "density");
  expect(n_rho == 8, "density bins");
  seams_topo *fp = NULL;
  expect(seams_fingerprint(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 0, &fp) == SEAMS_OK,
         "fingerprint");
  size_t key_len = 0;
  expect(seams_fingerprint_key(fp, NULL, 0, &key_len) == SEAMS_OK, "frame key");
  expect(key_len > 0, "frame key length");
  char *key = calloc(key_len + 1, 1);
  expect(seams_fingerprint_key(fp, key, key_len + 1, &key_len) == SEAMS_OK, "frame key text");
  expect(seams_fingerprint_atom(fp, 0, key, key_len + 1, &key_len) == SEAMS_OK, "atom key");
  int fp_rings[8];
  expect(seams_fingerprint_rings(fp, fp_rings, 8) == SEAMS_OK, "fingerprint rings");
  seams_fingerprint_destroy(fp);
  int kinds[64];
  expect(seams_site_table(sys, "1=waterO", kinds, (size_t)n) == SEAMS_OK, "site table");
  seams_site_kind kind = (seams_site_kind)kinds[0];
  expect(kind == SEAMS_SITE_WATER_O, "water oxygen kind");
  size_t n_pairs = 0;
  int n_cation = 0;
  int n_anion = 0;
  expect(seams_contact_pairs(sys, "1=cationHead", &n_pairs, &n_cation, &n_anion) == SEAMS_OK,
         "contact pairs");
  int n_subset = 0;
  int largest = 0;
  double percolation = 0.0;
  expect(seams_domains(sys, 3.5, "1=waterO", 0, &n_subset, &largest, &percolation) == SEAMS_OK,
         "domains");
  int perfect[7];
  int imperfect[7];
  expect(seams_prisms(sys, 3.5, 0, perfect, imperfect, 7) == SEAMS_OK, "prisms");
  int sheet[7];
  expect(seams_monolayer(sys, 3.5, 0, 6, 100.0, sheet, 7) == SEAMS_OK, "monolayer");
  expect(sheet[6] > 0, "monolayer hexagons");
  seams_cage_set *again_cages = NULL;
  expect(seams_cages_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 4, 0, &again_cages) == SEAMS_OK,
         "cages find");
  seams_cages_destroy(again_cages);

  seams_system_destroy(sys);
  free(xyz);
  free(labels);
  free(q4);
  free(q6);
  free(sizes);
  free(sym);
  free(off_copy);
  free(col_copy);
  free(q3);
  free(q3bar);
  free(ql4);
  free(ql6);
  free(q6bar);
  free(verts);
  free(per_atom);
  free(shell);
  free(ice);
  free(z_mid);
  free(molality);
  free(ice_fraction);
  free(rdf_r);
  free(rdf_g);
  free(key);
  if (g_fails != 0) {
    fprintf(stderr, "%d checks failed\n", g_fails);
    return 1;
  }
  printf("seams_c_api ok %s\n", seams_version_string());
  return 0;
}
