/* Copyright (c) d-SEAMS developers
 * Call every seams.h export on the in-tree fixtures and print one line each.
 * docs/check_c_api_tour.py requires those lines to appear in docs/api/c.md.
 */
#include "seams.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die(const char *what, seams_status status, const void *handle) {
  fprintf(stderr, "%s status %s last %s handle %s\n", what, seams_status_string(status),
          seams_last_error(), handle ? seams_handle_error(handle) : "");
  exit(1);
}

static void need(const char *what, seams_status status, const void *handle) {
  if (status != SEAMS_OK) {
    die(what, status, handle);
  }
}

static double mean_of(const double *values, int n) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += values[i];
  }
  return n > 0 ? sum / n : 0.0;
}

static int count_label(const int *labels, int n, int want) {
  int found = 0;
  for (int i = 0; i < n; ++i) {
    if (labels[i] == want) {
      ++found;
    }
  }
  return found;
}

static void diamond8(double *xyz, double *edge) {
  const double a = 2.75 * 4.0 / sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  int n = 0;
  for (int f = 0; f < 4; ++f) {
    for (int s = 0; s < 2; ++s, ++n) {
      xyz[3 * n] = (fcc[f][0] + basis[s][0]) * a;
      xyz[3 * n + 1] = (fcc[f][1] + basis[s][1]) * a;
      xyz[3 * n + 2] = (fcc[f][2] + basis[s][2]) * a;
    }
  }
  *edge = a;
}

int main(int argc, char **argv) {
  const char *root = argc > 1 ? argv[1] : ".";
  char ice[512], si[512], tip[512];
  snprintf(ice, sizeof ice, "%s/input/traj/ice1h.lammpstrj", root);
  snprintf(si, sizeof si, "%s/input/traj/genice_sI.lammpstrj", root);
  snprintf(tip, sizeof tip, "%s/input/traj/exampleTraj.lammpstrj", root);

  printf("capi seams_abi %d\n", seams_abi());
  printf("capi seams_version %d\n", seams_version());
  printf("capi seams_version_major %d\n", seams_version_major());
  printf("capi seams_version_minor %d\n", seams_version_minor());
  printf("capi seams_version_patch %d\n", seams_version_patch());
  printf("capi seams_version_string %s\n", seams_version_string());
  printf("capi seams_status_string %s\n", seams_status_string(SEAMS_OK));
  printf("capi seams_label_name %s\n", seams_label_name(SEAMS_LABEL_HEXAGONAL));
  printf("capi seams_ring_symmetry_name %s\n", seams_ring_symmetry_name(SEAMS_RING_HEXAGON));

  need("threads", seams_set_num_threads(1), NULL);
  printf("capi seams_set_num_threads %d\n", seams_get_num_threads());
  need("threads0", seams_set_num_threads(0), NULL);
  printf("capi seams_get_num_threads %d\n", seams_get_num_threads() > 0 ? 1 : 0);
  need("threads1", seams_set_num_threads(1), NULL);

  printf("capi seams_system_count %d\n", seams_system_count(NULL));
  seams_system *empty = seams_system_create();
  if (empty == NULL) {
    die("create", SEAMS_ERR_NOMEM, NULL);
  }
  printf("capi seams_system_create %d\n", empty != NULL);
  int scratch[8] = {0};
  seams_status bad = seams_chill_plus(empty, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                                       scratch, NULL);
  printf("capi seams_last_error %s\n", seams_last_error());
  printf("capi seams_last_error_message %s\n", seams_last_error_message());
  printf("capi seams_handle_error %s\n",
         strstr(seams_handle_error(empty), "no frame") ? "no-frame" : "other");
  seams_system_destroy(empty);
  printf("capi seams_system_destroy ok\n");
  seams_system_destroy(NULL);
  seams_neighbours_destroy(NULL);
  printf("capi seams_neighbours_destroy ok\n");
  seams_rings_destroy(NULL);
  printf("capi seams_rings_destroy ok\n");
  seams_cages_destroy(NULL);
  printf("capi seams_cages_destroy ok\n");
  seams_brine_destroy(NULL);
  printf("capi seams_brine_destroy ok\n");
  seams_fingerprint_destroy(NULL);
  printf("capi seams_fingerprint_destroy ok\n");

  double xyz[24], edge = 0;
  diamond8(xyz, &edge);
  double box9[9] = {0, 0, 0, edge, edge, edge, 0, 0, 0};
  seams_system *cell = seams_system_create();
  need("set_atoms", seams_system_set_atoms(cell, xyz, 8, box9, 9, NULL, NULL), cell);
  printf("capi seams_system_set_atoms %d\n", seams_system_count(cell));
  need("set_frame", seams_system_set_frame(cell, xyz, 8, box9, NULL, NULL), cell);
  printf("capi seams_system_set_frame %d\n", seams_system_count(cell));
  need("set_ortho",
       seams_system_set_ortho(cell, xyz, 8, 0, 0, 0, edge, edge, edge, NULL, NULL), cell);
  printf("capi seams_system_set_ortho %d\n", seams_system_count(cell));
  int labels8[8];
  int counts9[9];
  need("chill8",
       seams_chill_plus(cell, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels8, counts9),
       cell);
  printf("capi seams_chill_plus cubic %d\n", counts9[SEAMS_LABEL_CUBIC]);

  size_t nframes = 0;
  need("nframes", seams_path_nframes(ice, &nframes), NULL);
  printf("capi seams_path_nframes %zu\n", nframes);
  seams_system *ih = seams_system_create();
  need("read", seams_system_read(ih, ice, 1, 1), ih);
  const int n = seams_system_count(ih);
  printf("capi seams_system_read %d\n", n);
  size_t natoms = 0;
  need("natoms", seams_system_natoms(ih, &natoms), ih);
  printf("capi seams_system_natoms %zu\n", natoms);

  double moved[24];
  for (int k = 0; k < 24; ++k) {
    moved[k] = xyz[k] + 0.25;
  }
  need("set_positions", seams_system_set_positions(cell, moved, 8, NULL, 0), cell);
  printf("capi seams_system_set_positions %d\n", seams_system_count(cell));
  need("stage_positions", seams_system_stage_positions(cell, 1), cell);
  need("set_positions", seams_system_set_positions(cell, moved, 8, NULL, 0), cell);
  need("stage_positions", seams_system_stage_positions(cell, 0), cell);
  printf("capi seams_system_stage_positions %d\n", seams_system_count(cell));
  double cell9[9];
  need("box", seams_system_box(cell, cell9, 9), cell);
  printf("capi seams_system_box %.3f %.3f %.3f %.3f %.3f %.3f\n", cell9[0], cell9[1], cell9[2],
         cell9[3], cell9[4], cell9[5]);
  double back[24];
  need("positions", seams_system_positions(cell, back, 24), cell);
  printf("capi seams_system_positions %.3f %.3f %.3f\n", back[0], back[1], back[2]);
  int types8[8];
  need("types", seams_system_types(cell, types8, 8), cell);
  int ones = 0;
  for (int i = 0; i < 8; ++i) {
    ones += types8[i] == 1;
  }
  printf("capi seams_system_types %d of type 1\n", ones);
  int ids8[8];
  need("ids", seams_system_ids(cell, ids8, 8), cell);
  printf("capi seams_system_ids %d %d\n", ids8[0], ids8[7]);

  seams_neighbours *nb = NULL;
  need("neigh", seams_neighbours_build(ih, 3.5, 1, SEAMS_GRAPH_CUTOFF, &nb), ih);
  printf("capi seams_neighbours_build %d\n", seams_neighbours_count(nb));
  printf("capi seams_neighbours_count %d\n", seams_neighbours_count(nb));
  size_t rows = 0;
  need("nrows", seams_neighbours_nrows(nb, &rows), nb);
  printf("capi seams_neighbours_nrows %zu\n", rows);
  size_t n_off = 0, nnz = 0;
  need("copyq", seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz), nb);
  int *offsets = calloc(n_off, sizeof(int));
  int *cols = calloc(nnz ? nnz : 1, sizeof(int));
  need("copy", seams_neighbours_copy(nb, offsets, n_off, &n_off, cols, nnz, &nnz), nb);
  printf("capi seams_neighbours_copy %zu %zu\n", n_off, nnz);
  const int *boff = NULL, *bcol = NULL;
  int bn = 0, bnnz = 0;
  need("csr", seams_neighbours_csr(nb, &boff, &bn, &bcol, &bnnz), nb);
  printf("capi seams_neighbours_csr %d %d\n", bn, bnnz);

  seams_neighbours *knn = NULL;
  need("knn", seams_neighbours_knn(ih, 4, 3.5, 1, 1, &knn), ih);
  printf("capi seams_neighbours_knn %d\n", seams_neighbours_count(knn));

  int *labels = calloc((size_t)n, sizeof(int));
  need("chill",
       seams_chill(ih, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, (size_t)n,
                   counts9),
       ih);
  printf("capi seams_chill hexagonal %d\n", counts9[SEAMS_LABEL_HEXAGONAL]);
  need("plus",
       seams_chill_plus(ih, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts9), ih);
  printf("capi seams_chill_plus hexagonal %d\n", counts9[SEAMS_LABEL_HEXAGONAL]);
  need("plus_n",
       seams_chill_plus_n(ih, 3.5, 1, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, natoms,
                          counts9),
       ih);
  printf("capi seams_chill_plus_n hexagonal %d\n", counts9[SEAMS_LABEL_HEXAGONAL]);

  double *q4 = calloc((size_t)n, sizeof(double));
  double *q6 = calloc((size_t)n, sizeof(double));
  double *q6bar = calloc((size_t)n, sizeof(double));
  need("q4q6", seams_q4q6(ih, 3.5, 1, q4, q6), ih);
  printf("capi seams_q4q6 %.6f\n", mean_of(q6, n));
  need("q4q6_n", seams_q4q6_n(ih, 3.5, 1, q4, natoms, q6, natoms), ih);
  printf("capi seams_q4q6_n %.6f\n", mean_of(q6, n));
  need("ql", seams_ql(ih, 3.5, 1, 6, SEAMS_GRAPH_CUTOFF, q6, (size_t)n, q6bar, (size_t)n), ih);
  printf("capi seams_ql %.6f\n", mean_of(q6bar, n));
  double mean = 0, median = 0;
  int n_ice = 0;
  need("q6bar", seams_q6bar(ih, 3.5, 1, 0.36, q6bar, (size_t)n, &mean, &median, &n_ice), ih);
  printf("capi seams_q6bar %.6f %d\n", mean, n_ice);
  need("stein", seams_steinhardt(ih, 3.5, 1, 6, q6, (size_t)n, q6bar, (size_t)n), ih);
  printf("capi seams_steinhardt %.6f\n", mean_of(q6bar, n));

  seams_rings *rings = NULL;
  need("rings", seams_rings_find(ih, 3.5, 1, SEAMS_GRAPH_CUTOFF, 6, &rings), ih);
  printf("capi seams_rings_find %d\n", seams_rings_count(rings));
  printf("capi seams_rings_count %d\n", seams_rings_count(rings));
  size_t nrings = 0;
  need("nrings", seams_rings_nrings(rings, &nrings), rings);
  printf("capi seams_rings_nrings %zu\n", nrings);
  int hist[13] = {0};
  need("hist", seams_rings_histogram(rings, hist, 13), rings);
  printf("capi seams_rings_histogram %d\n", hist[6]);
  need("hist_n", seams_rings_histogram_n(rings, hist, 13), rings);
  printf("capi seams_rings_histogram_n %d\n", hist[6]);
  int nr = seams_rings_count(rings);
  int *sizes = calloc((size_t)(nr > 0 ? nr : 1), sizeof(int));
  int *sym = calloc((size_t)(nr > 0 ? nr : 1), sizeof(int));
  need("rlab", seams_rings_labels(rings, sizes, sym), rings);
  printf("capi seams_rings_labels %d\n", nr > 0 ? sizes[0] : 0);
  need("rlab_n", seams_rings_labels_n(rings, sizes, nrings, sym, nrings), rings);
  printf("capi seams_rings_labels_n %d\n", nr > 0 ? sizes[0] : 0);
  size_t nv = 0;
  need("vq", seams_rings_vertices(rings, 0, NULL, 0, &nv), rings);
  int *verts = calloc(nv ? nv : 1, sizeof(int));
  need("vv", seams_rings_vertices(rings, 0, verts, nv, &nv), rings);
  printf("capi seams_rings_vertices %zu\n", nv);
  double fraction = 0;
  need("frac", seams_rings_fraction(rings, 6, &fraction), rings);
  printf("capi seams_rings_fraction %.6f\n", fraction);

  seams_cage_set *cages = NULL;
  need("cfind", seams_cages_find(ih, 3.5, 1, SEAMS_GRAPH_CUTOFF, 4, 0, &cages), ih);
  int n_hex = 0, n_cubic = 0, n_water = 0, n_both = 0;
  need("cc", seams_cages_counts(cages, &n_hex, &n_cubic, &n_water, &n_both), cages);
  printf("capi seams_cages_find %d %d %d %d\n", n_hex, n_cubic, n_water, n_both);
  printf("capi seams_cages_counts %d\n", n_hex + n_cubic + n_water + n_both);
  int *per = calloc((size_t)n, sizeof(int));
  /* per[] holds seams_cage_label values, one per atom. */
  need("clab", seams_cages_labels(cages, per, (size_t)n), cages);
  printf("capi seams_cages_labels %d\n", count_label(per, n, SEAMS_CAGE_WATER));
  seams_cages_destroy(cages);
  need("cages", seams_cages(ih, 3.5, 1, SEAMS_GRAPH_CUTOFF, 4, 0, &cages), ih);
  need("cc2", seams_cages_counts(cages, &n_hex, &n_cubic, &n_water, &n_both), cages);
  printf("capi seams_cages %d\n", n_water);
  seams_cages_destroy(cages);
  printf("capi seams_cages_destroy ok\n");

  seams_system *guest = seams_system_create();
  need("si", seams_system_read(guest, si, 1, 1), guest);
  size_t n_cages = 0, n_atoms = 0;
  need("sig",
       seams_cages_signature(guest, 3.5, 1, SEAMS_GRAPH_CUTOFF, "512", &n_cages, &n_atoms), guest);
  printf("capi seams_cages_signature %zu %zu\n", n_cages, n_atoms);
  n_cages = 0;
  n_atoms = 0;
  need("sig2", seams_signature_cages(guest, 3.5, 1, SEAMS_GRAPH_CUTOFF, "512", &n_cages, &n_atoms),
       guest);
  printf("capi seams_signature_cages %zu %zu\n", n_cages, n_atoms);

  int *shell = calloc((size_t)n, sizeof(int));
  need("sym", seams_system_symmetry(ih, 0, shell, (size_t)n), ih);
  printf("capi seams_system_symmetry %d\n", count_label(shell, n, 0));

  int *ice_flag = calloc((size_t)n, sizeof(int));
  for (int i = 0; i < n / 2; ++i) {
    ice_flag[i] = 1;
  }
  seams_brine *front = NULL;
  need("brine", seams_brine_analyse(ih, ice_flag, (size_t)n, 1, 0, 0, 2, &front), ih);
  printf("capi seams_brine_analyse ok\n");
  double z0 = 0, z1 = 0, thick = 0;
  int found = 0;
  need("fronts", seams_brine_fronts(front, &z0, &z1, &thick, &found), front);
  printf("capi seams_brine_fronts %d\n", found);
  double c_ice = 0, c_liq = 0, kval = 0;
  int defined = 0;
  need("rej", seams_brine_rejection(front, &c_ice, &c_liq, &kval, &defined), front);
  printf("capi seams_brine_rejection %d\n", defined);
  size_t nbins = 0;
  need("profq", seams_brine_profile(front, NULL, NULL, NULL, 0, &nbins), front);
  printf("capi seams_brine_profile %zu\n", nbins);
  int pockets = 0, reservoir = 0;
  need("pock", seams_brine_pockets(front, &pockets, &reservoir), front);
  printf("capi seams_brine_pockets %d %d\n", pockets, reservoir);
  seams_brine_destroy(front);

  size_t rdf_n = 0;
  double volume = 0;
  need("rdfq", seams_rdf(ih, 6.0, 20, 1, 1, NULL, NULL, NULL, 0, &rdf_n, &volume), ih);
  double *rr = calloc(rdf_n ? rdf_n : 1, sizeof(double));
  double *gg = calloc(rdf_n ? rdf_n : 1, sizeof(double));
  int *rc = calloc(rdf_n ? rdf_n : 1, sizeof(int));
  need("rdf", seams_rdf(ih, 6.0, 20, 1, 1, rr, gg, rc, rdf_n, &rdf_n, &volume), ih);
  printf("capi seams_rdf %zu %.3f\n", rdf_n, volume);
  size_t cn_n = 0;
  need("rcnq", seams_running_cn(ih, 6.0, 20, 1, 1, NULL, NULL, 0, &cn_n), ih);
  double *cnr = calloc(cn_n ? cn_n : 1, sizeof(double));
  double *cnv = calloc(cn_n ? cn_n : 1, sizeof(double));
  need("rcn", seams_running_cn(ih, 6.0, 20, 1, 1, cnr, cnv, cn_n, &cn_n), ih);
  printf("capi seams_running_cn %.3f\n", cn_n ? cnv[cn_n - 1] : 0);
  double coord = 0;
  need("coord", seams_coordination(ih, 3.5, 20, 1, 1, &coord), ih);
  printf("capi seams_coordination %.3f\n", coord);

  seams_system *mix = seams_system_create();
  need("tip", seams_system_read(mix, tip, 1, -1), mix);
  int bonds = 0;
  need("hb", seams_hbonds(mix, 3.5, 2, 1, 2.42, 30.0, &bonds), mix);
  printf("capi seams_hbonds %d\n", bonds);
  size_t dens_n = 0;
  need("densq", seams_density_axis(mix, 2, 4, 2, NULL, NULL, 0, &dens_n), mix);
  double *zz = calloc(dens_n ? dens_n : 1, sizeof(double));
  double *rho = calloc(dens_n ? dens_n : 1, sizeof(double));
  need("dens", seams_density_axis(mix, 2, 4, 2, zz, rho, dens_n, &dens_n), mix);
  printf("capi seams_density_axis %zu\n", dens_n);

  seams_topo *fp = NULL;
  need("fp", seams_fingerprint(ih, 3.5, 1, 4, SEAMS_GRAPH_CUTOFF, 2, 0, &fp), ih);
  printf("capi seams_fingerprint ok\n");
  size_t key_n = 0;
  need("keyq", seams_fingerprint_key(fp, NULL, 0, &key_n), fp);
  char *key = calloc(key_n + 1, 1);
  need("key", seams_fingerprint_key(fp, key, key_n + 1, &key_n), fp);
  printf("capi seams_fingerprint_key %zu\n", key_n);
  size_t atom_n = 0;
  need("atomq", seams_fingerprint_atom(fp, 0, NULL, 0, &atom_n), fp);
  printf("capi seams_fingerprint_atom %zu\n", atom_n);
  int fhist[13] = {0};
  need("fpr", seams_fingerprint_rings(fp, fhist, 13), fp);
  printf("capi seams_fingerprint_rings %d\n", fhist[6]);
  seams_fingerprint_destroy(fp);

  size_t pairs = 0;
  int ncat = 0, nanion = 0;
  need("pairs", seams_contact_pairs(mix, "1=cationHead,2=anion", &pairs, &ncat, &nanion), mix);
  printf("capi seams_contact_pairs %zu %d\n", pairs, nanion);
  int nsub = 0, largest = 0;
  double perc = 0;
  need("dom", seams_domains(mix, 3.5, "1=cationHead,2=anion", 1, &nsub, &largest, &perc), mix);
  printf("capi seams_domains %d %d\n", nsub, largest);

  int perfect[8] = {0}, imperfect[8] = {0};
  need("prisms", seams_prisms(ih, 3.5, 1, perfect, imperfect, 8), ih);
  printf("capi seams_prisms %d %d\n", perfect[6], imperfect[6]);
  int mono[9] = {0};
  need("mono", seams_monolayer(ih, 3.5, 1, 6, 100.0, mono, 9), ih);
  printf("capi seams_monolayer %d\n", mono[6]);
  int *kinds = calloc((size_t)seams_system_count(mix), sizeof(int));
  /* kinds[] holds seams_site_kind values, one per atom. */
  need("sites", seams_site_table(mix, "1=waterH,2=waterO", kinds, (size_t)seams_system_count(mix)),
       mix);
  printf("capi seams_site_table %d\n", count_label(kinds, seams_system_count(mix), SEAMS_SITE_WATER_O));

  seams_rings_destroy(rings);
  seams_neighbours_destroy(nb);
  seams_neighbours_destroy(knn);
  seams_system_destroy(cell);
  seams_system_destroy(ih);
  seams_system_destroy(guest);
  seams_system_destroy(mix);
  free(offsets);
  free(cols);
  free(labels);
  free(q4);
  free(q6);
  free(q6bar);
  free(sizes);
  free(sym);
  free(verts);
  free(per);
  free(shell);
  free(ice_flag);
  free(rr);
  free(gg);
  free(rc);
  free(cnr);
  free(cnv);
  free(zz);
  free(rho);
  free(key);
  free(kinds);
  return 0;
}
