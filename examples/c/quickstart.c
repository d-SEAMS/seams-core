#include "seams.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* 2x2x2 cubic ice, then the analyses libseams exposes.
 * An optional path reads frame 1 of that trajectory as well. */

static int fail(const char *step, const void *handle) {
  const char *detail = handle == NULL ? seams_last_error() : seams_handle_error(handle);
  fprintf(stderr, "%s: %s\n", step, detail);
  return 1;
}

int main(int argc, char **argv) {
  if (seams_version() != SEAMS_VERSION_NUMBER) {
    return fail("version", NULL);
  }
  /* 2x2x2 cubic ice. One cell is thinner than twice the 3.5 A cutoff,
   * so the cage census refuses it. */
  const double a = 2.75 * 4.0 / sqrt(3.0);
  const double fcc[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  const double basis[2][3] = {{0, 0, 0}, {0.25, 0.25, 0.25}};
  enum { NREP = 2, N = NREP * NREP * NREP * 8 };
  double xyz[N * 3];
  double box[9] = {0, 0, 0, NREP * a, NREP * a, NREP * a, 0, 0, 0};
  int n = 0;
  for (int i = 0; i < NREP; i++) {
    for (int j = 0; j < NREP; j++) {
      for (int k = 0; k < NREP; k++) {
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

  seams_system *sys = seams_system_create();
  if (sys == NULL) {
    return fail("create", NULL);
  }
  if (seams_system_set_atoms(sys, xyz, (size_t)n, box, 9, NULL, NULL) != SEAMS_OK) {
    return fail("frame", sys);
  }

  double q3[N], q4[N], q6[N], q6bar[N];
  if (seams_ql(sys, 3.5, 0, 3, SEAMS_GRAPH_CUTOFF, q3, (size_t)n, NULL, 0) != SEAMS_OK ||
      seams_ql(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, q4, (size_t)n, NULL, 0) != SEAMS_OK ||
      seams_ql(sys, 3.5, 0, 6, SEAMS_GRAPH_CUTOFF, q6, (size_t)n, q6bar, (size_t)n) != SEAMS_OK) {
    return fail("ql", sys);
  }
  for (int i = 0; i < n; i++) {
    if (!isfinite(q3[i]) || !isfinite(q4[i]) || !isfinite(q6[i])) {
      return fail("ql value", sys);
    }
  }
  double mean = 0, median = 0;
  int n_ice = 0;
  if (seams_q6bar(sys, 3.5, 0, 0.36, q6bar, (size_t)n, &mean, &median, &n_ice) != SEAMS_OK) {
    return fail("q6bar", sys);
  }

  seams_neighbours *nb = NULL;
  if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) != SEAMS_OK) {
    return fail("neighbours", sys);
  }
  size_t n_off = 0;
  size_t nnz = 0;
  if (seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz) != SEAMS_OK) {
    return fail("csr", nb);
  }
  seams_neighbours_destroy(nb);

  int labels[N], counts[9];
  if (seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, (size_t)n,
                  counts) != SEAMS_OK ||
      seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_OVITO, labels, (size_t)n, counts) !=
          SEAMS_OK ||
      seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, counts) !=
          SEAMS_OK) {
    return fail("chill", sys);
  }

  seams_rings *rings = NULL;
  if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK) {
    return fail("rings", sys);
  }
  size_t n_vert = 0;
  if (seams_rings_count(rings) > 0 &&
      seams_rings_vertices(rings, 0, NULL, 0, &n_vert) != SEAMS_OK) {
    return fail("vertices", rings);
  }
  seams_rings_destroy(rings);

  seams_cage_set *cages = NULL;
  if (seams_cages(sys, 3.5, 0, SEAMS_GRAPH_SEEDED, 0, 0, &cages) != SEAMS_OK) {
    return fail("cages", sys);
  }
  int n_hex = 0, n_cubic = 0, n_water = 0, n_both = 0;
  int per_atom[N];
  seams_cage_label tag = SEAMS_CAGE_WATER;
  if (seams_cages_counts(cages, &n_hex, &n_cubic, &n_water, &n_both) != SEAMS_OK ||
      seams_cages_labels(cages, per_atom, (size_t)n) != SEAMS_OK) {
    return fail("cage labels", cages);
  }
  tag = (seams_cage_label)per_atom[0];
  seams_cages_destroy(cages);
  size_t n_sig = 0, n_sig_atoms = 0;
  if (seams_signature_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", &n_sig, &n_sig_atoms) !=
      SEAMS_OK) {
    return fail("signature", sys);
  }

  int shell[N];
  if (seams_system_symmetry(sys, 0.0, shell, (size_t)n) != SEAMS_OK) {
    return fail("symmetry", sys);
  }

  int ice[N];
  for (int i = 0; i < n; i++) {
    ice[i] = 1;
  }
  seams_brine *front = NULL;
  if (seams_brine_analyse(sys, ice, (size_t)n, 0, 0, 4, 2, &front) != SEAMS_OK) {
    return fail("brine", sys);
  }
  double z_enter = 0, z_exit = 0, thickness = 0, c_ice = 0, c_liquid = 0, rejection = 0;
  int found = 0, defined = 0, n_pockets = 0, reservoir = 0;
  size_t n_bins = 0;
  if (seams_brine_fronts(front, &z_enter, &z_exit, &thickness, &found) != SEAMS_OK ||
      seams_brine_rejection(front, &c_ice, &c_liquid, &rejection, &defined) != SEAMS_OK ||
      seams_brine_profile(front, NULL, NULL, NULL, 0, &n_bins) != SEAMS_OK ||
      seams_brine_pockets(front, &n_pockets, &reservoir) != SEAMS_OK) {
    return fail("brine report", front);
  }
  seams_brine_destroy(front);

  if (argc > 1) {
    size_t nframes = 0;
    if (seams_path_nframes(argv[1], &nframes) != SEAMS_OK || nframes < 1) {
      return fail("nframes", NULL);
    }
    if (seams_system_read(sys, argv[1], 1, -1) != SEAMS_OK) {
      return fail("read", sys);
    }
  }

  printf("%s cubic %d q3 %.3f q6bar %.3f cages %d shell %d tag %d front %d\n",
         seams_version_string(), counts[SEAMS_LABEL_CUBIC], q3[0], mean,
         n_hex + n_cubic + n_both, shell[0], (int)tag, found);
  seams_system_destroy(sys);
  return 0;
}
