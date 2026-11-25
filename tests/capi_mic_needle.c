/* Copyright (c) d-SEAMS developers
 * Minimum image on needle cells through the C API.
 *
 * A needle has an a edge far longer than b and c. Here b and c are tilted
 * so that the fractional wrap misses the nearest image inside the cutoff,
 * and the cutoff is past the Smith ball radius, so every entry point below
 * needs the refined minimum image. Each runs at one thread and at eight.
 *
 * With lx from 2^16 to 2^50 every result must be the exact minimum image.
 * The cell numbers and the coordinates are short dyadic fractions, so the
 * brute-force references here are exact in binary64, and the atoms and
 * references are the same for every such lx. q4, q6 and CHILL+ are
 * compared with the same atoms in a large orthorhombic box, with the
 * images that the needle's minimum image reaches copied in explicitly.
 * The mutual graph may instead report SEAMS_ERR_RANGE when the cell needs
 * more linked-cell bins than linkcell allows.
 *
 * With lx = 2^57, 2e16, 1e17, 1e160 and 2^532 (ly 1, lz 0.4 and tilts 0.3
 * for the decimal ones) the lattice coefficients are past binary64, and
 * each entry point must return SEAMS_ERR_RANGE: never SEAMS_OK with wrong
 * pairs, never SEAMS_ERR_INTERNAL and never an abort.
 */
#include <seams.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { N_ATOMS = 640, N_HALF = N_ATOMS / 2, WIDE = 8, MAX_GHOSTS = 16384, KNN = 4 };

typedef struct {
  const char *name;
  double lx, ly, lz, xy, xz, yz;
} cell_t;

#define DYADIC(name, lx) {name, lx, 1.0, 0.375, 0.3125, 0.3125, 0.3125}
#define DECIMAL(name, lx) {name, lx, 1.0, 0.4, 0.3, 0.3, 0.3}

static const cell_t kInRange[] = {
    DYADIC("lx 2^16", 65536.0),
    DYADIC("lx 2^30", 1073741824.0),
    DYADIC("lx 2^40", 1099511627776.0),
    DYADIC("lx 2^50", 1125899906842624.0),
};
enum { N_IN_RANGE = sizeof kInRange / sizeof kInRange[0] };

/* ldexp(1, 57) and ldexp(1, 532). */
static const cell_t kPast[] = {
    DYADIC("lx 2^57", 144115188075855872.0),
    DECIMAL("lx 2e16", 2e16),
    DECIMAL("lx 1e17", 1e17),
    DECIMAL("lx 1e160", 1e160),
    DYADIC("lx 2^532", 0x1p532),
};
enum { N_PAST = sizeof kPast / sizeof kPast[0] };

/* The dyadic cells' shortest lattice vector is c, of length 0.58, so at
 * most one image of an atom is inside 0.27, and their Smith ball radius is
 * lz / 2 = 0.1875. */
static const double kCutoff = 0.27;
static const double kGateCutoff = 0.354858;

static int failures = 0;

static void fail(const char *fmt, const char *where, const char *what, long a, long b) {
  printf("FAIL %s %s: ", where, what);
  printf(fmt, a, b);
  printf("\n");
  failures++;
}

/* Exact squared distances lie on a 2^-28 grid. This square is half way
 * between two grid points, so no pair is within rounding of the cutoff. */
static double off_grid_square(double r) {
  return (floor(r * r * 268435456.0) + 0.5) / 268435456.0;
}

static uint64_t rng_state;

static double draw(void) {
  rng_state ^= rng_state << 13;
  rng_state ^= rng_state >> 7;
  rng_state ^= rng_state << 17;
  return (double)(rng_state >> 54) / 1024.0;
}

/* fb b + fc c + xa along x, with fb, fc and xa / 16 on a 2^-10 grid: the
 * atoms fill the b-c cross-section over a strip 16 long. Types alternate
 * and each atom is its own molecule. */
static void make_atoms(const cell_t *cell, double *xyz, int *types, int *mols) {
  rng_state = 0x9e3779b97f4a7c15ull;
  for (int i = 0; i < N_ATOMS; i++) {
    const double fb = draw();
    const double fc = draw();
    const double xa = 16.0 * draw();
    xyz[3 * i] = xa + fb * cell->xy + fc * cell->xz;
    xyz[3 * i + 1] = fb * cell->ly + fc * cell->yz;
    xyz[3 * i + 2] = fc * cell->lz;
    types[i] = 1 + i % 2;
    mols[i] = i + 1;
  }
}

/* Smallest squared distance over the images of q within `bound` of p, or
 * INFINITY when there is none. For the dyadic cells every coordinate is a
 * multiple of 2^-14 below 2^5 and every cell number but lx a multiple of
 * 2^-4, and the coefficients below stay under 2^4 except along a, where
 * any nonzero coefficient lands far outside the bound. Every difference,
 * product and sum that can reach the minimum is therefore exact. The
 * coefficient ranges are widened by one each way. */
static double exact_d2(const cell_t *c, const double *p, const double *q, double bound) {
  const double dx = p[0] - q[0];
  const double dy = p[1] - q[1];
  const double dz = p[2] - q[2];
  double best = INFINITY;
  const double nc_lo = floor((dz - bound) / c->lz) - 1.0;
  const double nc_hi = ceil((dz + bound) / c->lz) + 1.0;
  for (double nc = nc_lo; nc <= nc_hi; nc += 1.0) {
    const double ry = dy - nc * c->yz;
    const double nb_lo = floor((ry - bound) / c->ly) - 1.0;
    const double nb_hi = ceil((ry + bound) / c->ly) + 1.0;
    for (double nb = nb_lo; nb <= nb_hi; nb += 1.0) {
      const double rx = dx - nb * c->xy - nc * c->xz;
      const double na_lo = floor((rx - bound) / c->lx) - 1.0;
      const double na_hi = ceil((rx + bound) / c->lx) + 1.0;
      for (double na = na_lo; na <= na_hi; na += 1.0) {
        const double vx = rx - na * c->lx;
        const double vy = ry - nb * c->ly;
        const double vz = dz - nc * c->lz;
        const double d2 = vx * vx + vy * vy + vz * vz;
        if (d2 < best && d2 < bound * bound) {
          best = d2;
        }
      }
    }
  }
  return best;
}

/* The image the fractional wrap picks, evaluated exactly. */
static double wrap_d2(const cell_t *c, const double *p, const double *q) {
  const double dx = p[0] - q[0];
  const double dy = p[1] - q[1];
  const double dz = p[2] - q[2];
  const double sc = dz / c->lz;
  const double sb = (dy - sc * c->yz) / c->ly;
  const double sa = (dx - sb * c->xy - sc * c->xz) / c->lx;
  const double na = nearbyint(sa);
  const double nb = nearbyint(sb);
  const double nc = nearbyint(sc);
  const double vx = dx - na * c->lx - nb * c->xy - nc * c->xz;
  const double vy = dy - nb * c->ly - nc * c->yz;
  const double vz = dz - nc * c->lz;
  return vx * vx + vy * vy + vz * vz;
}

static seams_system *make_system(const cell_t *c, const double *xyz, const int *types,
                                 const int *mols, size_t n) {
  const double box9[9] = {0, 0, 0, c->lx, c->ly, c->lz, c->xy, c->xz, c->yz};
  seams_system *sys = seams_system_create();
  const seams_status st = seams_system_set_atoms(sys, xyz, n, box9, 9, types, mols);
  if (sys == NULL || st != SEAMS_OK) {
    printf("FAIL %s set_atoms status %d: %s\n", c->name, (int)st, seams_last_error_message());
    exit(1);
  }
  return sys;
}

typedef struct {
  unsigned char *pair;  /* N_ATOMS x N_ATOMS cutoff flags */
  unsigned char *knn;   /* N_ATOMS x N_ATOMS mutual four-nearest flags */
  long n_pairs;
  long n_knn;
  long wrap_misses;     /* cutoff pairs whose wrapped image is outside */
  int nearest[N_HALF];  /* per guest, in input order */
  double nearest_d2[N_HALF];
  int contacts[N_HALF];
  long mutual;
} exact_ref;

/* The unique nearest atom of type `want` to atom i, or -1 on a tie. */
static int nearest_of(const cell_t *c, const double *xyz, const int *types, int i, int want,
                      double *d2_out) {
  double best = INFINITY;
  double second = INFINITY;
  int arg = -1;
  for (int j = 0; j < N_ATOMS; j++) {
    if (j == i || types[j] != want) {
      continue;
    }
    const double d2 = exact_d2(c, &xyz[3 * i], &xyz[3 * j], 1.0);
    if (d2 < best) {
      second = best;
      best = d2;
      arg = j;
    } else if (d2 < second) {
      second = d2;
    }
  }
  *d2_out = best;
  return isfinite(best) && best < second ? arg : -1;
}

/* The KNN nearest atoms of i, of any type, when the KNN-th is strictly
 * nearer than the next; otherwise 0. */
static int nearest_k(const cell_t *c, const double *xyz, int i, int *out) {
  double d2[KNN + 1];
  int arg[KNN + 1];
  for (int m = 0; m <= KNN; m++) {
    d2[m] = INFINITY;
    arg[m] = -1;
  }
  for (int j = 0; j < N_ATOMS; j++) {
    if (j == i) {
      continue;
    }
    const double v = exact_d2(c, &xyz[3 * i], &xyz[3 * j], 1.0);
    if (!(v < d2[KNN])) {
      continue;
    }
    int m = KNN;
    while (m > 0 && v < d2[m - 1]) {
      d2[m] = d2[m - 1];
      arg[m] = arg[m - 1];
      m--;
    }
    d2[m] = v;
    arg[m] = j;
  }
  memcpy(out, arg, KNN * sizeof(int));
  return isfinite(d2[KNN - 1]) && d2[KNN - 1] < d2[KNN];
}

static void build_reference(const cell_t *c, const double *xyz, const int *types, double rc2,
                            exact_ref *ref) {
  memset(ref->pair, 0, (size_t)N_ATOMS * N_ATOMS);
  memset(ref->knn, 0, (size_t)N_ATOMS * N_ATOMS);
  ref->n_pairs = 0;
  ref->wrap_misses = 0;
  for (int i = 0; i < N_ATOMS; i++) {
    for (int j = 0; j < N_ATOMS; j++) {
      if (i == j) {
        continue;
      }
      const double d2 = exact_d2(c, &xyz[3 * i], &xyz[3 * j], 1.0);
      if (d2 == 0.0) {
        printf("FAIL test data: atoms %d and %d coincide\n", i, j);
        exit(1);
      }
      if (d2 < rc2) {
        ref->pair[(size_t)i * N_ATOMS + j] = 1;
        ref->n_pairs++;
        ref->wrap_misses += wrap_d2(c, &xyz[3 * i], &xyz[3 * j]) >= rc2;
      }
    }
  }
  static int near4[N_ATOMS][KNN];
  for (int i = 0; i < N_ATOMS; i++) {
    if (!nearest_k(c, xyz, i, near4[i])) {
      printf("FAIL test data: atom %d has no unique %d nearest atoms within 1\n", i, KNN);
      exit(1);
    }
  }
  ref->n_knn = 0;
  for (int i = 0; i < N_ATOMS; i++) {
    for (int m = 0; m < KNN; m++) {
      const int j = near4[i][m];
      for (int back = 0; back < KNN; back++) {
        if (near4[j][back] == i) {
          ref->knn[(size_t)i * N_ATOMS + j] = 1;
          ref->n_knn++;
        }
      }
    }
  }
  int partner[N_ATOMS];
  for (int i = 0; i < N_ATOMS; i++) {
    double d2 = 0.0;
    partner[i] = nearest_of(c, xyz, types, i, types[i] == 1 ? 2 : 1, &d2);
    if (partner[i] < 0) {
      printf("FAIL test data: atom %d has no unique nearest unlike atom within 1\n", i);
      exit(1);
    }
    if (types[i] == 2) {
      ref->nearest[i / 2] = partner[i];
      ref->nearest_d2[i / 2] = d2;
      int count = 0;
      for (int j = 0; j < N_ATOMS; j++) {
        count += types[j] == 1 && ref->pair[(size_t)i * N_ATOMS + j];
      }
      ref->contacts[i / 2] = count;
    }
  }
  ref->mutual = 0;
  for (int i = 0; i < N_ATOMS; i++) {
    ref->mutual += types[i] == 1 && partner[partner[i]] == i;
  }
}

/* The atoms plus every image within twice the cutoff of one of them, in an
 * orthorhombic box wide enough that its own wrap is exact. CHILL+ reads
 * the neighbours' neighbours, hence twice the cutoff. Returns the count. */
static size_t explicit_images(const cell_t *c, const double *xyz, double reach, double *out,
                              int *types_out) {
  size_t n = 0;
  for (int i = 0; i < N_ATOMS; i++) {
    memcpy(&out[3 * n], &xyz[3 * i], 3 * sizeof(double));
    types_out[n++] = 1 + i % 2;
  }
  for (int nb = -3; nb <= 3; nb++) {
    for (int nc = -4; nc <= 4; nc++) {
      if (nb == 0 && nc == 0) {
        continue;
      }
      for (int j = 0; j < N_ATOMS; j++) {
        const double g[3] = {xyz[3 * j] + nb * c->xy + nc * c->xz,
                             xyz[3 * j + 1] + nb * c->ly + nc * c->yz, xyz[3 * j + 2] + nc * c->lz};
        int near = 0;
        for (int i = 0; i < N_ATOMS && !near; i++) {
          const double vx = g[0] - xyz[3 * i];
          const double vy = g[1] - xyz[3 * i + 1];
          const double vz = g[2] - xyz[3 * i + 2];
          near = vx * vx + vy * vy + vz * vz < reach * reach;
        }
        if (!near) {
          continue;
        }
        if (nb == -3 || nb == 3 || nc == -4 || nc == 4 || n == MAX_GHOSTS) {
          printf("FAIL test data: images reach translation (%d, %d)\n", nb, nc);
          exit(1);
        }
        memcpy(&out[3 * n], g, sizeof g);
        types_out[n++] = 1 + j % 2;
      }
    }
  }
  return n;
}

typedef struct {
  int *offsets;
  int *cols;
  size_t nnz;
  int *knn_offsets;
  int *knn_cols;
  size_t knn_nnz;
  int knn_built;
  double q4[N_ATOMS];
  double q6[N_ATOMS];
  int labels[N_ATOMS];
  int counts9[9];
  int guest[N_HALF];
  int nearest[N_HALF];
  int nearest_type[N_HALF];
  double nearest_distance[N_HALF];
  int contacts[N_HALF];
  size_t mutual;
  int n_cation;
  int n_anion;
} results;

static void free_results(results *r) {
  free(r->offsets);
  free(r->cols);
  free(r->knn_offsets);
  free(r->knn_cols);
  r->offsets = NULL;
  r->cols = NULL;
  r->knn_offsets = NULL;
  r->knn_cols = NULL;
}

static seams_status neighbours(seams_system *sys, double rc, seams_graph graph, int **offsets,
                               int **cols, size_t *nnz) {
  seams_neighbours *nb = NULL;
  const seams_status st = seams_neighbours_build(sys, rc, 0, graph, &nb);
  if (st != SEAMS_OK) {
    if (nb != NULL) {
      printf("FAIL neighbours: a failed build returned a handle\n");
      failures++;
      seams_neighbours_destroy(nb);
    }
    return st;
  }
  size_t n_off = 0;
  seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, nnz);
  *offsets = malloc((n_off + 1) * sizeof(int));
  *cols = malloc((*nnz + 1) * sizeof(int));
  if (*offsets == NULL || *cols == NULL || n_off != N_ATOMS + 1 ||
      seams_neighbours_copy(nb, *offsets, n_off, &n_off, *cols, *nnz, nnz) != SEAMS_OK) {
    printf("FAIL neighbours: copy\n");
    exit(1);
  }
  seams_neighbours_destroy(nb);
  return SEAMS_OK;
}

static void compare_pairs(const char *where, const char *what, const int *offsets,
                          const int *cols, const unsigned char *want) {
  long missing = 0;
  long extra = 0;
  static unsigned char seen[N_ATOMS];
  for (int i = 0; i < N_ATOMS; i++) {
    memset(seen, 0, sizeof seen);
    for (int k = offsets[i]; k < offsets[i + 1]; k++) {
      const int j = cols[k];
      if (j < 0 || j >= N_ATOMS || seen[j]) {
        extra++;
        continue;
      }
      seen[j] = 1;
      extra += !want[(size_t)i * N_ATOMS + j];
    }
    for (int j = 0; j < N_ATOMS; j++) {
      missing += want[(size_t)i * N_ATOMS + j] && !seen[j];
    }
  }
  if (missing != 0 || extra != 0) {
    fail("%ld missing and %ld extra directed pairs", where, what, missing, extra);
  }
}

static seams_status host_guest(seams_system *sys, double rc, results *r) {
  size_t n_guest = 0;
  size_t n_types = 0;
  int framework_types[4];
  const seams_status st = seams_host_guest(
      sys, 2, 1, rc, r->guest, r->nearest, r->nearest_type, r->nearest_distance, N_HALF,
      &n_guest, framework_types, 4, &n_types, r->contacts, N_HALF);
  if (st == SEAMS_OK && (n_guest != N_HALF || n_types != 1 || framework_types[0] != 1)) {
    printf("FAIL host_guest: %zu guests and %zu framework types\n", n_guest, n_types);
    exit(1);
  }
  return st;
}

static seams_status contact_pairs(seams_system *sys, results *r) {
  return seams_contact_pairs(sys, "1=cationHead,2=anion", &r->mutual, &r->n_cation, &r->n_anion);
}

/* SEAMS_ERR_RANGE whose message starts with `prefix`. */
static int is_range(seams_status st, const char *prefix) {
  const char *msg = seams_last_error_message();
  return st == SEAMS_ERR_RANGE && msg != NULL && strncmp(msg, prefix, strlen(prefix)) == 0;
}

static void expect_range(const char *where, const char *what, seams_status st,
                         const char *prefix) {
  if (!is_range(st, prefix)) {
    const char *msg = seams_last_error_message();
    printf("FAIL %s %s: status %d (%s), message '%s'; expected SEAMS_ERR_RANGE from '%s'\n",
           where, what, (int)st, seams_status_string(st), msg == NULL ? "" : msg, prefix);
    failures++;
  }
}

static const char kMinimumImage[] = "minimum-image";
static const char kMesh[] = "linked-cell mesh";
static const char kSingular[] = "singular or non-positive cell";

static void run_exact(const cell_t *c, seams_system *sys, double rc, const exact_ref *ref,
                      const results *expl, results *r, int width) {
  char where[64];
  snprintf(where, sizeof where, "%s threads %d", c->name, width);
  seams_status st = neighbours(sys, rc, SEAMS_GRAPH_CUTOFF, &r->offsets, &r->cols, &r->nnz);
  if (st != SEAMS_OK) {
    fail("status %ld, expected %ld", where, "neighbours", st, SEAMS_OK);
  } else {
    compare_pairs(where, "neighbours", r->offsets, r->cols, ref->pair);
  }

  st = neighbours(sys, rc, SEAMS_GRAPH_MUTUAL, &r->knn_offsets, &r->knn_cols, &r->knn_nnz);
  r->knn_built = st == SEAMS_OK;
  if (st == SEAMS_OK) {
    compare_pairs(where, "mutual neighbours", r->knn_offsets, r->knn_cols, ref->knn);
  } else if (is_range(st, kMesh)) {
    printf("%s: mutual neighbours, SEAMS_ERR_RANGE: %s\n", where, seams_last_error_message());
  } else {
    expect_range(where, "mutual neighbours", st, kMesh);
  }

  st = seams_q4q6_n(sys, rc, 0, r->q4, N_ATOMS, r->q6, N_ATOMS);
  if (st != SEAMS_OK) {
    fail("status %ld, expected %ld", where, "q4q6", st, SEAMS_OK);
  } else {
    for (int i = 0; i < N_ATOMS; i++) {
      if (!(fabs(r->q4[i] - expl->q4[i]) <= 1e-12) || !(fabs(r->q6[i] - expl->q6[i]) <= 1e-12)) {
        fail("atom %ld differs from the explicit images (%ld)", where, "q4q6", i, 0);
        break;
      }
    }
  }

  st = seams_chill_plus_n(sys, rc, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, r->labels,
                          N_ATOMS, r->counts9);
  if (st != SEAMS_OK) {
    fail("status %ld, expected %ld", where, "chill_plus", st, SEAMS_OK);
  } else {
    for (int i = 0; i < N_ATOMS; i++) {
      if (r->labels[i] != expl->labels[i]) {
        fail("atom %ld label differs from the explicit images, %ld", where, "chill_plus", i,
             r->labels[i]);
        break;
      }
    }
  }

  st = host_guest(sys, rc, r);
  if (st != SEAMS_OK) {
    fail("status %ld, expected %ld", where, "host_guest", st, SEAMS_OK);
  } else {
    for (int g = 0; g < N_HALF; g++) {
      const double want = sqrt(ref->nearest_d2[g]);
      if (r->guest[g] != 2 * g + 1 || r->nearest[g] != ref->nearest[g] ||
          r->nearest_type[g] != 1 || !(fabs(r->nearest_distance[g] - want) <= 1e-12 * want) ||
          r->contacts[g] != ref->contacts[g]) {
        fail("guest %ld nearest %ld differs from the exact image", where, "host_guest", g,
             r->nearest[g]);
        break;
      }
    }
  }

  st = contact_pairs(sys, r);
  if (st != SEAMS_OK) {
    fail("status %ld, expected %ld", where, "contact_pairs", st, SEAMS_OK);
  } else if ((long)r->mutual != ref->mutual || r->n_cation != N_HALF || r->n_anion != N_HALF) {
    fail("%ld mutual pairs, exact %ld", where, "contact_pairs", (long)r->mutual, ref->mutual);
  }
}

static void same_across_widths(const char *name, const results *a, const results *b) {
  char where[64];
  snprintf(where, sizeof where, "%s threads 1 and %d", name, WIDE);
  if (a->offsets != NULL && b->offsets != NULL &&
      (a->nnz != b->nnz || memcmp(a->offsets, b->offsets, (N_ATOMS + 1) * sizeof(int)) != 0 ||
       memcmp(a->cols, b->cols, a->nnz * sizeof(int)) != 0)) {
    fail("CSR differs (%ld, %ld)", where, "neighbours", (long)a->nnz, (long)b->nnz);
  }
  if (a->knn_built != b->knn_built ||
      (a->knn_built && (a->knn_nnz != b->knn_nnz ||
                        memcmp(a->knn_offsets, b->knn_offsets, (N_ATOMS + 1) * sizeof(int)) != 0 ||
                        memcmp(a->knn_cols, b->knn_cols, a->knn_nnz * sizeof(int)) != 0))) {
    fail("CSR differs (%ld, %ld)", where, "mutual neighbours", (long)a->knn_nnz,
         (long)b->knn_nnz);
  }
  if (memcmp(a->q4, b->q4, sizeof a->q4) != 0 || memcmp(a->q6, b->q6, sizeof a->q6) != 0) {
    fail("values differ (%ld, %ld)", where, "q4q6", 0, 0);
  }
  if (memcmp(a->labels, b->labels, sizeof a->labels) != 0 ||
      memcmp(a->counts9, b->counts9, sizeof a->counts9) != 0) {
    fail("labels differ (%ld, %ld)", where, "chill_plus", 0, 0);
  }
  if (memcmp(a->nearest, b->nearest, sizeof a->nearest) != 0 ||
      memcmp(a->nearest_distance, b->nearest_distance, sizeof a->nearest_distance) != 0 ||
      memcmp(a->contacts, b->contacts, sizeof a->contacts) != 0) {
    fail("columns differ (%ld, %ld)", where, "host_guest", 0, 0);
  }
  if (a->mutual != b->mutual) {
    fail("%ld and %ld mutual pairs", where, "contact_pairs", (long)a->mutual, (long)b->mutual);
  }
}

static void run_past(const cell_t *c, double rc, int width) {
  double xyz[3 * N_ATOMS];
  int types[N_ATOMS];
  int mols[N_ATOMS];
  make_atoms(c, xyz, types, mols);
  seams_system *sys = make_system(c, xyz, types, mols, N_ATOMS);
  char where[64];
  snprintf(where, sizeof where, "%s threads %d", c->name, width);
  results r;
  memset(&r, 0, sizeof r);

  expect_range(where, "neighbours",
               neighbours(sys, rc, SEAMS_GRAPH_CUTOFF, &r.offsets, &r.cols, &r.nnz),
               kMinimumImage);
  /* linkcell's cell construction rejects the cells near 1e160 as singular. */
  const seams_status st =
      neighbours(sys, rc, SEAMS_GRAPH_MUTUAL, &r.knn_offsets, &r.knn_cols, &r.knn_nnz);
  if (!is_range(st, kMinimumImage) && !is_range(st, kSingular)) {
    expect_range(where, "mutual neighbours", st, kMesh);
  }
  free_results(&r);
  expect_range(where, "q4q6", seams_q4q6_n(sys, rc, 0, r.q4, N_ATOMS, r.q6, N_ATOMS),
               kMinimumImage);
  expect_range(where, "chill_plus",
               seams_chill_plus_n(sys, rc, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, r.labels,
                                  N_ATOMS, r.counts9),
               kMinimumImage);
  expect_range(where, "host_guest", host_guest(sys, rc, &r), kMinimumImage);
  expect_range(where, "contact_pairs", contact_pairs(sys, &r), kMinimumImage);
  seams_system_destroy(sys);
}

int main(void) {
  static double xyz[3 * N_ATOMS];
  static int types[N_ATOMS];
  static int mols[N_ATOMS];
  static double images[3 * MAX_GHOSTS];
  static int image_types[MAX_GHOSTS];
  static results expl;
  static results run[2];
  /* An abort must not discard the lines before it. */
  setvbuf(stdout, NULL, _IONBF, 0);

  const double rc2 = off_grid_square(kCutoff);
  const double rc = sqrt(rc2);
  const cell_t *first = &kInRange[0];
  make_atoms(first, xyz, types, mols);
  exact_ref ref;
  ref.pair = malloc((size_t)N_ATOMS * N_ATOMS);
  ref.knn = malloc((size_t)N_ATOMS * N_ATOMS);
  if (ref.pair == NULL || ref.knn == NULL) {
    return 1;
  }
  build_reference(first, xyz, types, rc2, &ref);
  printf("%ld directed cutoff pairs, %ld of them missed by the fractional wrap, %ld mutual "
         "four-nearest pairs, %ld mutual nearest unlike pairs\n",
         ref.n_pairs, ref.wrap_misses, ref.n_knn, ref.mutual);
  if (ref.wrap_misses == 0) {
    printf("FAIL test data: the fractional wrap is exact on this cell\n");
    return 1;
  }

  /* The explicit images, with the same width at every step. */
  seams_set_num_threads(1);
  const size_t n_images = explicit_images(first, xyz, 2.0 * kCutoff + 0.0625, images,
                                          image_types);
  const cell_t wide_box = {"explicit images", 128.0, 128.0, 128.0, 0.0, 0.0, 0.0};
  {
    static double shifted[3 * MAX_GHOSTS];
    for (size_t k = 0; k < 3 * n_images; k++) {
      shifted[k] = images[k] + 32.0;
    }
    seams_system *sys = make_system(&wide_box, shifted, image_types, NULL, n_images);
    static double q4[MAX_GHOSTS];
    static double q6[MAX_GHOSTS];
    static int labels[MAX_GHOSTS];
    int counts9[9];
    if (seams_q4q6_n(sys, rc, 0, q4, n_images, q6, n_images) != SEAMS_OK ||
        seams_chill_plus_n(sys, rc, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels,
                           n_images, counts9) != SEAMS_OK) {
      printf("FAIL explicit images: %s\n", seams_last_error_message());
      return 1;
    }
    memcpy(expl.q4, q4, sizeof expl.q4);
    memcpy(expl.q6, q6, sizeof expl.q6);
    memcpy(expl.labels, labels, sizeof expl.labels);
    seams_system_destroy(sys);
  }
  printf("%s: %zu atoms with their images\n", wide_box.name, n_images);

  const int widths[2] = {1, WIDE};
  for (int c = 0; c < N_IN_RANGE; c++) {
    double cell_xyz[3 * N_ATOMS];
    make_atoms(&kInRange[c], cell_xyz, types, mols);
    if (memcmp(cell_xyz, xyz, sizeof cell_xyz) != 0) {
      printf("FAIL test data: %s has other atoms than %s\n", kInRange[c].name, first->name);
      return 1;
    }
    for (int w = 0; w < 2; w++) {
      seams_set_num_threads(widths[w]);
      seams_system *sys = make_system(&kInRange[c], xyz, types, mols, N_ATOMS);
      run_exact(&kInRange[c], sys, rc, &ref, &expl, &run[w], widths[w]);
      seams_system_destroy(sys);
    }
    same_across_widths(kInRange[c].name, &run[0], &run[1]);
    free_results(&run[0]);
    free_results(&run[1]);
  }

  for (int c = 0; c < N_PAST; c++) {
    const double cutoff = kPast[c].lz == 0.4 ? kGateCutoff : rc;
    for (int w = 0; w < 2; w++) {
      seams_set_num_threads(widths[w]);
      run_past(&kPast[c], cutoff, widths[w]);
    }
  }
  seams_set_num_threads(0);
  free(ref.pair);
  free(ref.knn);

  if (failures != 0) {
    printf("%d failures\n", failures);
    return 1;
  }
  printf("ok\n");
  return 0;
}
