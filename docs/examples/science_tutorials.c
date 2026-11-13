/* Science tutorials beside docs/examples/science_tutorials.py.
 *
 * interface: CHILL+ and layer cubicity on the committed slab and nucleus.
 * cages: sII and sH signature counts.
 * brine: two fetched XYZ frames, oxygen labels and salt rejection.
 *
 * Layer period and layer kind follow scripts/science_render.py at
 * 702d85b2454e459dd14d9ec5fa3bd6b331e73862.
 */
#include "seams.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die(const char *what) {
  const char *err = seams_last_error();
  fprintf(stderr, "%s: %s\n", what, err != NULL ? err : "");
  exit(1);
}

static int read_dump(const char *path, double **xyz, int *n_out) {
  FILE *handle = fopen(path, "r");
  char line[512];
  int reading = 0;
  int cap = 0;
  int n = 0;
  double *pos = NULL;
  if (handle == NULL) {
    return 1;
  }
  while (fgets(line, sizeof line, handle) != NULL) {
    if (strncmp(line, "ITEM: ATOMS", 11) == 0) {
      reading = 1;
      continue;
    }
    if (!reading) {
      continue;
    }
    int id = 0;
    int type = 0;
    double x = 0;
    double y = 0;
    double z = 0;
    if (sscanf(line, "%d %d %lf %lf %lf", &id, &type, &x, &y, &z) != 5) {
      break;
    }
    if (n == cap) {
      cap = cap == 0 ? 256 : cap * 2;
      pos = realloc(pos, (size_t)cap * 3 * sizeof(double));
    }
    pos[3 * n] = x;
    pos[3 * n + 1] = y;
    pos[3 * n + 2] = z;
    n += 1;
  }
  fclose(handle);
  *xyz = pos;
  *n_out = n;
  return n == 0;
}

static void chill_file(const char *path, int **labels_out, int counts[9], int *n_out) {
  seams_system *sys = seams_system_create();
  if (seams_system_read(sys, path, 1, -1) != SEAMS_OK) {
    die(path);
  }
  int n = seams_system_count(sys);
  int *labels = calloc((size_t)n, sizeof(int));
  if (labels == NULL || seams_chill_plus(sys, 3.5, 0, 0, 0, labels, counts) != SEAMS_OK) {
    die("chill+");
  }
  seams_system_destroy(sys);
  *labels_out = labels;
  *n_out = n;
}

static double layer_period(const double *projection, int n) {
  double mean = 0;
  for (int i = 0; i < n; i++) {
    mean += projection[i];
  }
  mean /= (double)n;
  double best_score = -1;
  double best = 3.6;
  for (int step = 0; step < 16; step++) {
    double period = 3.30 + (4.05 - 3.30) * (double)step / 15.0;
    double re = 0;
    double im = 0;
    for (int i = 0; i < n; i++) {
      double angle = (2.0 * 3.14159265358979323846 / period) * (projection[i] - mean);
      re += cos(angle);
      im += sin(angle);
    }
    double score = hypot(re, im) / (double)n;
    if (score > best_score) {
      best_score = score;
      best = period;
    }
  }
  return best;
}

static double layer_score(const double *projection, int n, double period) {
  double mean = 0;
  for (int i = 0; i < n; i++) {
    mean += projection[i];
  }
  mean /= (double)n;
  double re = 0;
  double im = 0;
  for (int i = 0; i < n; i++) {
    double angle = (2.0 * 3.14159265358979323846 / period) * (projection[i] - mean);
    re += cos(angle);
    im += sin(angle);
  }
  return hypot(re, im) / (double)n;
}

/* Returns the stacking string. axis is 0, 1, or 2. */
static void layer_string(const double *xyz, const int *chill, int n, int axis, char *out) {
  int crystal = 0;
  for (int i = 0; i < n; i++) {
    if (chill[i] == 0 || chill[i] == 1) {
      crystal += 1;
    }
  }
  double *projection = malloc((size_t)crystal * sizeof(double));
  int *kinds = malloc((size_t)crystal * sizeof(int));
  int cursor = 0;
  double lo = 0;
  for (int i = 0; i < n; i++) {
    if (chill[i] != 0 && chill[i] != 1) {
      continue;
    }
    projection[cursor] = xyz[3 * i + axis];
    kinds[cursor] = chill[i];
    if (cursor == 0 || projection[cursor] < lo) {
      lo = projection[cursor];
    }
    cursor += 1;
  }
  double period = layer_period(projection, crystal);
  int bins_max = 0;
  int *bin = malloc((size_t)crystal * sizeof(int));
  for (int i = 0; i < crystal; i++) {
    bin[i] = (int)floor((projection[i] - lo) / period);
    if (bin[i] > bins_max) {
      bins_max = bin[i];
    }
  }
  int length = 0;
  for (int index = 0; index <= bins_max; index++) {
    int hexagonal = 0;
    int cubic = 0;
    for (int i = 0; i < crystal; i++) {
      if (bin[i] != index) {
        continue;
      }
      if (kinds[i] == 1) {
        hexagonal += 1;
      } else {
        cubic += 1;
      }
    }
    if (hexagonal + cubic == 0) {
      continue;
    }
    out[length++] = hexagonal >= cubic ? 'H' : 'C';
  }
  out[length] = '\0';
  free(projection);
  free(kinds);
  free(bin);
}

static int strongest_axis(const double *xyz, const int *chill, int n) {
  int crystal = 0;
  for (int i = 0; i < n; i++) {
    if (chill[i] == 0 || chill[i] == 1) {
      crystal += 1;
    }
  }
  double *projection = malloc((size_t)crystal * sizeof(double));
  int best = 0;
  double best_score = -1;
  for (int axis = 0; axis < 3; axis++) {
    int cursor = 0;
    for (int i = 0; i < n; i++) {
      if (chill[i] == 0 || chill[i] == 1) {
        projection[cursor++] = xyz[3 * i + axis];
      }
    }
    double period = layer_period(projection, crystal);
    double score = layer_score(projection, crystal, period);
    if (score > best_score) {
      best_score = score;
      best = axis;
    }
  }
  free(projection);
  return best;
}

static int cmd_interface(void) {
  const char *interface = "input/traj/interface_f0.lammpstrj";
  const char *nucleus = "input/traj/nucleation_f1.lammpstrj";
  double *xyz = NULL;
  int n = 0;
  int *labels = NULL;
  int counts[9] = {0};
  char stack[64];
  if (read_dump(interface, &xyz, &n) != 0) {
    die("interface dump");
  }
  chill_file(interface, &labels, counts, &n);
  layer_string(xyz, labels, n, 2, stack);
  printf("interface hexagonal %d water %d interfacial %d cubic %d\n", counts[1], counts[2],
         counts[3], counts[0]);
  printf("interface layers %s\n", stack);
  free(xyz);
  free(labels);
  if (read_dump(nucleus, &xyz, &n) != 0) {
    die("nucleus dump");
  }
  chill_file(nucleus, &labels, counts, &n);
  int axis = strongest_axis(xyz, labels, n);
  layer_string(xyz, labels, n, axis, stack);
  const char *axes = "xyz";
  printf("nucleus cubic %d hexagonal %d water %d interfacial %d\n", counts[0], counts[1],
         counts[2], counts[3]);
  printf("nucleus axis %c layers %s\n", axes[axis], stack);
  free(xyz);
  free(labels);
  return 0;
}

static int cages_of(const char *path, const char *spec) {
  seams_system *sys = seams_system_create();
  if (seams_system_read(sys, path, 1, -1) != SEAMS_OK) {
    die(path);
  }
  size_t cages = 0;
  size_t atoms = 0;
  if (seams_cages_signature(sys, 3.5, 0, 0, spec, &cages, &atoms) != SEAMS_OK) {
    die(spec);
  }
  seams_system_destroy(sys);
  return (int)cages;
}

static int cmd_cages(void) {
  printf("cages sII 512 %d\n", cages_of("input/traj/genice_sII.lammpstrj", "512"));
  printf("cages sII 51264 %d\n", cages_of("input/traj/genice_sII.lammpstrj", "51264"));
  printf("cages sH 512 %d\n", cages_of("input/traj/genice_sH.lammpstrj", "512"));
  printf("cages sH sh %d\n", cages_of("input/traj/genice_sH.lammpstrj", "sh"));
  printf("cages sH 51268 %d\n", cages_of("input/traj/genice_sH.lammpstrj", "51268"));
  return 0;
}

static int read_xyz(const char *path, double **xyz, int **types, char chill[][16], int *n_out,
                    double box9[9]) {
  FILE *handle = fopen(path, "r");
  char line[1024];
  if (handle == NULL || fgets(line, sizeof line, handle) == NULL) {
    return 1;
  }
  int n = atoi(line);
  if (fgets(line, sizeof line, handle) == NULL) {
    fclose(handle);
    return 1;
  }
  double lattice[9] = {0};
  double origin[3] = {0};
  char *lat = strstr(line, "Lattice=\"");
  char *org = strstr(line, "origin=\"");
  if (lat == NULL || org == NULL) {
    fclose(handle);
    return 1;
  }
  sscanf(lat + 9, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &lattice[0], &lattice[1], &lattice[2],
         &lattice[3], &lattice[4], &lattice[5], &lattice[6], &lattice[7], &lattice[8]);
  sscanf(org + 8, "%lf %lf %lf", &origin[0], &origin[1], &origin[2]);
  box9[0] = origin[0];
  box9[1] = origin[1];
  box9[2] = origin[2];
  box9[3] = lattice[0];
  box9[4] = lattice[4];
  box9[5] = lattice[8];
  box9[6] = box9[7] = box9[8] = 0;
  double *pos = calloc((size_t)n * 3, sizeof(double));
  int *kind = calloc((size_t)n, sizeof(int));
  for (int i = 0; i < n; i++) {
    char species[8];
    if (fgets(line, sizeof line, handle) == NULL) {
      fclose(handle);
      return 1;
    }
    if (sscanf(line, "%7s %lf %lf %lf %15s", species, &pos[3 * i], &pos[3 * i + 1],
               &pos[3 * i + 2], chill[i]) != 5) {
      fclose(handle);
      return 1;
    }
    if (strcmp(species, "O") == 0) {
      kind[i] = 1;
    } else if (strcmp(species, "H") == 0) {
      kind[i] = 2;
    } else if (strcmp(species, "Na") == 0) {
      kind[i] = 3;
    } else if (strcmp(species, "Cl") == 0) {
      kind[i] = 4;
    }
  }
  fclose(handle);
  *xyz = pos;
  *types = kind;
  *n_out = n;
  return 0;
}

static const char *LABELS[] = {"cubic",     "hexagonal", "water",
                               "interfacial", "clathrate", "interClathrate",
                               "unclassified", "reCubic",  "reHex"};

static void one_brine(const char *path, const char *tag) {
  double *xyz = NULL;
  int *types = NULL;
  int n = 0;
  double box9[9];
  char(*chill)[16] = NULL;
  /* 4608 is the largest frame in this tutorial. Grow if the header says more. */
  FILE *probe = fopen(path, "r");
  char first[64];
  if (probe == NULL || fgets(first, sizeof first, probe) == NULL) {
    die(path);
  }
  fclose(probe);
  int declared = atoi(first);
  chill = calloc((size_t)declared, 16);
  if (read_xyz(path, &xyz, &types, chill, &n, box9) != 0) {
    die(path);
  }
  seams_system *sys = seams_system_create();
  if (seams_system_set_frame(sys, xyz, n, box9, types, NULL) != SEAMS_OK) {
    die("set_frame");
  }
  int *labels = calloc((size_t)n, sizeof(int));
  int counts[9] = {0};
  if (seams_chill_plus(sys, 3.5, 1, 0, 0, labels, counts) != SEAMS_OK) {
    die("chill+");
  }
  int agree = 0;
  int oxygen = 0;
  int histogram[9] = {0};
  for (int i = 0; i < n; i++) {
    if (types[i] != 1) {
      continue;
    }
    oxygen += 1;
    histogram[labels[i]] += 1;
    if (strcmp(chill[i], LABELS[labels[i]]) == 0) {
      agree += 1;
    }
  }
  int *ice = calloc((size_t)n, sizeof(int));
  for (int i = 0; i < n; i++) {
    ice[i] = types[i] == 1 && (labels[i] == 0 || labels[i] == 1);
  }
  seams_brine *front = NULL;
  if (seams_brine_analyse(sys, ice, (size_t)n, 1, 3, 20, 2, &front) != SEAMS_OK) {
    die("brine");
  }
  double thickness = 0;
  int found = 0;
  double c_ice = 0;
  double coeff = 0;
  int defined = 0;
  if (seams_brine_fronts(front, NULL, NULL, &thickness, &found) != SEAMS_OK) {
    die("fronts");
  }
  if (seams_brine_rejection(front, &c_ice, NULL, &coeff, &defined) != SEAMS_OK) {
    die("rejection");
  }
  printf("brine %s hexagonal %d water %d interfacial %d cubic %d agree %d\n", tag,
         histogram[1], histogram[2], histogram[3], histogram[0], agree);
  printf("brine %s thickness %.4f k %.1f found %d defined %d\n", tag, thickness, coeff, found,
         defined);
  seams_brine_destroy(front);
  seams_system_destroy(sys);
  free(xyz);
  free(types);
  free(labels);
  free(ice);
  free(chill);
  (void)oxygen;
}

static int cmd_brine(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr, "usage: science_tutorials brine T0 TEND\n");
    return 1;
  }
  one_brine(argv[2], "t0");
  one_brine(argv[3], "tend");
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: science_tutorials interface|cages|brine\n");
    return 1;
  }
  if (seams_set_num_threads(1) != SEAMS_OK) {
    die("threads");
  }
  if (strcmp(argv[1], "interface") == 0) {
    return cmd_interface();
  }
  if (strcmp(argv[1], "cages") == 0) {
    return cmd_cages();
  }
  if (strcmp(argv[1], "brine") == 0) {
    return cmd_brine(argc, argv);
  }
  fprintf(stderr, "unknown mode %s\n", argv[1]);
  return 1;
}
