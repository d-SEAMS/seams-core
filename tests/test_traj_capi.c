/* Copyright (c) d-SEAMS developers
 * Analytic-oracle C11 tests for the trajectory and time-correlation C API.
 */
#include "seams.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int check_close(double got, double want, double tol, const char *label) {
  if (fabs(got - want) > tol) {
    fprintf(stderr, "FAIL %s: got %.12g, want %.12g (tol %.3g)\n", label, got, want, tol);
    return 1;
  }
  return 0;
}

int main(void) {
  int fails = 0;

  /* 1. Multi-crossing triclinic unwrapping in C */
  const double box9[9] = {0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 2.5, -1.0, 1.5};
  double wrapped[15];
  double truth[15];
  for (int f = 0; f < 5; ++f) {
    double rx = 8.0 + 3.5 * f;
    double ry = 7.0 - 3.2 * f;
    double rz = 9.0 + 4.1 * f;
    truth[3 * f + 0] = rx;
    truth[3 * f + 1] = ry;
    truth[3 * f + 2] = rz;
    double sz = rz / box9[5];
    double sy = (ry - box9[8] * sz) / box9[4];
    double sx = (rx - box9[6] * sy - box9[7] * sz) / box9[3];
    sx -= floor(sx);
    sy -= floor(sy);
    sz -= floor(sz);
    wrapped[3 * f + 0] = box9[3] * sx + box9[6] * sy + box9[7] * sz;
    wrapped[3 * f + 1] = box9[4] * sy + box9[8] * sz;
    wrapped[3 * f + 2] = box9[5] * sz;
  }
  double unwrapped[15] = {0};
  if (seams_traj_unwrap_positions(wrapped, 5, 1, box9, 1, unwrapped, 15) != SEAMS_OK) {
    fprintf(stderr, "FAIL unwrap status: %s\n", seams_last_error_message());
    return 1;
  }
  const double shift[3] = {unwrapped[0] - truth[0], unwrapped[1] - truth[1],
                           unwrapped[2] - truth[2]};
  for (int i = 0; i < 15; ++i) {
    fails += check_close(unwrapped[i], truth[i] + shift[i % 3], 1e-11, "triclinic_unwrap");
  }

  /* 1b. Skewed LAMMPS-legal cell where rounding each fractional coordinate
   * picks an image one a-vector longer than the Euclidean one; the third frame
   * scales the cell by 3%. */
  {
    const double c0[9] = {0.0, 0.0, 0.0, 18.03850526, 10.63623724, 11.40838967,
                          -8.77049545, 6.89505882, 2.90161048};
    const double scale = 1.03;
    double boxes[27];
    for (int f = 0; f < 3; ++f) {
      for (int i = 0; i < 9; ++i) {
        boxes[9 * f + i] = (f == 2 && i >= 3) ? scale * c0[i] : c0[i];
      }
    }
    const double e[3] = {3.42737868, 3.53897334, -3.29810805};
    double u[9] = {4.1, 2.7, 5.3};
    for (int k = 0; k < 3; ++k) {
      u[3 + k] = u[k] + e[k];
      u[6 + k] = u[3 + k] + scale * e[k];
    }
    double w[9];
    for (int f = 0; f < 3; ++f) {
      const double *b = boxes + 9 * f;
      double sz = u[3 * f + 2] / b[5];
      double sy = (u[3 * f + 1] - b[8] * sz) / b[4];
      double sx = (u[3 * f] - b[6] * sy - b[7] * sz) / b[3];
      sx -= floor(sx);
      sy -= floor(sy);
      sz -= floor(sz);
      w[3 * f + 0] = b[3] * sx + b[6] * sy + b[7] * sz;
      w[3 * f + 1] = b[4] * sy + b[8] * sz;
      w[3 * f + 2] = b[5] * sz;
    }
    double uw[9] = {0};
    if (seams_traj_unwrap_positions(w, 3, 1, boxes, 3, uw, 9) != SEAMS_OK) {
      fprintf(stderr, "FAIL skewed unwrap status: %s\n", seams_last_error_message());
      return 1;
    }
    for (int i = 0; i < 9; ++i) {
      fails += check_close(uw[i], u[i], 1e-9, "skewed_cell_euclidean_unwrap");
    }
  }

  /* 2. Streaming MSD with axial vs lateral split */
  const double z_axis[3] = {0.0, 0.0, 1.0};
  seams_msd *msd = NULL;
  if (seams_msd_create(1, 0, 1, 0.5, z_axis, &msd) != SEAMS_OK) {
    return 1;
  }
  for (int f = 0; f < 10; ++f) {
    const double t = 0.5 * f;
    const double p[3] = {2.0 * t, 0.0, 3.0 * t};
    if (seams_msd_push(msd, p, 1, NULL, 0) != SEAMS_OK) {
      return 1;
    }
  }
  double lt[10], mtot[10], maxl[10], mlat[10];
  size_t norig[10], nlags = 0;
  if (seams_msd_curve(msd, lt, mtot, maxl, mlat, norig, 10, &nlags) != SEAMS_OK || nlags != 10) {
    return 1;
  }
  fails += check_close(maxl[4], 9.0 * lt[4] * lt[4], 1e-11, "msd_axial_lag4");
  fails += check_close(mlat[4], 4.0 * lt[4] * lt[4], 1e-11, "msd_lateral_lag4");
  fails += check_close(mtot[4], 13.0 * lt[4] * lt[4], 1e-11, "msd_total_lag4");
  double alpha = 0.0;
  if (seams_msd_fit(msd, 0.5, 4.0, NULL, NULL, NULL, &alpha) != SEAMS_OK) {
    return 1;
  }
  fails += check_close(alpha, 2.0, 1e-10, "ballistic_alpha");
  seams_msd_destroy(msd);

  /* 2b. Skipped log-log fit: alpha is NaN, the Einstein slope stays exact. One
   * origin; lag 1 moves 1e-11 (MSD 1e-22 < 1e-20), lag 2 moves 2. */
  {
    const double x_axis[3] = {1.0, 0.0, 0.0};
    seams_msd *skip = NULL;
    if (seams_msd_create(1, 0, 10, 1.0, x_axis, &skip) != SEAMS_OK) {
      return 1;
    }
    const double xs[3] = {0.0, 1e-11, 2.0};
    for (int f = 0; f < 3; ++f) {
      const double p[3] = {xs[f], 0.0, 0.0};
      if (seams_msd_push(skip, p, 1, NULL, 0) != SEAMS_OK) {
        return 1;
      }
    }
    double dtot = -1.0, dax = -1.0, dlat = -1.0, a = 0.0;
    if (seams_msd_fit(skip, 0.0, 2.0, &dtot, &dax, &dlat, &a) != SEAMS_OK) {
      fprintf(stderr, "FAIL skipped-alpha status: %s\n", seams_last_error_message());
      return 1;
    }
    if (!isnan(a)) {
      fprintf(stderr, "FAIL skipped_alpha_nan: got %.12g\n", a);
      ++fails;
    }
    fails += check_close(dtot, 1.0 / 3.0, 1e-12, "skipped_alpha_d_total");
    fails += check_close(dax, 1.0, 1e-12, "skipped_alpha_d_axial");
    fails += check_close(dlat, 0.0, 1e-12, "skipped_alpha_d_lateral");
    seams_msd_destroy(skip);

    /* Two frames leave one positive lag: the line fits, alpha does not. */
    seams_msd *two = NULL;
    if (seams_msd_create(1, 0, 1, 1.0, x_axis, &two) != SEAMS_OK) {
      return 1;
    }
    const double q0[3] = {0.0, 0.0, 0.0}, q1[3] = {1.0, 0.0, 0.0};
    if (seams_msd_push(two, q0, 1, NULL, 0) != SEAMS_OK ||
        seams_msd_push(two, q1, 1, NULL, 0) != SEAMS_OK) {
      return 1;
    }
    a = 0.0;
    if (seams_msd_fit(two, 0.0, 1.0, &dtot, NULL, NULL, &a) != SEAMS_OK) {
      return 1;
    }
    if (!isnan(a)) {
      fprintf(stderr, "FAIL two_point_alpha_nan: got %.12g\n", a);
      ++fails;
    }
    fails += check_close(dtot, 1.0 / 6.0, 1e-12, "two_point_d_total");
    seams_msd_destroy(two);
  }

  /* Counts whose 3 * n or max_lag + 1 wraps size_t are refused at create.
   * Accepted, 3 * (SIZE_MAX / 3 + 1) sizes the coordinate store at 2 and
   * max_lag = SIZE_MAX sizes the lag bins at 0, and push writes past both. */
  {
    const double z_axis[3] = {0.0, 0.0, 1.0};
    const size_t wraps3 = SIZE_MAX / 3 + 1;
    struct {
      const char *label;
      size_t n;
      size_t max_lag;
    } msd_cases[] = {
        {"msd natoms 3n wraps", wraps3, 4},
        {"msd natoms SIZE_MAX", SIZE_MAX, 4},
        {"msd max_lag SIZE_MAX", 1, SIZE_MAX},
        {"msd max_lag SIZE_MAX - 1", 1, SIZE_MAX - 1},
    };
    for (size_t c = 0; c < sizeof msd_cases / sizeof msd_cases[0]; ++c) {
      seams_msd *m = (seams_msd *)&fails;
      const seams_status st =
          seams_msd_create(msd_cases[c].n, msd_cases[c].max_lag, 1, 1.0, z_axis, &m);
      if (st != SEAMS_ERR_RANGE || m != NULL) {
        fprintf(stderr, "FAIL %s: status %d, handle %s\n", msd_cases[c].label, (int)st,
                m == NULL ? "NULL" : "set");
        ++fails;
        if (st == SEAMS_OK) {
          seams_msd_destroy(m);
        }
      } else {
        printf("ok %s: %s\n", msd_cases[c].label, seams_last_error_message());
      }
    }
    struct {
      const char *label;
      size_t n;
      size_t max_lag;
    } reo_cases[] = {
        {"reorient n_vectors 3n wraps", wraps3, 4},
        {"reorient n_vectors SIZE_MAX", SIZE_MAX, 4},
        {"reorient max_lag SIZE_MAX", 1, SIZE_MAX},
    };
    for (size_t c = 0; c < sizeof reo_cases / sizeof reo_cases[0]; ++c) {
      seams_reorient *r = (seams_reorient *)&fails;
      const seams_status st =
          seams_reorient_create(reo_cases[c].n, reo_cases[c].max_lag, 1, 1.0, &r);
      if (st != SEAMS_ERR_RANGE || r != NULL) {
        fprintf(stderr, "FAIL %s: status %d, handle %s\n", reo_cases[c].label, (int)st,
                r == NULL ? "NULL" : "set");
        ++fails;
        if (st == SEAMS_OK) {
          seams_reorient_destroy(r);
        }
      } else {
        printf("ok %s: %s\n", reo_cases[c].label, seams_last_error_message());
      }
    }
    seams_hbond_tcf *h = (seams_hbond_tcf *)&fails;
    const seams_status st = seams_hbond_tcf_create(SIZE_MAX, 1, 1.0, 1, &h);
    if (st != SEAMS_ERR_RANGE || h != NULL) {
      fprintf(stderr, "FAIL hbond max_lag SIZE_MAX: status %d, handle %s\n", (int)st,
              h == NULL ? "NULL" : "set");
      ++fails;
      if (st == SEAMS_OK) {
        seams_hbond_tcf_destroy(h);
      }
    } else {
      printf("ok hbond max_lag SIZE_MAX: %s\n", seams_last_error_message());
    }

    const size_t too_many_boxes = (size_t)PTRDIFF_MAX / (9 * sizeof(double)) + 1;
    const size_t wraps9 = SIZE_MAX / 9 + 1;
    const double dummy_xyz[3] = {1.0, 2.0, 3.0};
    double dummy_out[3] = {0.0, 0.0, 0.0};
    const int dummy_mols[1] = {1};
    const size_t bad_boxes[] = {too_many_boxes, wraps9};
    for (size_t c = 0; c < sizeof bad_boxes / sizeof bad_boxes[0]; ++c) {
      const size_t nb = bad_boxes[c];
      const seams_status st_pos = seams_traj_unwrap_positions(
          dummy_xyz, nb, 1, box9, nb, dummy_out, nb * 3);
      if (st_pos != SEAMS_ERR_RANGE) {
        fprintf(stderr, "FAIL unwrap_positions nboxes=%zu: status %d\n", nb, (int)st_pos);
        ++fails;
      } else {
        printf("ok unwrap_positions nboxes=%zu: %s\n", nb, seams_last_error_message());
      }
      size_t nm = 99;
      const seams_status st_com = seams_traj_unwrap_centres(
          dummy_xyz, nb, 1, box9, nb, dummy_mols, 1, NULL, 0, NULL, 0, NULL, 0, &nm);
      if (st_com != SEAMS_ERR_RANGE) {
        fprintf(stderr, "FAIL unwrap_centres nboxes=%zu: status %d\n", nb, (int)st_com);
        ++fails;
      } else {
        printf("ok unwrap_centres nboxes=%zu: %s\n", nb, seams_last_error_message());
      }
    }
  }

  return fails == 0 ? 0 : 1;
}
