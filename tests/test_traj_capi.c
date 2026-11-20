/* Copyright (c) d-SEAMS developers
 * Analytic-oracle C11 tests for the trajectory and time-correlation C API.
 */
#include "seams.h"

#include <math.h>
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

  return fails == 0 ? 0 : 1;
}
