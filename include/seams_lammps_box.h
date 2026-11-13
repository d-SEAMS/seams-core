#ifndef SEAMS_LAMMPS_BOX_H
#define SEAMS_LAMMPS_BOX_H

/* LAMMPS BOX BOUNDS to ox, oy, oz, lx, ly, lz, xy, xz, yz.
 * lo/hi are the bound column. tilt[0..2] are xy, xz, yz (0 when the
 * dump is orthorhombic). lx is xhi - xlo, not the bound span. */

static inline void seams_lammps_bounds(const double lo[3], const double hi[3],
                                      const double tilt[3], double box9[9]) {
  const double xy = tilt[0];
  const double xz = tilt[1];
  const double yz = tilt[2];
  const double xmin = xy < xz ? (xy < 0.0 ? xy : 0.0) : (xz < 0.0 ? xz : 0.0);
  const double xy_xz = xy + xz;
  const double xmin_all = xy_xz < xmin ? xy_xz : xmin;
  const double xmax = xy > xz ? (xy > 0.0 ? xy : 0.0) : (xz > 0.0 ? xz : 0.0);
  const double xmax_all = xy_xz > xmax ? xy_xz : xmax;
  const double ymin = yz < 0.0 ? yz : 0.0;
  const double ymax = yz > 0.0 ? yz : 0.0;
  box9[0] = lo[0] - xmin_all;
  box9[1] = lo[1] - ymin;
  box9[2] = lo[2];
  box9[3] = (hi[0] - lo[0]) - xmax_all + xmin_all;
  box9[4] = (hi[1] - lo[1]) - ymax + ymin;
  box9[5] = hi[2] - lo[2];
  box9[6] = xy;
  box9[7] = xz;
  box9[8] = yz;
}

#endif
