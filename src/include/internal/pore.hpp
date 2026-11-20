//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#ifndef SEAMS_PORE_H_
#define SEAMS_PORE_H_

#include <mol_sys.hpp>

#include <stdexcept>
#include <vector>

/** @file pore.hpp
 *  @brief Triclinic fractional density, pore-frame profiles, and host-guest
 *         contacts.
 *
 *  A LAMMPS cell has edge vectors a, b, c recovered from the dump box.
 *  Fractional coordinates live in the unit cube. A slab of fractional
 *  thickness 1/n along a cell edge has volume V/n, with V = det(H).
 *  An arbitrary Cartesian direction bins the projection of that cube.
 *  Its bin volumes are the exact slice measures of the unit cube, times V,
 *  so a uniform gas is flat at N/V even when the projected slabs are not
 *  equally thick in the projection coordinate's density of states.
 *
 *  Cylindrical shells are Euclidean distances to a pore axis parallel to a
 *  cell edge. The shell volume is the area of the annulus inside the
 *  perpendicular face, times the length of that edge.
 */

namespace pore {

struct BoxError : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

struct Profile {
  std::vector<double> coordinate;
  std::vector<double> number_density;
  std::vector<double> mass_density;
  std::vector<double> volume;
};

struct Frame {
  double axis_unit[3];
  double centre[3];
  double centre_frac[3];
  int axis_index;
};

struct Profiles {
  Profile radial;
  Profile azimuthal;
  Profile axial;
};

struct HostGuest {
  std::vector<int> guest;
  std::vector<int> nearest;
  std::vector<int> nearest_type;
  std::vector<double> distance;
  std::vector<int> framework_types;
  std::vector<int> contacts;
};

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

//! Fractional density along cell axis 0, 1, or 2. Origin is in [0, 1) after wrapping.
Profile fractionalAxis(const Cloud &cloud, int type, int bins, int axis, double origin,
                       const double *masses);

//! Density of the reduced cell projected on a Cartesian direction.
Profile alongDirection(const Cloud &cloud, int type, int bins, const double direction[3],
                       double origin, const double *masses);

//! Pore axis and centre. axis < 0 and direction == nullptr selects the longest edge.
//! centre_frac == nullptr selects the framework centroid in the perpendicular plane.
Frame poreFrame(const Cloud &cloud, int framework_type, int axis, const double *direction,
                const double *centre_frac);

//! Radial, azimuthal, and axial profiles. The axis must be a cell edge.
//! r_max <= 0 covers the perpendicular face. Masses follow fractionalAxis.
Profiles poreProfiles(const Cloud &cloud, int type, int radial_bins, int azimuthal_bins,
                      int axial_bins, double r_max, int axis, const double *direction,
                      const double *centre_frac, int framework_type, double origin,
                      const double *masses);

//! Nearest framework site and contact counts inside cutoff, one row per guest.
//! framework_type 0 keeps every atom whose type is not guest_type.
//! distance is negative when the guest has no framework atom.
HostGuest hostGuest(const Cloud &cloud, int guest_type, int framework_type, double cutoff);

} // namespace pore

#endif // SEAMS_PORE_H_
