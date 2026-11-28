<!-- https://docs.dseams.info/api/namespace_pore.md -->

# namespace pore {#namespace-pore}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Functions {#global-functions}

Fractional density along cell axis 0, 1, or 2. Origin is in \[0, 1) after wrapping.

Density of the reduced cell projected on a Cartesian direction.

Pore axis and centre. axis -1 and direction == nullptr selects the longest edge. Any other negative axis is rejected. centre_frac == nullptr selects the framework centroid in the perpendicular plane.

Radial, azimuthal, and axial profiles. The axis must be a cell edge. r_max of zero covers the perpendicular face. A negative value is rejected. Azimuth zero points along the perpendicular projection of cell edge (axis + 1) mod 3. Masses follow fractionalAxis.

Area of a convex polygon clipped to the disk of `radius` centred at the origin. Vertices are an (x, y) ring. A non-positive radius returns 0. Once `radius` reaches `cover`, the result is `full`, which is the polygon’s own area when the disk contains it.

Nearest framework site and contact counts inside cutoff, one row per guest. framework_type 0 keeps every atom whose type is not guest_type. framework_type equal to guest_type is rejected. distance is negative when the guest has no framework atom.
