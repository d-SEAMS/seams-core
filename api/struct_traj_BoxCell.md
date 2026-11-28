<!-- https://docs.dseams.info/api/struct_traj_BoxCell.md -->

# struct traj::BoxCell {#struct-traj-boxcell}

## Overview {#overview}

Restricted-triclinic periodic simulation cell. [More…](#details-structtraj-1-1-box-cell)

## Detailed Documentation {#detailed-documentation}

Restricted-triclinic periodic simulation cell.

Edge matrix H has columns a = (lx, 0, 0), b = (xy, ly, 0), c = (xz, yz, lz) and origin (ox, oy, oz).

### Methods {#methods}

Euclidean minimum-image Cartesian displacement for raw delta (dx, dy, dz).

The fractional image inside the Smith ball, otherwise micDeltaSearch. Independent rounding of fractional coordinates is not the minimum image on a tilted cell.

Euclidean minimum image by exhaustive search over the lattice shifts that can beat the fractional image.
