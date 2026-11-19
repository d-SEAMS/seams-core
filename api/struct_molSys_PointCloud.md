<!-- https://docs.dseams.info/api/struct_molSys_PointCloud.md -->

# template struct molSys::PointCloud {#template-struct-molsys-pointcloud}

## Overview {#overview}

This contains a collection of points; contains information for a particular frame. [More…](#details-structmol-sys-1-1-point-cloud)

## Detailed Documentation {#detailed-documentation}

This contains a collection of points; contains information for a particular frame.

Specifically

- A vector of [Point](struct_molSys_Point.md#doxid-structmol-sys-1-1-point) structs
- The current frame number
- The number of particles in the current frame
- A vector for the simulation box lengths in each dimension
- A vector containing the absolute lower box coordinates

### Fields {#fields}

Collection of points.

Current frame number.

Number of atoms.

Periodic box lengths.

xlo, ylo, zlo
