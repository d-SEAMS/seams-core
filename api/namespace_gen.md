<!-- https://docs.dseams.info/api/namespace_gen.md -->

# namespace gen {#namespace-gen}

## Overview {#overview}

Small generic functions that are shared by all namespaces. [More…](#details-namespacegen)

## Detailed Documentation {#detailed-documentation}

Small generic functions that are shared by all namespaces.

These are general functions (eg. for finding the periodic distance) which are required by several namespaces.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Nov 14, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Variables {#global-variables}

Uses Boost to get the value of pi.

#### Global Functions {#global-functions}

Inline function for converting radians->degrees.

**Parameters:**

| angle | The input angle, in radians |
| --- | --- |

**Returns:**

The input angle, in degrees

Eigen function for getting the angle (in radians) between the O O and O-H vectors.

Get the average, after excluding the outliers, using quartiles.

Inline generic function for calculating the median given a vector of the values.

**Parameters:**

| yCloud | The input PointCloud, which contains the particle coordinates, simulation box lengths etc. |
| --- | --- |
| input | The input vector with the values |

**Returns:**

The median value

Inline generic function for obtaining the unwrapped periodic distance between two particles, whose indices (not IDs) have been given.

**Parameters:**

| yCloud | The input PointCloud, which contains the particle coordinates, simulation box lengths etc. |
| --- | --- |
| iatom | The index of the i^{th} atom. |
| jatom | The index of the j^{th} atom. |

**Returns:**

The unwrapped periodic distance.

Inline generic function for obtaining the unwrapped periodic distance between two particles, whose indices (not IDs) have been given.

**Parameters:**

| yCloud | The input PointCloud, which contains the particle coordinates, simulation box lengths etc. |
| --- | --- |
| iatom | The index of the i^{th} atom. |
| jatom | The index of the j^{th} atom. |

**Returns:**

The unwrapped periodic distance.

Inline generic function for obtaining the unwrapped periodic distance between one particle and another point, whose index has been given.

**Parameters:**

| b | The box of yCloud, from makeFracBox. |
| --- | --- |
| yCloud | The input PointCloud, which contains the particle coordinates, simulation box lengths etc. |
| iatom | The index of the i^{th} atom. |
| xj |  |
| yj |  |
| zj | The other point. |

**Returns:**

The unwrapped periodic distance.

Inline generic function for obtaining the wrapped distance between two particles WITHOUT applying PBCs, whose indices (not IDs) have been given.

**Parameters:**

| yCloud | The input PointCloud, which contains the particle coordinates, simulation box lengths etc. |
| --- | --- |
| iatom | The index of the i^{th} atom. |
| jatom | The index of the j^{th} atom. |

**Returns:**

The wrapped distance.

Inline generic function for getting the relative unwrapped distance between two particles for each dimension.

The indices (not IDs) of the particles have been given.

**Parameters:**

| yCloud | The input PointCloud, which contains the particle coordinates, simulation box lengths etc. |
| --- | --- |
| iatom | The index of the i^{th} atom. |
| jatom | The index of the j^{th} atom. |

**Returns:**

The unwrapped relative distances for each dimension.

Inline generic function for sorting or comparing two particles, according to the atom ID when the entire Point objects have been passed.

**Parameters:**

| a | The input Point for A. |
| --- | --- |
| b | The input Point for B. |

**Returns:**

True if the atom ID of A is less than the atom ID of B

Generic function for printing all the struct information.

Shift particles (unwrapped coordinates).

Function for getting the angular distance between two quaternions. Returns the result in degrees

Function for tokenizing line strings into words (strings) delimited by whitespace.

This returns a vector with the words in it.

**Parameters:**

| line | The string containing the line to be tokenized |
| --- | --- |

Function for tokenizing line strings into a vector of doubles.

**Parameters:**

| line | The string containing the line to be tokenized |
| --- | --- |

Function for tokenizing line strings into a vector of ints.

**Parameters:**

| line | The string containing the line to be tokenized |
| --- | --- |

Function for checking if a file exists or not.

**Parameters:**

| name | The name of the file |
| --- | --- |

Calculates the complex vector, normalized by the number of nearest neighbours, of length $2l+1$.

**Parameters:**

| v | The complex vector to be normalized, of length 2l+1 |
| --- | --- |
| l | A free integer parameter |
| neigh | The number of nearest neighbours |

**Returns:**

length $2l+1$, normalized by the number of nearest neighbours

Given a pointCloud containing certain atom types, this returns a pointCloud containing atoms of only the desired type

Given a pointCloud set the inSlice bool for every atom, if the molecules are inside the specified (single) region. If even one atom of a molecule is inside the region, then all atoms of that molecule will be inside the region (irrespective of type)

Given a pointCloud set the inSlice bool for every atom, if the atoms are inside the specified (single) region. Does not handle atoms in molecules straddling the boundary

Given a particular molecule ID and a pointCloud set the inSlice bool for all atoms, with that molecule ID
