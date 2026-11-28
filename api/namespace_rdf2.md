<!-- https://docs.dseams.info/api/namespace_rdf2.md -->

# namespace rdf2 {#namespace-rdf2}

## Overview {#overview}

Defines RDF-specific functions. [More…](#details-namespacerdf2)

## Detailed Documentation {#detailed-documentation}

Defines RDF-specific functions.

The radial distribution function $g(r)$, or pair distribution function. It can be used to illuminate features of short-range and long-range order. The RDF is the probability of finding a particle at a distance of $r$ from a tagged reference particle, relative to that of an ideal gas. For a system of $N$ particles, the pair correlation function for $N(N-1)$ pairs is:

$$
\rho_N^{(2)}(r,r') = ⟨\sum_{i=1}^{N} \sum_{j=1,j \neq i}^{N} \delta
(r-r_i) \delta (r'-r_j) \langle
$$

The code essentially bins distances between pairs of particles, and normalizes the resulting histogram is normalized with respect to an ideal gas. The algorithm for the calculation of $g(r)$ may be divided into the following steps:

1. **Initialization:** The $g(r)$ array is initialized to zero.
2. **Sampling:** The histogram is added to for a particular bin, if the distance of a pair of atoms falls within the $r$ associated with the bin.
3. **Normalization:** Every bin of the $g(r)$ array is normalized by the product of the number of ideal gas particles in that bin, and the number of particles and number of frames.

To account for the geometry of the slab-like quasi-two-dimensional system, the RDF must be additionally normalized by a form factor $f$ such that $f = 1 - \frac{r}{2h}$ for $r \le h$ and $f = \frac{h}{2r}$ for $r > h$; where $h$ is the height of the quasi-two-dimensional system. That is what `rdf2::normalizeRDF` implements.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Dec 10, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Main function for calculating the RDF for the same type of particle: calls other functions for initializing, sampling and normalizing the RDF

Samples the RDF histogram at every step.

Normalize the histogram.

Gets the lengths of the volume slice of the quasi-two-dimensional system. box.size()>=6 returns dump H {lx, ly, lz}; a length-3 box stays the particle AABB so slab tests keep the occupied extents.

Gets the plane area from the volume lengths vector.
