<!-- https://docs.dseams.info/api/namespace_bond.md -->

# namespace bond {#namespace-bond}

## Overview {#overview}

Functions for bond-related analyses. [More…](#details-namespacebond)

## Detailed Documentation {#detailed-documentation}

Functions for bond-related analyses.

This namespace contains functions that are used for determining the hydrogen bonds, or lists of bonded atoms for write-outs to data files and visualization.

These functions are distinct from the neighbour-list building functions. Bonds can be built from cages, ring lists etc.

Hydrogen bonds are determined using a strict geometric criterion.

A hydrogen bond between two water molecules exists when:

1. The distance between the acceptor oxygen atom and the donor hydrogen atom is less than 2.42 Angstrom
2. The angle between the O-O vector and the O-H vector should be less than 30 degrees

Those numbers are the water defaults on populateHbonds and populateHbondsWithInputClouds. populateHbondsFromDonors takes the same defaults so a water donor-H set stays comparable; other chemistries pass their own cutoffs. They are not an ionic-liquid default.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Nov 13, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Create a vector of vectors (similar to the neighbour list conventions) containing the hydrogen bond connectivity information. Decides the existence of the hydrogen bond depending on the O O and O H vectors from the neighbour list already constructed. The geometric criterion defaults to the water values (acceptor H distance below 2.42 Angstrom, O-O-H angle below 30 degrees); other hydrogen-bonded systems pass their own thresholds

Create a vector of vectors (similar to the neighbour list conventions) containing the hydrogen bond connectivity information. Decides the existence of the hydrogen bond depending on the O O and O H vectors from the neighbour list already constructed, taking a PointCloud for the H atoms as input

Geometric test for one donor-acceptor assignment. donorHs are hCloud indices on the donor. O-O is gen::relDist(yCloud, acceptor, donor). O-H uses [gen::relDistFromPoint](namespace_gen.md#doxid-namespacegen-1a2fbd5127b1fdc22c99f5cd6d57d4bf86) (dump MIC).

Hydrogen-bond network from an explicit donor-H set. donorHs are hCloud indices; each is paired with the heavy atom that shares its molID. Water callers build donorHs from hAtomMolList (two H per molecule).

Directed donor-hydrogen-acceptor edges on one cloud.

Atoms whose type equals hydrogenType are hydrogens. Every other atom is a heavy atom. Each hydrogen is given to the nearest heavy atom of a donor type that shares its molID; the lower index wins an equal distance. A hydrogen with no such heavy atom donates nothing. An edge is kept when the acceptor-hydrogen minimum image is strictly shorter than distCutoff and the acceptor-first angle is less than or equal to angleCutoff. The neighbour list keeps a donor-acceptor pair only when its minimum image is strictly shorter than pairCutoff. excludeSameMol drops a pair that shares a molID. An empty donor or acceptor type list returns no edges. Results are sorted by donor, hydrogen, acceptor.

Calculates the distance of the hydrogen bond between O and H (of different atoms), given the respective pointClouds and the indices to each atom

Create a vector of vectors containing bond connectivity information. May contain duplicates! Gets the bond information from the vector of vectors containing the neighbour list by index

Create a vector of vectors containing bond connectivity information Gets the bond information from the vector of vectors containing the neighbour list by index. Bonds between dummy atoms are not filled.

Creates a vector of vectors containing bond connectivity information from the rings vector of vectors and cage information

Remove duplicate bonds.
