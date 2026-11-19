<!-- https://docs.dseams.info/api/namespace_pntToPnt.md -->

# namespace pntToPnt {#namespace-pnttopnt}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Functions {#global-functions}

Fills up an eigen matrix point set a reference ring, which is a regular n-gonal polygon; where n is the number of nodes in the ring

Creates an eigen matrix for the points of a prism block, constructed from the points of a perfect polygon of radius 1, given the basal rings and axial dimension

Calculate the average radial distance for the basal rings, calculated from the centroid of each basal ring

Calculate the average height of the prism block, calculated using the basal rings of the prism and the axial dimension

Get the relative ordering of a pair of basal rings for a deformed prism/perfect prism. Outputs a vector of vectors of indices, such that the first vector is for the first basal ring, and the second vector is for the second basal ring. The input neighbour list is with respect to indices, not IDs

Get the relative ordering of a pair of basal rings for a deformed prism/perfect prism. Outputs a vector of vectors of indices, such that the first vector is for the first basal ring, and the second vector is for the second basal ring.

Fill up an Eigen Matrix of a prism basal ring from an input vector of atom indices

Fill up an Eigen matrix of a prism block (two basal rings) from input vectors for the basal rings

Fills up an eigen matrix point set a reference ring, which is a regular n-gonal polygon, constructed with radius 1 by default; where n is the number of nodes in the ring

Matches the order of the basal rings of an HC or a potential HC.

Matches the order of the basal rings of an DDC or a potential HC.

As relOrderDDC, but deriving the connectivity from an explicit equatorial traversal instead of the stored ring row; passing the reversed traversal yields the mirror correspondence for the shape-matching search

Fills up an eigen matrix point set using the basal rings basal1 and basal2, changing the order of the point set by filling up from the startingIndex (starting from 0 to 5)

Fills up an eigen matrix point set using information of the equatorial ring and peripheral rings, embedded in a vector, already filled in relOrderDDC.
