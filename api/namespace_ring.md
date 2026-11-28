<!-- https://docs.dseams.info/api/namespace_ring.md -->

# namespace ring {#namespace-ring}

## Overview {#overview}

Topological network criteria functions. [More…](#details-namespacering)

## Detailed Documentation {#detailed-documentation}

Topological network criteria functions.

This namespace contains functions for the topological network criteria for bulk and confined systems both.

Although the namespace is shared by bulk and confined topological network criteria, the functions are split into files specific to each set of criteria.

All the topological network criteria are based on the identification of primitive rings (using the Franzblau algorithm for primitive rings), following which specific rules for combinations of primitive rings are used to identify quasi-one-dimensional, quasi-two-dimensional and bulk ices.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Nov 13, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Batch computation of the affiliation predicates for six-membered rings. nList is by index with the leading self entry, as neighbourListByIndex produces.

Seeded affiliation over two graphs on the same atoms: the strict graph (typically the mutual k-nearest bonds) supplies seeds, the permissive supergraph (typically the union bonds) supplies completions, and a permissively affiliated atom is accepted only when its bonded component of affiliated atoms contains a seed.

Specificity on structureless input is structural rather than statistical: when the strict pass affiliates nothing, nothing is accepted regardless of what the permissive graph builds. Where both graphs label an atom, the strict labels win.

The fifth argument is the ring-adjacent completion flag. The default is false. When the flag is true the accepted labels pass through [ringAdjacentCompletion()](#doxid-namespacering-1a52d5bc312612884ce62bd430d94fdb17) on the permissive six-rings. HC and DDC complete on separate flag vectors. seams cages complete sets the flag. tests/walk_compare leaves the flag false.

**Parameters:**

| strictRings | Six-membered rings of the strict graph. |
| --- | --- |
| strictNList | Strict graph, by index with leading self entries. |
| permissiveRings | Six-membered rings of the permissive graph. |
| permissiveNList | Permissive graph, same conventions. |
| ringAdjacentCompletion | When true, call ringAdjacentCompletion() on the accepted labels. Default false. |

The same labels from the affiliation of each graph’s six-rings, computed by the caller.

A frame walk passes [AffiliationUpdater](class_ring_AffiliationUpdater.md#doxid-classring-1-1-affiliation-updater) results here, so both graphs carry their per-ring answers across frames.

**Parameters:**

| strictRings | Six-membered rings of the strict graph. |
| --- | --- |
| strictAffiliation | cageAffiliation of strictRings. |
| permissiveRings | Six-membered rings of the permissive graph. |
| permissiveAffiliation | cageAffiliation of permissiveRings. |
| permissiveNList | Permissive graph, by index with leading self entries. |
| ringAdjacentCompletion | As in the overload above. |

Ring completion of seeded labels.

A permissive six-ring whose vertices all carry a cage label but one is a cage ring with a vacancy. The last vertex takes that label. The walk repeats until a fixed point. HC and DDC complete on separate flag vectors.

A liquid ring that only touches a nucleus has at most a few labeled vertices, so the all-but-one rule cannot walk into the liquid. An empty seed stays empty: no ring then has five labeled vertices.

The result is the least fixed point above the seed. Visiting order does not change the labeled set. lean/DseamsProofs/Completion.lean states that claim. Catch2 covers the C++ identity on a crystal, the refill of a single vacancy, and the structural zero.

**Parameters:**

| labels | Seeded per-atom flags. |
| --- | --- |
| permissiveRings | Six-membered rings of the permissive graph. |

Returns a vector of vectors of rings of a single size.

Check to see if two vectors have common elements or not True, if common elements are present and false if there are no common elements

Compares two disordered vectors and checks to see if they contain the same elements

Searches a particular ring for a triplet.

Common elements in 3 rings.

Returns the common elements of two rings.

Erases memory for a vector of vectors for a list of rings.

Assign an atomType (equal to the number of nodes in the ring) given n-membered rings.

Select edge molecules and atoms which are part of rings, such that rings formed with even one atom in the slice will be included in the selection Modifies the inSlice bool of a given PointCloud (this may be the same) as the given oxygen atom PointCloud which was used to construct the neighbour list used to construct the rings vector of vectors. We assume that the PointCloud structs have the inSlice bool values set according to the presence of the atom in the slice (this can be done using the [gen::moleculesInSingleSlice](namespace_gen.md#doxid-namespacegen-1aa6a3bfc360df06a888d89bbd1496f4c7) function.

Master function for selecting edge molecules and atoms which are part of rings, such that rings formed with even one atom in the slice will be included in the selection Modifies the inSlice bool of a given PointCloud (this may be the same) as the given oxygen atom PointCloud which was used to construct the neighbour list used to construct the rings vector of vectors (calls [ring::getEdgeMoleculesInRings](#doxid-namespacering-1af72db90faa20a1dfdfaf8c0e868b305d)) Prints out molecule IDs individually of molecules in the slice, and also prints out a LAMMPS data file of just the molecules and atoms in the slice

Find out rings in the bulk, looping through all ring sizes upto the maxDepth The input ringsAllSizes array has rings of every size.

Find out which rings are DDCs or HCs, which are comprised of 6-membered primitive rings. Start with a neighbour list (by index) and a vector of vectors of rings (also by index). TODO: try ‘square’ ice and ice0

Builds the inverted atom-to-rings index used by the cage searches. With *withRows* false only the packed rows are filled, and ringsContainingAtom stays empty.

As buildRingSearchIndex, ranking the atoms by position ([ring::positionRank](#doxid-namespacering-1a51764966eab7aea2d90649bb6cc607ea)).

Rank of every atom of *yCloud* when the atoms are sorted by position wrapped into the box, compared by x, then y, then z. Throws std::invalid_argument when two atoms share a wrapped position

Oriented integer winding of every ring. Component k is the sum of the minimum-image bond steps, traversed in ring order, in units of periodic cell vector k: a ring that crosses the cell once along +x is +1 on axis

1. A contractible ring has {0, 0, 0}. Throws std::invalid_argument when a vertex is outside the cloud. A cloud without a valid box returns zeros.

The rings of *rings* whose minimum-image bonds have integer lattice images summing to zero. Any other ring winds through the periodic cell and bounds no face. A cell narrower than a ring’s span has such rings. Without a box every ring is kept

Same rings, compacted in place, so a kept ring is moved and not copied.

Faces for a cage search: contractibleRings, after refusing a cutoff that can contain more than one periodic image of a pair. The narrowest face separation must be longer than twice *cutoff*. A cloud without a box keeps every ring. Throws std::invalid_argument when the cell is too thin for *cutoff*

Find out which hexagonal rings are DDC (Double Diamond Cages) rings. Returns a vector containing all the ring IDs which are DDC rings

Peripheral rings of one equatorial candidate. n == 0 when the ring is not a candidate. The first six sit inline. A defect that names more conflicts stores the rest beside the row, so a crystal of candidates stays a plain block of integers.

Maximum set of mutually non-conflicting equatorial rings. *periph* \[i\].n is 0 when ring i is not a candidate. The returned ids are a maximum packing; *certificate* is filled when non-null.

As findDDC, reusing an index the caller has already built. When *certificate* is given it receives the matching or search bound.

Find out which hexagonal rings are both DDCs (Double Diamond Cages) and HCs (Hexagonal Cages). Returns a vector containing all the ring IDs which are of this type

Find out which hexagonal rings are HC rings. Returns a vector containing all the ring IDs which are HC rings

As findHC, reusing an index the caller has already built.

First condition for the DDC: There must be at least 3 other rings in which each element of the equatorial ring is present

As conditionOneDDC, reusing an index the caller has already built.

Equatorial test of a DDC. Holds when every atom of ring iring lies in at least three other rings and, for some starting atom and direction of iring, six distinct rings P_k through its consecutive triplets satisfy the even/odd common-atom and alternate-pair conditions. When *peripheralUnion* is given it receives every ring that serves as a P_k in some valid assignment, sorted. When *chosen* and *rank* are given, *chosen* receives the six rings of the valid assignment whose ranks, sorted, are lexicographically least, in ascending rank.

Tests whether two rings are basal rings (true) or not (false): some alternate triplet of basal1, starting at its first or second atom, is bonded one-to-one to an alternate triplet of basal2, and no atom of the other triplet of basal2 is bonded to basal1

As basalConditions, answering bonds from a prebuilt [NeighBond](class_ring_NeighBond.md#doxid-classring-1-1-neigh-bond).

Tests whether the last two elements of a triplet are neighbours of two atom IDs passed in

Tests to check that elements of a triplet are not neighbours of a ring (vector) passed

Finds the prismatic rings from basal rings iring and jring.

As findPrismatic, drawing candidate rings from the inverted index.

Assigns a type of enum class iceType, to every atom, using information from ringType, which has the information of every ring

Determines the number of HCs, DDCs, Mixed rings, prismatic and basal rings.

Conflict build, clique ids, matching and branch-and-bound from the last pack.

Order, equatorial test and packing from the last findDDC call.

Find out which rings are prisms. Returns a vector containing all the ring IDs which are prisms

Tests whether two rings are basal rings (true) or not (false) for a prism (strict criterion)

Reduced criterion: Two candidate basal rings of a prism block should have at least one bond between them

Checks whether two 4-membered rings are parallel in one dimension or not to prevent overcounting

Find out which rings are prisms, looping through all ring sizes upto the maxDepth The input ringsAllSizes array has rings of every size.

Assign an atomType (equal to the number of nodes in the ring) given a vector with a list of indices of rings comprising the prisms

Get the atom type values for deformed prisms.

Shift the entire ice nanotube and remove axial translations.

Find out which rings are prisms, looping through all ring sizes upto the maxDepth The input ringsAllSizes array has rings of every size.
