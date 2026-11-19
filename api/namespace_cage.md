<!-- https://docs.dseams.info/api/namespace_cage.md -->

# namespace cage {#namespace-cage}

## Overview {#overview}

Functions for topological network criteria cage types. [More…](#details-namespacecage)

## Detailed Documentation {#detailed-documentation}

Functions for topological network criteria cage types.

This namespace contains structs and enums for cage types

Type definitions for atoms, rings and cages. HexC and DoubleDiaC are the TUM ice cages. A [Signature](struct_cage_Signature.md#doxid-structcage-1-1-signature) names any polyhedron by its ring-size census (sodalite is {4:6, 6:8}); the enumerator in `cage_enum.hpp` grows face-sharing ring sets that match that census.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Nov 13, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

True when every edge of the listed faces is used by exactly two of those faces.

An empty face list is not closed.

True when the listed faces are closed (isClosedPolyhedron), at least three faces around each of their vertices form a single cycle, and V - E + F = 2.

Paired edges alone also admit pinched vertices, two faces that share two edges, and surfaces of higher genus.

Face-sharing rings whose size census equals `signature` and whose edges close.

The input may hold rings of every size; only sizes that appear in the signature take part. Named `hc` and `ddc` without a neighbour list use the geometric census (`4:6,6:2` and `6:7`).

True when `cycle` is `pattern` rotated or reversed.

As findBySignature.

A ring is a face only when its species sequence matches one pattern of the same length, up to rotation and reversal. An empty pattern list keeps every ring whose size is in the census. `species` is one class per atom index. A vertex past the end of `species` keeps its ring out.

As above.

Named `hc` and `ddc` call findHC / findDDC on the six-membered rings so the vertex sets match those finders.

As above, with findHC / findDDC claiming rings in the order of `atomRank` ([ring::positionRank](namespace_ring.md#doxid-namespacering-1a51764966eab7aea2d90649bb6cc607ea) of the frame), so the `hc` and `ddc` cages do not depend on the atom labels.

`none`, `off`, or `false` disables the window.

`LO,HI`, `LO:HI`, or `LO/HI` sets the inclusive bounds. An empty string is the default window. `lo` is >= 0 and strictly less than `hi`.

`none` when the window is off, otherwise `lo,hi`.

findBySignature for each of `signatures`, one list per entry.

The ring-size censuses share one search, so a face set that fits several of them is grown once; `hc` and `ddc` go to findHC / findDDC as in findBySignature. Each list equals findBySignature on its entry. `window` keeps a face only when every minimum-image edge lies in \[lo, hi\]. A disabled window, or a null `xyz`, leaves that test out and the cage bound is the face census. `box` is the cloud box (three lengths, or lengths plus xy, xz, yz). `boxLow` is the dump origin; a tilted cell uses it so a wrapped oxygen-oxygen contact is measured on the lattice, not on the orthorhombic bounds.

Connected face sets that stay inside the signature budget, have at least `minFaces` faces, and are not closed.

Cups and incomplete cages during hydrate nucleation. Closed polyhedra are omitted. `minFaces <= 0` uses `max(1, signature.faceCount() / 2)`. A cup whose every face already belongs to a closed cage is omitted.

Runs Search with the closed seed filter left on.

A ring whose edge sits in fewer than two rings is skipped.

Coordinates of `vertices` only.

An index past `all` is skipped.

Point group of this vertex set.

Does not see any other atom.

Overlay `cageXyz` on `ref`.

Row counts must agree. Does not see any atom outside the two sets.

`formerSpecies < 0` keeps every atom.

True when the two tables have the same per-atom species, coordination, homopolar count, and ring-size multiset.
