<!-- https://docs.dseams.info/api/namespace_seams_domain.md -->

# namespace seams::domain {#namespace-seams-domain}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Variables {#global-variables}

Bonds of ghost skin around an owned atom. A primitive ring of length at most *ringSize* lies in the closed ball of radius floor(ringSize/2) about any of its members (formal/Formal/Halo.lean, primitive_ring_subset_ball), so this many bonds hold the ring and any shortcut the primitivity test can see.

Tags for the measured cell weights. Zero keeps the point-count cut.

### Global Functions {#global-functions}

The summary of *ringList* over a frame of *nAtoms* atoms.

Position of cell (x, y, z) along the Hilbert curve over a 2^bits grid.

Cutoffs of coordinates a rank needs when it builds every neighbour row itself. One bond is at most one cutoff, so the coordinate skin is [ringHaloBonds()](#doxid-namespaceseams_1_1domain_1a0a3449f33e010e244fcd5cac2d7b0127) cutoffs.

Cutoffs of neighbour rows the ring search reads past the owned block. The coordinate halo on the distributed path is one cutoff; rows out to this reach are computed by their owners and exchanged. Edges that leave the reach are dropped, which is the graph a coordinate halo of this many cutoffs would have built.

The atoms rank *rank* of *nRanks* owns, and every atom in a cell painted within *halo* of an owned cell. Without a periodic box rank 0 owns every atom.

Ranks that own a halo atom of *share*, sorted, without this rank. The relation is symmetric: a cell within *halo* of both sides is listed by each owner.

One where a local atom is owned or lies within *inner* of the owned block. *inner* below zero marks only the owned atoms.

The listed atoms of *cloud*, in that order, in the same box.

Primitive rings of the cutoff graph up to *maxDepth* members whose lowest-indexed member rank *rank* owns, by frame index, as [primitive::ringNetwork](namespace_primitive.md#doxid-namespaceprimitive-1a69f1b870183f33e58426c88ee93f0635) lists them for the graph with ascending rows. Each rank needs only its share, and the ranks’ rings partition the frame’s.

[Rings](class_seams_Rings.md#doxid-classseams-1-1-rings) of a cloud that is already one share, in share order. *owned* selects the sources. Member indices are local.

Primitive rings from neighbour rows already in hand. *frame* \[i\] is the frame index of a vertex whose neighbours are *neigh* \[i\], also frame indices. A neighbour that is not one of *frame* is dropped. *sources* lists the frames that may enumerate a ring. Members in the result are frame indices.
