<!-- https://docs.dseams.info/api/class_nneigh_SkinNeighborList.md -->

# class nneigh::SkinNeighborList {#class-nneigh-skinneighborlist}

## Overview {#overview}

Persistent neighbour list with a LAMMPS skin. [More…](#details-classnneigh-1-1-skin-neighbor-list)

## Detailed Documentation {#detailed-documentation}

Persistent neighbour list with a LAMMPS skin.

Vesin (or the brute-force fallback) builds candidates at cutoff+skin, which is the ghost halo: periodic images already sit in that shell. The cell list is rebuilt only when some atom has moved more than skin/2 from the last rebuild, the Verlet trigger. The bond graph for rings is the candidate pairs whose current distance is inside the analysis cutoff.

### Construction {#construction}

graph selects cutoff vs k-nearest (mutual or union). k is the neighbour count for the knn graphs (default 4).

### Methods {#methods}

Refresh from a new frame. The returned list is ID-keyed with a leading self entry, the same shape as neighListO.

Atoms whose cutoff bond set changed on the last update.
