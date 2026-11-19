<!-- https://docs.dseams.info/api/class_seams_Neighbours.md -->

# class seams::Neighbours {#class-seams-neighbours}

## Overview {#overview}

Owned neighbour list. [More…](#details-classseams-1-1-neighbours)

### Inherited Members {#inherited-members}

## Detailed Documentation {#detailed-documentation}

Owned neighbour list.

Move-only. Destruction frees the C handle. The self entry is omitted. A neighbour list is not safe to share across threads.

### Methods {#methods}

Number of neighbour pairs, omitting self.

**Returns:**

The pair count.

Degree of each atom.

**Returns:**

One integer per atom. The self entry is not counted.

Directed pairs `(atom, neighbour)`.

**Returns:**

One pair per stored column. Self is omitted.

Chordless cycles of one length on this graph.

**Parameters:**

| size | Cycle length, from 3 through 12. |
| --- | --- |

**Returns:**

The rows and the per-atom counts of seams_rings_chordless.
