<!-- https://docs.dseams.info/api/namespace_chordless.md -->

# namespace chordless {#namespace-chordless}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Variables {#global-variables}

Largest cycle the meet-in-the-middle search takes. A half is at most six steps.

### Global Functions {#global-functions}

Chordless cycles of one size on a CSR graph.

A cycle is chordless when no edge joins two of its vertices that are not neighbours on the cycle. Each cycle is listed once, starting at its lowest vertex, in the direction whose second vertex is the smaller of that vertex’s two cycle neighbours. The list is sorted lexicographically, so it does not depend on the thread count.

**Parameters:**

| offsets | n + 1 row offsets. |
| --- | --- |
| cols | Column indices. Self entries must not appear. |
| size | Cycle length, from 3 through kMaxSize. |
| rings | When not null, receives count * size vertex indices. |
| per_atom | When not null, n counts of cycles through each vertex. |

**Returns:**

The number of cycles.
