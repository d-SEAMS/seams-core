<!-- https://docs.dseams.info/api/namespace_primitive.md -->

# namespace primitive {#namespace-primitive}

## Overview {#overview}

Functions for generating primitive rings. [More…](#details-namespaceprimitive)

## Detailed Documentation {#detailed-documentation}

Functions for generating primitive rings.

This namespace contains struct definitions and functions that are used for generating primitive (shortest-path) rings (directed cyclic graphs).

The [Vertex](struct_primitive_Vertex.md#doxid-structprimitive-1-1-vertex) object is a collection of elements for each point, required for graph traversal. The [Graph](struct_primitive_Graph.md#doxid-structprimitive-1-1-graph) object is an object for the whole frame, containing the information of all vertices, and a row-ordered vector of vector of the rings generated.

The [Franzblau shortest-path criterion](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.44.4925) has been used. The SP (shortest-path) criterion is midway between the least restrictive and most restrictive criteria in the hierarchy.

The following is the procedure for finding primitive rings:

1. All possible rings (including non-SP) rings are found, in the [primitive::countAllRingsFromIndex](#doxid-namespaceprimitive-1ac81b1dd4b67e8200de61e52c498c1b2f) function, using the [backtracking algorithm](https://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf). This is a recursive algorithm.
2. The non-SP rings are then removed from the list of all rings, using the Franzblau shortest path criterion ([primitive::removeNonSPrings](#doxid-namespaceprimitive-1a592086d3afae4fab5aa5ff6aa5c1d215)), answered by a hop-bounded breadth-first sweep per vertex.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Nov 14, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Returns a vector of vectors containing the rings (of all sizes), by atom index, given the neighbour list also by index (preferably the hydrogen-bonded neighbour list). Each ring starts at its lowest-indexed member and runs toward the lesser of that member’s ring neighbours, and the rings are sorted lexicographically.

Per-atom primitive-ring counts, without storing the rings. `counts[atom * (maxDepth + 1) + size]` is how many rings of that size contain the atom.

Primitive rings whose size is at least *minDepth* and at most *maxDepth*. ringNetwork(nList, maxDepth) is minDepth 3. A separate name keeps ringNetwork(nList, maxDepth, {}) the empty source mask.

The rings ringNetwork returns whose lowest-indexed member is flagged in *sources* (by index; empty flags every vertex), in the same order

Rings of a CSR neighbour graph. `cols[offsets[i] .. offsets[i + 1]]` are the neighbours of vertex i, with no self entry. The ring set is the one ringNetwork returns for those rows with the vertex written in front of each row.

Primitive rings of an explicit caller graph after converting every directed edge to one undirected edge. Self edges and duplicate columns are ignored. Throws std::invalid_argument for malformed CSR or an out-of-range column. The search itself is ringNetwork(offsets, cols, maxDepth).

Histogram of a ring list. counts\[s\] is the number of size-s rings. The returned vector has maxDepth + 1 entries.

Ring sizes of the same search as ringNetwork(offsets, cols, maxDepth). counts\[s\] is the number of primitive rings of size s, for s below countsLen. The rings themselves are not gathered.

The same ring set, built with a bounded ball at every vertex and a primitivity test of every pair. ringNetwork leaves that table out when the largest ring is at most 7, because a strict shortcut is then a walk of at most two edges. The timing harness compares the two.

0 when the last ringNetwork call answered every primitivity query from the adjacency rows, 1 when it filled the sorted pair table, 2 when it filled the occupancy bitset of a radius-4 ball, 3 when a radius-2 pair table answers distance 4 by a meet in the middle. Rings of size at most 7 report 0. Rings of size 8 and 9 report 1. Rings of size 10 report 3. ringNetworkBalls always reports 1.

Creates a graph object and fills it with the information from a neighbour list and pointCloud created before. NOTE: the neighbourListIndex contains the indices and NOT the atom IDs as in the neighbour list

Creates a graph object and fills it with the information from a neighbour list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex contains the indices and NOT the atom IDs as in the neighbour list

Re-fills the neighbour lists of a graph object from a neighbour list of INDICES NOT ATOM IDs created before. NOTE: the neighbourListIndex contains the indices and NOT the atom IDs as in the neighbour list

Creates a vector of vectors of all possible rings.

Removes the non-SP rings, using the Franzblau shortest path criterion.

Main function that searches for all rings.

Calculates the shortest path.

Function for clearing vectors in [Graph](struct_primitive_Graph.md#doxid-structprimitive-1-1-graph) after multiple usage.

Balls, source enumeration and the gather from the last ringNetwork call.

When set, ringNetwork records the primitive-test time in the split.

Primitive-test milliseconds from the last timed ringNetwork, then zero.

When non-null, ringNetwork adds each source’s enumeration time, in microseconds, at that source’s index. The buffer must hold one entry per neighbour-list row. A null buffer skips the clock.
