<!-- https://docs.dseams.info/api/class_primitive_RingUpdater.md -->

# class primitive::RingUpdater {#class-primitive-ringupdater}

## Overview {#overview}

Exact frame-to-frame maintenance of the primitive ring network. [More…](#details-classprimitive-1-1-ring-updater)

## Detailed Documentation {#detailed-documentation}

Exact frame-to-frame maintenance of the primitive ring network.

Consecutive trajectory frames share almost all of their bond topology. A ring of at most maxDepth members can appear or vanish between frames only if it passes within (maxDepth/2 - 2)/2 bonds of a vertex whose neighbour row changed, which is on the ring itself up to maxDepth

1. On a new frame the stored rings that meet that set are dropped and the rings through it are enumerated afresh. The output equals ringNetwork element by element, at a cost that scales with the number of changed rows rather than with system size.

### Methods {#methods}

Primitive rings for this frame’s neighbour list (row-ordered, by index, first element of each row the vertex itself)

Vertices the last update enumerated rings through; the full vertex count on a whole-frame pass, zero when nothing changed

Vertices whose bounded balls were rebuilt on the last update.

Undirected bonds whose endpoints differ from the previous frame. Zero on the first frame and when the bond set is unchanged. The fraction below is this count divided by the molecule count, maintained while the rows are compared and readable in constant time afterwards.

Changed bonds per molecule on the last update.

Kept for callers that used to force a whole-frame pass. The update enumerates the sources in T and no others, so a whole-frame pass happens only when T is every vertex. The stored fraction does not replace that search. The result equals a rebuild at every fraction.

Kept for callers that used to force a whole-frame pass from the changed-bond rate. The default is negative. The rate does not replace the search of T.
