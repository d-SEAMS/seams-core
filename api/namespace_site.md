<!-- https://docs.dseams.info/api/namespace_site.md -->

# namespace site {#namespace-site}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Variables {#global-variables}

Relative tolerance of the inside test, in units of the largest distance from c_F to a cage vertex.

### Global Functions {#global-functions}

A trapping event is a change of [IonState](enum_site_IonState.md#doxid-namespacesite-1a7cdb003964d40a8454361db4b34b4791) into ice.

A release is a change out of ice. Other changes are recorded as well; the caller filters on `to`.

Histogram Point::{x,y,z} for typeI (0 = every atom). axis is 0, 1, or 2.

Same histogram, restricted to atoms whose [Table](struct_site_Table.md#doxid-structsite-1-1-table) kind matches.

`iceFlag` is indexed like `yCloud.pts`.

An empty `waterTypes` accepts every atom that is not an ion. A non-empty list keeps every named type, the same union `cutoffWater` uses. `cutoff` is the first shell radius in the cloud’s length unit.

`waterType` 0 accepts every non-ion. Any other value is that one type.

Rings of the water network that pass through a first shell: the census by size of every ring with at least one vertex in `shell`.

An ion is not a vertex of the network, so the rings it would have closed are gone and the rings its shell still carries measure how far the network survives around it (the hydration-shell ring census of a brine or an ionic solution). `census[s]` counts rings of size s up to `maxRingSize`.

Capped rings: the water-only primitive rings of `rings` (the ring network of the water graph, ions not vertices) of at most `maxRingSize` members, every member of which lies in `shell`.

The ion sits above a complete water ring.

Broken rings: the primitive rings, of at most `maxRingSize` members, of the water graph `nList` (by index, leading self entry, ions not bonded) plus one virtual ion vertex bonded to every water of `shell`, that contain the ion.

These are the rings the water-only network lost. Each ring is in ring order and writes the virtual ion as `static_cast<int>(nList.size())`.

`cages` are vertex index lists into `yCloud.pts`; `radius` in the cloud’s length unit (half the cage diameter, about 4 A for a 5^12 cage).

Periodic centroid of a set of atoms: every atom is unwrapped to its minimum image about the first, and the mean is taken there.

The fixed, ordered unit directions the inside test casts rays along.

The inside predicate for one cage.

The face vertices and the guest are unwrapped about c_F (`center`, or the periodic centroid of the distinct face vertices when null) and each face is fan-triangulated from its own centroid. The guest is inside when it lies on a triangle within `kInsideRelTol`; otherwise a ray is cast along the first direction of `insideRayDirections()` that passes no triangle edge or vertex within that tolerance, and the guest is inside when the ray meets an odd number of triangles. When every direction grazes, the one with the widest clearance decides.

Occupancy by the inside predicate (`pointInFaces`) on every cage.

Among cages that contain the guest the nearest centroid wins. `rings` is the ring list the cage faces index; `cageFaces[c]` is the face indices of cage c.
