<!-- https://docs.dseams.info/api/struct_ring_RingSearchIndex.md -->

# struct ring::RingSearchIndex {#struct-ring-ringsearchindex}

## Overview {#overview}

Inverted index from an atom to the rings that contain it. [More…](#details-structring-1-1-ring-search-index)

## Detailed Documentation {#detailed-documentation}

Inverted index from an atom to the rings that contain it.

The cage searches repeatedly ask which rings a given atom belongs to. Answering that by scanning every ring makes the searches quadratic in the ring count. One pass over the ring network answers it in constant time thereafter, and the ring indices in each row stay ascending, which is the order the searches previously visited them in.

findHC visits rings in the order of a key that does not use absolute coordinates: the sorted atom indices, then the unwrapped shift of the cycle. findDDC reports a maximum set of cages under the same key, so the cage count is independent of the order in which rings are claimed and does not change when the frame is translated.

### Fields {#fields}

Same rows as ringsContainingAtom, packed. ringCol\[ringOff\[a\] .. ringOff\[a+1\]) is the ascending ring list of atom a. Empty when the index was not built.

Integer minimum-image winding of each ring. Empty if no box was supplied.

### Methods {#methods}

Atoms the packed rows cover.
