<!-- https://docs.dseams.info/api/namespace_water.md -->

# namespace water {#namespace-water}

## Overview {#overview}

Water-oxygen selection for a mixed frame. [More…](#details-namespacewater)

## Detailed Documentation {#detailed-documentation}

Water-oxygen selection for a mixed frame.

Neighbours, CHILL, CHILL+, rings and cages run on the oxygen atoms. Ions and hydrogens stay out of that graph. A requested type or molecule that is not in the frame is an error, as is a selection that keeps nothing. The messages are the strings the tests match.

### Global Functions {#global-functions}

Historical dump guess: type 2 when any atom has it, otherwise type 1.

A miss leaves an empty cloud and keeps the box.

oxygen is the labeled water subset.

front is the nearest ice oxygen, in the box units, or -1 when no ion sits next to ice.
