<!-- https://docs.dseams.info/api/class_ring_AffiliationUpdater.md -->

# class ring::AffiliationUpdater {#class-ring-affiliationupdater}

## Overview {#overview}

Exact incremental affiliation across frames. [More…](#details-classring-1-1-affiliation-updater)

## Detailed Documentation {#detailed-documentation}

Exact incremental affiliation across frames.

Stores, per six-ring, its passing basal pairs with their prismatic rings and its equatorial result. On a new frame it matches the rings to the previous frame’s and evaluates again only the pairs of rings that are new, hold an atom whose neighbour row changed, or hold or bond an atom of an added or removed ring, and the equatorial test of rings that are new or hold an atom of an added or removed ring. The flags are the unions over the stored answers, identical to cageAffiliation on the new frame.

### Methods {#methods}

Classify the frame, reusing the previous frame’s answers outside the locality closure of what changed. The reference stays valid until the next update call.

Rings whose pairs or equatorial test the last update evaluated again; the first call reports every ring.

Above this fraction of dirty rings, update runs one sweep instead of remapping stored pairs. The default is 1, so a transition keeps every ring outside the changed region. A smaller fraction opts into the sweep. An unchanged graph returns the stored flags.
