<!-- https://docs.dseams.info/api/struct_seams_domain_RingSummary.md -->

# struct seams::domain::RingSummary {#struct-seams-domain-ringsummary}

## Overview {#overview}

What the ring analyses read off a frame’s primitive rings. [More…](#details-structseams-1-1domain-1-1-ring-summary)

## Detailed Documentation {#detailed-documentation}

What the ring analyses read off a frame’s primitive rings.

### Fields {#fields}

[Rings](class_seams_Rings.md#doxid-classseams-1-1-rings) of each size, indexed by member count, 0 through maxDepth.

Largest ring size each atom belongs to, 1 when it is in none, the rule [ring::assignPolygonType](namespace_ring.md#doxid-namespacering-1afa6c593ce7cd865c74e11d39fa27f04b) applies. Filled on the root rank only.
