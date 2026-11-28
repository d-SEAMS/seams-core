<!-- https://docs.dseams.info/api/struct_site_GuestOccupancy.md -->

# struct site::GuestOccupancy {#struct-site-guestoccupancy}

## Overview {#overview}

Guests read against enumerated cages. [More…](#details-structsite-1-1-guest-occupancy)

## Detailed Documentation {#detailed-documentation}

Guests read against enumerated cages.

A cage is a set of vertex atoms; its center is the periodic centroid of the vertices. A guest (methane, THF, an ion, any atom outside the network) sits in the cage whose center is nearest to it when that center is within `radius`; a guest with no center within `radius` is free. This is the occupancy of a clathrate hydrate (fraction of 5^12 and 5^12 6^2 cages filled) and the test of whether an ion sits inside the network rather than at a vertex.

### Fields {#fields}

per cage, guests inside

per guest in input order, cage index or -1

per guest, distance to its cage center (-1 when free)

cages with at least one guest

cages with more than one guest

guests in no cage

hist\[k\] = cages with k guests
