<!-- https://docs.dseams.info/api/struct_chill_BondOrder.md -->

# struct chill::BondOrder {#struct-chill-bondorder}

## Overview {#overview}

Unnormalised real l = 3 bond-order vector of the atoms *which* selects. [More…](#details-structchill-1-1-bond-order)

## Detailed Documentation {#detailed-documentation}

Unnormalised real l = 3 bond-order vector of the atoms *which* selects. *which* null computes every atom. Ghost slots stay zero so a neighbour rank can write the vector it computed for that atom.

### Fields {#fields}

atom-major, 7 channels

squared norm of q, with the m != 0 factor

keepMax > 0, local indices, -1 empty
