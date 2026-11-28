<!-- https://docs.dseams.info/api/struct_seams_hbond_edge.md -->

# struct seams_hbond_edge {#struct-seams-hbond-edge}

## Overview {#overview}

Hydrogen-bond edge record for one donor-acceptor pair. [More…](#details-structseams-hbond-edge)

## Detailed Documentation {#detailed-documentation}

Hydrogen-bond edge record for one donor-acceptor pair.

### Fields {#fields}

Caller-set sizeof(seams_hbond_edge) on every element.

Donor heavy-atom index in the frame.

Acceptor heavy-atom index in the frame.

Donor hydrogen index, or UINT32_MAX when omitted.

Donor-acceptor distance in Angstrom, or NaN when omitted.

O-O-H angle in degrees, or NaN when omitted.
