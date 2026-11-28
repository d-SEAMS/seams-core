<!-- https://docs.dseams.info/api/struct_seams_DirectedBond.md -->

# struct seams::DirectedBond {#struct-seams-directedbond}

## Overview {#overview}

One directed hydrogen bond. [More…](#details-structseams-1-1-directed-bond)

## Detailed Documentation {#detailed-documentation}

One directed hydrogen bond.

Indices are frame indices. Heavies from set_hydrogen_frame occupy `[0, n_heavy)` and hydrogens follow. `angle` is acceptor-first degrees: 0 is a linear donor-hydrogen-acceptor geometry.

### Fields {#fields}

Donor heavy-atom index.

Hydrogen index.

Acceptor heavy-atom index.

Acceptor-hydrogen distance, Angstrom.

Donor-acceptor distance, Angstrom.

Acceptor-first angle, degrees.
