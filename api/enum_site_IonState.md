<!-- https://docs.dseams.info/api/enum_site_IonState.md -->

# enum site::IonState {#enum-site-ionstate}

## Overview {#overview}

Ions read against a per-atom ice assignment. [More…](#details-namespacesite-1a7cdb003964d40a8454361db4b34b4791)

## Detailed Documentation {#detailed-documentation}

Ions read against a per-atom ice assignment.

Ions are not part of the hydrogen-bond network; the assignment is computed on the water and each ion is classed by its first water shell: every shell molecule labeled is `ice`, none is `liquid`, otherwise `front`. An ion with an empty shell is liquid.
