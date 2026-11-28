<!-- https://docs.dseams.info/api/group_capi_pore.md -->

# Pore geometry {#pore-geometry}

## Overview {#overview}

Fractional density, cylindrical pore profiles, and host-guest contacts. [More…](#details-group-capi-pore)

## Detailed Documentation {#detailed-documentation}

Fractional density, cylindrical pore profiles, and host-guest contacts.

Slab and shell volumes use the triclinic cell, det(H). A fractional bin along a cell edge has volume V/n. An arbitrary direction bins the projection of the reduced cell, and each bin volume is the measure of that slice. Number density is atoms per cubic Angstrom. Mass density is the caller’s mass unit per cubic Angstrom, or mass 1 when masses is NULL. A mass given in amu is therefore amu per cubic Angstrom, and multiplying by 1.66054 converts that density to grams per cubic centimetre. The library does not apply the factor.
