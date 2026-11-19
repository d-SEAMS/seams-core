<!-- https://docs.dseams.info/api/struct_molSys_Point.md -->

# template struct molSys::Point {#template-struct-molsys-point}

## Overview {#overview}

This contains per-particle information. [More…](#details-structmol-sys-1-1-point)

## Detailed Documentation {#detailed-documentation}

This contains per-particle information.

Specifically

- Type ID
- Molecular ID
- Atom ID
- Coordinates
- Neighbourlist
- Bond correlation type of type [Result](struct_molSys_Result.md#doxid-structmol-sys-1-1-result)
- Type of [atom_state_type](enum_molSys_atom_state_type.md#doxid-namespacemol-sys-1a13ffe3a39a0816627553f195e102abb5) iceType
- In slice bool

### Fields {#fields}

type ID, molID, atomID

coordinates

Results (contains bond correlation type).

Type of ice/water etc based on cij.
