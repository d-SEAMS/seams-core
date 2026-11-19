<!-- https://docs.dseams.info/api/namespace_sout.md -->

# namespace sout {#namespace-sout}

## Overview {#overview}

## Detailed Documentation {#detailed-documentation}

### Global Functions {#global-functions}

Inline function for creating the desired directory (and parents).

Uses std::filesystem for portability.

**Parameters:**

| path | The path of the directory |
| --- | --- |

**Returns:**

0 on success, 1 on failure

Function for printing out ring info, when there is no volume slice.

Function for printing out the number of prism blocks, with or without slices. Be careful when using slices!

Function for printing out the coverage area and the number of rings of each type

Function for printing out the number of rings of each type in a bulk system

Function for printing out the RDF, given the filename.

Function for printing out the number of DDCs, HCs, mixed rings, basal and prismatic rings

Function for writing out each prism.

Function for writing out cluster statistics.

Function for printing out the molecule IDs present in the slice (compatible with the LAMMPS group command

Function for printing out the molecule IDs present in the slice (compatible with the OVITO Expression Select command

Write a data file for rings.

Write out a LAMMPS dump file containing the RMSD per atom.

Write out a LAMMPS dump file containing the inSlice value for each atom for a user-defined slice

Write out a LAMMPS dump file containing the RMSD per atom for bulk ice.

Write a data file for prisms of every type.

Write a data file for rings of every type for a monolayer.

Write a data file for a particular frame, writing out topological bulk ice structures (DDCs/HCs)

Write a data file for prisms of a single type.

Write out a lammps data file for DDCs or HCs, assuming that there is no slice

Write out all cages of all types into a folder called cages inside the output directory

Write out a particular cage to a file.

Write out the basal rings of a particular Hexagonal cage.

Write out the basal rings for a particular prism.

Generic function for writing out to a dump file.

Function for printing out Q6, Cij and averaged Q3 values as single columns to text files The file names are cij, q6, q3

Function for printing the largest ice cluster.

Function for writing out the XYZ files for each cluster.
