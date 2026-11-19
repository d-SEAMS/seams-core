<!-- https://docs.dseams.info/api/namespace_molSys.md -->

# namespace molSys {#namespace-molsys}

## Overview {#overview}

Bare-bones structs used throughout the architecture. [More…](#details-namespacemol-sys)

## Detailed Documentation {#detailed-documentation}

Bare-bones structs used throughout the architecture.

@detials This namespace defines [Point](struct_molSys_Point.md#doxid-structmol-sys-1-1-point) and [PointCloud](struct_molSys_PointCloud.md#doxid-structmol-sys-1-1-point-cloud) structs, alongwith other basic functions and enums.

[PointCloud](struct_molSys_PointCloud.md#doxid-structmol-sys-1-1-point-cloud) is a struct that contains the information of a collection of [Point](struct_molSys_Point.md#doxid-structmol-sys-1-1-point) structs, at every frame. [Point](struct_molSys_Point.md#doxid-structmol-sys-1-1-point) contains the coordinates and types of each point, along with information about whether the particular particle is within a slice or not.

### Changelog {#changelog}

- Amrita Goswami \[[amrita16thaug646@gmail.com](mailto:amrita16thaug646@gmail.com)\]; date modified: Sept 19, 2019
- Rohit Goswami \[[rog32@hi.is](mailto:rog32@hi.is)\]; date modified: Mar 20, 2021

#### Global Functions {#global-functions}

Creates an unordered map, with the atomIDs as keys and molecular IDs as the values

Creates an multimap with molecule IDs of the atoms as the keys and the atom IDs as the values. More than one atom can have the same molecule ID

Returns a vector of vectors, which contains the molIDs in the first column, and the hydrogen atom indices (not atom IDs) in the row

This function searches a vector of vectors molList, for a particular molecular ID, and returns the index in molList

//! Function for clearing vectors in [PointCloud](struct_molSys_PointCloud.md#doxid-structmol-sys-1-1-point-cloud) after multiple usage
