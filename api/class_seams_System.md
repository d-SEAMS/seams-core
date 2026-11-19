<!-- https://docs.dseams.info/api/class_seams_System.md -->

# class seams::System {#class-seams-system}

## Overview {#overview}

One frame: coordinates, a cell, and the analyses that read them. [More…](#details-classseams-1-1-system)

### Inherited Members {#inherited-members}

## Detailed Documentation {#detailed-documentation}

One frame: coordinates, a cell, and the analyses that read them.

Move-only. The destructor frees the C handle. Methods that read the frame require exclusive use of this object. The [Steinhardt](struct_seams_Steinhardt.md#doxid-structseams-1-1-steinhardt) scratch inside the library is process-global, so two frames must not run `steinhardt`, `ql`, `q4q6`, or `q6bar` at the same time. `set_num_threads` is also process-global.

### Construction {#construction}

Allocate an empty frame.

**Parameters:**

| Error | with status Memory when the allocation fails. |
| --- | --- |

### Methods {#methods}

Store points in an orthorhombic cell whose origin is zero.

**Parameters:**

| xyz | Interleaved x, y, z. Length is a multiple of 3. Angstrom. |
| --- | --- |
| lx | Edge length along x, Angstrom. |
| ly | Edge length along y, Angstrom. |
| lz | Edge length along z, Angstrom. |

Store points in an orthorhombic cell.

**Parameters:**

| xyz | Interleaved x, y, z. Length is a multiple of 3. Angstrom. |
| --- | --- |
| origin | Cell origin, Angstrom. |
| lengths | Edge lengths, Angstrom. Each must be positive. |
| types | One species per atom, or an empty span. |
| ids | One id per atom, or an empty span. |

Store points and a general nine-number cell.

**Parameters:**

| xyz | Interleaved coordinates, Angstrom. |
| --- | --- |
| box9 | ox, oy, oz, lx, ly, lz, xy, xz, yz. |
| types | One species per atom, or empty. |
| ids | One id per atom, or empty. |

Read one frame from a trajectory path.

**Parameters:**

| path | XYZ, CON, chemfiles, or a LAMMPS dump. |
| --- | --- |
| frame | 1-based frame index. |
| type_filter | Species to keep. Negative keeps every atom. |

Number of points in the frame.

**Returns:**

The atom count. An empty frame returns 0.

Original CHILL labels and nine population counts.

CHILL+ labels and nine population counts.

Local q_l on the cutoff shell.

**Parameters:**

| options | Cutoff, species, degree, and whether to average. |
| --- | --- |

Local q_l on a named graph.

Local q4 and q6 on the cutoff shell.

q6-bar, its mean and median, and the ice count above the cut.

Cutoff or mutual neighbour list.

Chordless cycles of one length on a neighbour graph.

**Parameters:**

| options | Graph of the cycles. |
| --- | --- |
| size | Cycle length, from 3 through 12. |

**Returns:**

[Neighbours::chordless](class_seams_Neighbours.md#doxid-classseams-1-1-neighbours-1a13edfc37618cf2779867c27eb72ed1fd) on that graph.

k-nearest neighbour list.

Degree of each atom on the cutoff graph.

**Parameters:**

| cutoff | Positive cutoff, Angstrom. |
| --- | --- |
| species | Species filter. Zero keeps every atom. |

Primitive rings up to a maximum size.

Hexagonal and double-diamond cage membership.

Closed cages of one ring census.

Partial radial distribution and the cell volume.

Running coordination number along the partial RDF.

Coordination number integrated to rmax.

Hydrogen-bond count.

Number density along one cell axis.

Topological fingerprint of the bonded graph.

Mutual nearest unlike contacts.

`spec` is site-table text.

Largest polar or apolar domain.

Perfect and imperfect prism counts by ring size.

Primitive-ring histogram of a monolayer.

**Returns:**

`counts[k]` is the number of rings of size k.

Per-atom site kind.

**Parameters:**

| spec | Site-table text, for example 1=waterO. |
| --- | --- |

**Returns:**

One kind code per atom.

Per-atom ice-shell label.

**Parameters:**

| shell_rmsd | Positive residual ceiling, or a non-positive value for the library ceiling. |
| --- | --- |

**Returns:**

0 ice Ih, 1 ice Ic, 2 stacking disorder, or -1 when the shell does not match.

Planar ice front from an ice assignment.

**Parameters:**

| ice | Non-zero when the atom is ice. Length must equal count(). |
| --- | --- |
| options | Water species, ion species, bins, and axis. |
