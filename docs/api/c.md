# C API

`seams.h` is ABI 1 of libseams. Handles are opaque. Coordinates are Angstrom.
A cell is nine numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`. A zero tilt is orthorhombic.
`lx, ly, lz` are edge lengths.

## Version and errors

`seams_abi` returns 1. `seams_version_major`, `seams_version_minor`, `seams_version_patch`
and `seams_version_string` identify the release. The `seams_version` prefix
is those four calls. The build writes those
macros from the project version (2.14.0).

`seams_status` is `SEAMS_OK`, `SEAMS_ERR_NULL`, `SEAMS_ERR_RANGE`, `SEAMS_ERR_BOX`,
`SEAMS_ERR_ARG`, `SEAMS_ERR_STATE`, `SEAMS_ERR_NOMEM`, or `SEAMS_ERR_INTERNAL`.
`seams_status_string` names a code. `seams_last_error` is the thread-local detail.
The enumerator types are `seams_status`, `seams_label`, `seams_ring_symmetry`,
`seams_chill_rule`, and `seams_graph`.

## Threads

`seams_set_num_threads` sets the OpenMP team used by the next call. Zero restores
the library default. A negative value returns `SEAMS_ERR_RANGE`.
`seams_get_num_threads` reports the width that will be used.

## Frames

`seams_system_create` and `seams_system_destroy` own one frame.
`seams_system_set_frame` copies `n` points and a nine-number cell. `types` and `ids`
may be null. `seams_system_set_ortho` is the same call with a zero tilt.
`seams_system_count` is the number of points, or zero when the handle is empty.

## Neighbours

`seams_neighbours_build` stores a cutoff graph (`SEAMS_GRAPH_CUTOFF`) or the mutual
four-nearest graph (`SEAMS_GRAPH_MUTUAL`). `type` less than or equal to zero keeps
every atom. `seams_neighbours_count` is the row count.
`seams_neighbours_csr` returns offsets of length `n + 1` and the neighbour columns.
The self entry is omitted. `seams_neighbours_destroy` frees the handle.

## CHILL+

`seams_chill_plus` writes one label per atom and, when `counts9` is not null, nine
population counts. `rule` is `SEAMS_CHILL_PUBLISHED` (Nguyen and Molinero) or
`SEAMS_CHILL_OVITO`. `seams_label_name` maps a label to cubic, hexagonal, water,
interfacial, clathrate, interClathrate, unclassified, reCubic, or reHex.

## q4 and q6

`seams_q4q6` writes the local Steinhardt parameters on the cutoff graph.
Either pointer may be null.

## Rings and symmetry

`seams_rings_find` builds primitive rings up to `max_depth` (3 through 12).
`seams_rings_count` is the number of rings. `seams_rings_histogram` fills
`counts[k]` with the number of rings of size `k`. `seams_rings_labels` writes
the size and a symmetry code for each ring. Six-rings that pass the cage tests
are `SEAMS_RING_HC_BASAL`, `SEAMS_RING_HC_PRISMATIC`, `SEAMS_RING_DDC`,
`SEAMS_RING_BOTH_BASAL`, or `SEAMS_RING_BOTH_PRISMATIC`. Other rings are
`SEAMS_RING_PENTAGON`, `SEAMS_RING_HEXAGON`, or `SEAMS_RING_OTHER`.
`seams_ring_symmetry_name` prints that code. `seams_rings_fraction` is the
fraction of atoms that sit on a ring; a size of zero counts every ring.
`seams_rings_destroy` frees the handle.

`seams_set_num_threads(0)` restores the OpenMP width captured at first use
(`omp_get_max_threads`), which is the library default. A negative width
returns `SEAMS_ERR_RANGE`.

## Runnable example

`docs/examples/minimal_c.c` reads
`input/traj/genice_sI.lammpstrj`, builds the cutoff neighbour list, and
computes q6, CHILL+, and primitive rings. `meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_c.c
:language: c
```

## Examples

`meson test seams_c_api` compiles and runs `tests/test_seams_api.c`.
Every function in `seams.h` is called there. `meson test docs_ice_example`
runs the eight-atom program in `docs/examples/ice.c`.

```{literalinclude} ../../../tests/test_seams_api.c
:language: c
:lines: 44-93
```

```{literalinclude} ../../../tests/test_seams_api.c
:language: c
:lines: 102-216
```

```{literalinclude} ../../../docs/examples/ice.c
:language: c
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### SEAMS_ABI

Purpose: ABI number compiled into the header.

Parameters: none

Units: dimensionless integer 1

Example: `docs/examples/minimal_c.c`

### SEAMS_CHILL_OVITO

Purpose: CHILL+ rule used by OVITO.

Parameters: none

Units: dimensionless code 1

Example: `docs/examples/minimal_c.c`

### SEAMS_CHILL_PUBLISHED

Purpose: CHILL+ rule from Nguyen and Molinero, 2015.

Parameters: none

Units: dimensionless code 0

Example: `docs/examples/minimal_c.c`

### SEAMS_ERR_ARG

Purpose: A caller argument is the wrong shape.

Parameters: none

Units: dimensionless status code 4

Example: `docs/examples/minimal_c.c`

### SEAMS_ERR_BOX

Purpose: The cell is missing or has no volume.

Parameters: none

Units: dimensionless status code 3

Example: `docs/examples/minimal_c.c`

### SEAMS_ERR_INTERNAL

Purpose: The library hit an internal failure.

Parameters: none

Units: dimensionless status code 7

Example: `docs/examples/minimal_c.c`

### SEAMS_ERR_NOMEM

Purpose: An allocation failed.

Parameters: none

Units: dimensionless status code 6

Example: `docs/examples/minimal_c.c`

### SEAMS_ERR_NULL

Purpose: A required pointer is null.

Parameters: none

Units: dimensionless status code 1

Example: `docs/examples/minimal_c.c`

### SEAMS_ERR_RANGE

Purpose: A number is outside the accepted range.

Parameters: none

Units: dimensionless status code 2

Example: `docs/examples/minimal_c.c`

### SEAMS_ERR_STATE

Purpose: The handle has no frame yet.

Parameters: none

Units: dimensionless status code 5

Example: `docs/examples/minimal_c.c`

### SEAMS_GRAPH_CUTOFF

Purpose: Neighbour graph of every pair inside the cutoff.

Parameters: none

Units: dimensionless code 0

Example: `docs/examples/minimal_c.c`

### SEAMS_GRAPH_MUTUAL

Purpose: Mutual four-nearest neighbour graph.

Parameters: none

Units: dimensionless code 1

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_CLATHRATE

Purpose: CHILL+ clathrate label.

Parameters: none

Units: dimensionless code 4

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_CUBIC

Purpose: CHILL+ cubic-ice label.

Parameters: none

Units: dimensionless code 0

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_HEXAGONAL

Purpose: CHILL+ hexagonal-ice label.

Parameters: none

Units: dimensionless code 1

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_INTERCLATHRATE

Purpose: CHILL+ interfacial-clathrate label.

Parameters: none

Units: dimensionless code 5

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_INTERFACIAL

Purpose: CHILL+ interfacial-ice label.

Parameters: none

Units: dimensionless code 3

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_RECUBIC

Purpose: CHILL+ reclassified cubic label.

Parameters: none

Units: dimensionless code 7

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_REHEX

Purpose: CHILL+ reclassified hexagonal label.

Parameters: none

Units: dimensionless code 8

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_UNCLASSIFIED

Purpose: CHILL+ label for an atom that matched no class.

Parameters: none

Units: dimensionless code 6

Example: `docs/examples/minimal_c.c`

### SEAMS_LABEL_WATER

Purpose: CHILL+ liquid-water label.

Parameters: none

Units: dimensionless code 2

Example: `docs/examples/minimal_c.c`

### SEAMS_OK

Purpose: The call succeeded.

Parameters: none

Units: dimensionless status code 0

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_BOTH_BASAL

Purpose: Six-ring that is a basal face of both cage types.

Parameters: none

Units: dimensionless code 4

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_BOTH_PRISMATIC

Purpose: Six-ring that is a prismatic face of both cage types.

Parameters: none

Units: dimensionless code 5

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_DDC

Purpose: Six-ring of a double-diamond cage.

Parameters: none

Units: dimensionless code 1

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_DEFORMED_PRISM

Purpose: Deformed prismatic six-ring.

Parameters: none

Units: dimensionless code 7

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_HC_BASAL

Purpose: Basal six-ring of a hexagonal cage.

Parameters: none

Units: dimensionless code 2

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_HC_PRISMATIC

Purpose: Prismatic six-ring of a hexagonal cage.

Parameters: none

Units: dimensionless code 3

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_HEXAGON

Purpose: Six-ring that is not a cage face.

Parameters: none

Units: dimensionless code 10

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_MIXED_PRISM

Purpose: Prismatic six-ring mixing the two cage types.

Parameters: none

Units: dimensionless code 8

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_OTHER

Purpose: Ring that is neither a pentagon nor a hexagon cage face.

Parameters: none

Units: dimensionless code 11

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_PENTAGON

Purpose: Five-membered ring.

Parameters: none

Units: dimensionless code 9

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_PRISM

Purpose: Prismatic six-ring.

Parameters: none

Units: dimensionless code 6

Example: `docs/examples/minimal_c.c`

### SEAMS_RING_UNCLASSIFIED

Purpose: Ring with no symmetry class yet.

Parameters: none

Units: dimensionless code 0

Example: `docs/examples/minimal_c.c`

### seams_abi

Purpose: Return the ABI number of this libseams.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### seams_chill_plus

Purpose: Write one CHILL+ label per atom.

Parameters: sys; cutoff; type; graph; rule; labels; counts9

Units: cutoff and coordinates in angstrom; labels and counts are dimensionless

Example: `docs/examples/minimal_c.c`

### seams_chill_rule

Purpose: Enumerator of the published and OVITO CHILL+ rules.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_c.c`

### seams_get_num_threads

Purpose: Report the OpenMP width the next call will use.

Parameters: none

Units: dimensionless thread count

Example: `docs/examples/minimal_c.c`

### seams_graph

Purpose: Enumerator of the cutoff graph and the mutual graph.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_c.c`

### seams_label

Purpose: Enumerator of the nine CHILL+ labels.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_c.c`

### seams_label_name

Purpose: Map a CHILL+ code to its name.

Parameters: label

Units: dimensionless code in, text out

Example: `docs/examples/minimal_c.c`

### seams_last_error

Purpose: Return the thread-local error detail.

Parameters: none

Units: text

Example: `docs/examples/minimal_c.c`

### seams_neighbours

Purpose: Opaque handle of one neighbour list.

Parameters: none

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_neighbours_build

Purpose: Build a cutoff or mutual neighbour list.

Parameters: sys; cutoff; type; graph; out

Units: cutoff in angstrom; type 0 keeps every atom

Example: `docs/examples/minimal_c.c`

### seams_neighbours_count

Purpose: Return the number of neighbour rows.

Parameters: nb

Units: dimensionless count

Example: `docs/examples/minimal_c.c`

### seams_neighbours_csr

Purpose: Return the neighbour list as CSR offsets and columns.

Parameters: nb; offsets; n_offsets; cols; nnz

Units: dimensionless indices; the self entry is omitted

Example: `docs/examples/minimal_c.c`

### seams_neighbours_destroy

Purpose: Free a neighbour-list handle.

Parameters: nb

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_q4q6

Purpose: Write local q4 and q6 on the cutoff shell.

Parameters: sys; cutoff; type; q4; q6

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_c.c`

### seams_ring_symmetry

Purpose: Enumerator of primitive-ring symmetry classes.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_c.c`

### seams_ring_symmetry_name

Purpose: Map a ring-symmetry code to its name.

Parameters: symmetry

Units: dimensionless code in, text out

Example: `docs/examples/minimal_c.c`

### seams_rings

Purpose: Opaque handle of one primitive-ring set.

Parameters: none

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_rings_count

Purpose: Return the number of primitive rings.

Parameters: rings

Units: dimensionless count

Example: `docs/examples/minimal_c.c`

### seams_rings_destroy

Purpose: Free a ring-set handle.

Parameters: rings

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_rings_find

Purpose: Find primitive rings up to a maximum size.

Parameters: sys; cutoff; type; graph; max_depth; out

Units: cutoff in angstrom; max_depth is a dimensionless ring size from 3

Example: `docs/examples/minimal_c.c`

### seams_rings_fraction

Purpose: Fraction of atoms that sit on a ring of one size.

Parameters: rings; size; fraction

Units: size 0 counts every ring; fraction is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_rings_histogram

Purpose: Count primitive rings by size.

Parameters: rings; counts; counts_len

Units: counts[k] is the dimensionless number of rings of size k

Example: `docs/examples/minimal_c.c`

### seams_rings_labels

Purpose: Write the size and symmetry of each ring.

Parameters: rings; sizes; symmetry

Units: dimensionless size and symmetry code

Example: `docs/examples/minimal_c.c`

### seams_set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Units: dimensionless thread count; 0 restores the library default

Example: `docs/examples/minimal_c.c`

### seams_status

Purpose: Enumerator of the status codes.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_c.c`

### seams_status_string

Purpose: Map a status code to its name.

Parameters: status

Units: dimensionless code in, text out

Example: `docs/examples/minimal_c.c`

### seams_system

Purpose: Opaque handle of one frame.

Parameters: none

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_system_count

Purpose: Return the number of points in the frame.

Parameters: sys

Units: dimensionless count

Example: `docs/examples/minimal_c.c`

### seams_system_create

Purpose: Allocate an empty frame handle.

Parameters: none

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_system_destroy

Purpose: Free a frame handle.

Parameters: sys

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_system_set_frame

Purpose: Store points and a general LAMMPS cell.

Parameters: sys; xyz; n; box9; types; ids

Units: coordinates and box lengths in angstrom; n is a dimensionless count

Example: `docs/examples/minimal_c.c`

### seams_system_set_ortho

Purpose: Store points and an orthorhombic cell.

Parameters: sys; xyz; n; ox; oy; oz; lx; ly; lz; types; ids

Units: origin and edge lengths in angstrom

Example: `docs/examples/minimal_c.c`

### seams_version

Purpose: Stem of the generated version header and the version calls.

Parameters: none

Units: text and dimensionless integers

Example: `docs/examples/minimal_c.c`

### seams_version_major

Purpose: Return the major version.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### seams_version_minor

Purpose: Return the minor version.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### seams_version_patch

Purpose: Return the patch version.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### seams_version_string

Purpose: Return the version as text.

Parameters: none

Units: text

Example: `docs/examples/minimal_c.c`
<!-- API_SYMBOLS_END -->
