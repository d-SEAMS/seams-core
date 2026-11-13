# Fortran

`bindings/fortran/seams.f90` is an `iso_c_binding` module. Every procedure is
`bind(C)` to the symbol in `seams.h`.

`seams_set_num_threads` and `seams_get_num_threads` are in the same module.
`bindings/fortran/test_seams.f90` is the program the language-front script
compiles. It calls every `bind(C)` procedure in the module on a 64-atom
cubic cell. `meson test docs_api_fortran` compiles that program with
gfortran and runs it.
`docs/examples/minimal_fortran.f90` reads
`input/traj/genice_sI.lammpstrj` and calls the same C entry points.
`meson test docs_examples` runs it.

The module exports `SEAMS_OK`, `SEAMS_CHILL_PUBLISHED`, `SEAMS_CHILL_OVITO`,
`SEAMS_GRAPH_CUTOFF`, `SEAMS_GRAPH_MUTUAL`, `SEAMS_LABEL_CUBIC`,
`seams_version_string`, `seams_last_error`, `seams_set_num_threads`,
`seams_get_num_threads`, `seams_system_create`, `seams_system_destroy`,
`seams_system_set_ortho`, `seams_system_set_frame`, `seams_chill_plus`,
`seams_q4q6`, `seams_neighbours_build`, `seams_neighbours_destroy`,
`seams_neighbours_csr`, `seams_rings_find`, `seams_rings_destroy`,
`seams_rings_histogram`, and `seams_rings_fraction`.

```{literalinclude} ../../../docs/examples/minimal_fortran.f90
:language: fortran
```

```{literalinclude} ../../../bindings/fortran/test_seams.f90
:language: fortran
:lines: 37-74
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### SEAMS_CHILL_OVITO

Purpose: CHILL+ rule used by OVITO.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CHILL_PUBLISHED

Purpose: CHILL+ rule from the 2015 paper.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_GRAPH_CUTOFF

Purpose: Neighbour graph of pairs inside the cutoff.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_GRAPH_MUTUAL

Purpose: Mutual four-nearest neighbour graph.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_CUBIC

Purpose: CHILL+ cubic-ice label. Fortran arrays of labels are 1-based, so clathrate is counts(5).

Parameters: none

Units: dimensionless code 0

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_OK

Purpose: The call succeeded.

Parameters: none

Units: dimensionless status code

Example: `docs/examples/minimal_fortran.f90`

### seams_chill_plus

Purpose: Write one CHILL+ label per atom.

Parameters: sys; cutoff; type; graph; rule; labels; counts

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_get_num_threads

Purpose: Report the OpenMP width the next call will use.

Parameters: none

Units: dimensionless thread count

Example: `docs/examples/minimal_fortran.f90`

### seams_last_error

Purpose: Return the thread-local error detail.

Parameters: none

Units: text

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_build

Purpose: Build a cutoff or mutual neighbour list.

Parameters: sys; cutoff; type; graph; out

Units: cutoff in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_csr

Purpose: Return the neighbour list as CSR.

Parameters: nb; offsets; n_offsets; cols; nnz

Units: dimensionless indices

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_destroy

Purpose: Free a neighbour-list handle.

Parameters: nb

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_q4q6

Purpose: Write local q4 and q6.

Parameters: sys; cutoff; type; q4; q6

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_destroy

Purpose: Free a ring-set handle.

Parameters: rings

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_find

Purpose: Find primitive rings up to a maximum size.

Parameters: sys; cutoff; type; graph; max_depth; out

Units: cutoff in angstrom; max_depth is a dimensionless ring size

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_fraction

Purpose: Fraction of atoms on a ring of one size.

Parameters: rings; size; fraction

Units: dimensionless fraction

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_histogram

Purpose: Count primitive rings by size.

Parameters: rings; counts; counts_len

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Units: dimensionless thread count

Example: `docs/examples/minimal_fortran.f90`

### seams_system_create

Purpose: Allocate an empty frame handle.

Parameters: none

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_system_destroy

Purpose: Free a frame handle.

Parameters: sys

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_system_set_frame

Purpose: Store points and a nine-number cell.

Parameters: sys; xyz; n; box; types; ids

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_system_set_ortho

Purpose: Store points and an orthorhombic cell.

Parameters: sys; xyz; n; ox; oy; oz; lx; ly; lz; types; ids

Units: origin and edge lengths in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_version_string

Purpose: Return the version as text.

Parameters: none

Units: text

Example: `docs/examples/minimal_fortran.f90`
<!-- API_SYMBOLS_END -->
