# Fortran

`bindings/fortran/seams.f90` is an `iso_c_binding` module. Every procedure is
`bind(C)` to the symbol in `seams.h`.

`seams_set_num_threads` and `seams_get_num_threads` are in the same module.
`bindings/fortran/test_seams.f90` is the program the language-front script
compiles. It calls every `bind(C)` procedure in the module on a 64-atom
cubic cell. `meson test docs_api_fortran` compiles that program with
gfortran and runs it. The run prints `fortran call cubic 64`,
`fortran call q6 0.628539361055`, `fortran call cages 0 64 0 0`,
`fortran call sig 0 0`, `fortran call cn 4.063492063492`,
`fortran call hbond 0`, `fortran call sheet 128`,
`fortran call vert 6 0`, and `fortran call read 128 46`.
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
The same module also binds `seams_handle_error`, `seams_system_count`,
`seams_system_read`, `seams_chill`, `seams_steinhardt`, `seams_cages`,
`seams_cages_destroy`, `seams_cages_counts`, `seams_signature_cages`,
`seams_neighbours_knn`, `seams_rdf`, `seams_running_cn`,
`seams_coordination`, `seams_hbonds`, `seams_density_axis`,
`seams_fingerprint`, `seams_fingerprint_destroy`, `seams_fingerprint_key`,
`seams_fingerprint_rings`, `seams_contact_pairs`, `seams_domains`,
`seams_prisms`, `seams_monolayer`, `seams_rings_count`,
`seams_rings_vertices`, and `seams_site_table`. Fortran wrappers
`read_frame`, `signature_count`, `site_kinds`, `pair_count`, and
`domain_size` turn Fortran strings into the C buffers those calls take.

```{literalinclude} ../../../docs/examples/minimal_fortran.f90
:language: fortran
```

```{literalinclude} ../../../bindings/fortran/test_seams.f90
:language: fortran
:lines: 45-141
```

```{literalinclude} ../../../bindings/fortran/seams.f90
:language: fortran
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

### SEAMS_GRAPH_SEEDED

Purpose: Cage graph of mutual seeds and union completions.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_fortran.f90`

### c_chars

Purpose: Copy a Fortran string into a null-terminated C buffer.

Parameters: text

Units: text

Example: `docs/examples/minimal_fortran.f90`

### domain_size

Purpose: Largest apolar domain for a site-table spec.

Parameters: sys; spec; n_subset; largest; percolation

Units: percolation is dimensionless; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### pair_count

Purpose: Mutual nearest unlike contacts for a site-table spec.

Parameters: sys; spec; n_pairs; n_cation; n_anion

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### read_frame

Purpose: Read one trajectory frame into a handle.

Parameters: sys; path; frame; type_filter

Units: frame index is 1-based; coordinates in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_cages

Purpose: Hexagonal and double-diamond cage membership.

Parameters: sys; cutoff; typ; graph; k; complete; out

Units: cutoff in angstrom; handle

Example: `docs/examples/minimal_fortran.f90`

### seams_cages_counts

Purpose: Disjoint hexagonal, cubic, water, and both counts.

Parameters: cages; n_hex; n_cubic; n_water; n_both

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_cages_destroy

Purpose: Free a cage-membership handle.

Parameters: cages

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_chill

Purpose: Original CHILL labels and nine population counts.

Parameters: sys; cutoff; typ; graph; rule; labels; labels_len; counts

Units: cutoff in angstrom; codes are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_contact_pairs

Purpose: Mutual nearest unlike contacts of a site table.

Parameters: sys; spec; n_pairs; n_cation; n_anion

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_coordination

Purpose: Coordination number integrated to rmax.

Parameters: sys; rmax; bins; type_i; type_j; cn

Units: rmax in angstrom; coordination is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_density_axis

Purpose: Number density along one cell axis.

Parameters: sys; typ; bins; axis; z; rho; cap; n_out

Units: density in inverse cubic angstrom; centres in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_domains

Purpose: Largest polar or apolar domain.

Parameters: sys; cutoff; spec; subset; n_subset; largest; percolation

Units: cutoff in angstrom; percolation is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint

Purpose: Topological fingerprint of the bonded graph.

Parameters: sys; cutoff; typ; k; graph; hops; colour_types; out

Units: cutoff in angstrom; handle

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint_destroy

Purpose: Free a fingerprint handle.

Parameters: fp

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint_key

Purpose: Frame fingerprint text.

Parameters: fp; buf; cap; n_out

Units: text

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint_rings

Purpose: Primitive-ring census stored on the fingerprint.

Parameters: fp; counts; counts_len

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_handle_error

Purpose: Detail string stored on a handle.

Parameters: handle

Units: text

Example: `docs/examples/minimal_fortran.f90`

### seams_hbonds

Purpose: Hydrogen-bond count.

Parameters: sys; cutoff; oxygen_type; hydrogen_type; dist_cutoff; angle_cutoff; n_bonds

Units: distances in angstrom; angle in degrees

Example: `docs/examples/minimal_fortran.f90`

### seams_monolayer

Purpose: Primitive-ring histogram of a monolayer.

Parameters: sys; cutoff; typ; max_depth; sheet_area; counts; counts_len

Units: cutoff in angstrom; area in square angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_knn

Purpose: k-nearest neighbour list.

Parameters: sys; k; cutoff; typ; mutual; out

Units: cutoff in angstrom; k is a dimensionless count

Example: `docs/examples/minimal_fortran.f90`

### seams_prisms

Purpose: Perfect and imperfect prism counts by ring size.

Parameters: sys; cutoff; typ; perfect; imperfect; len

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_rdf

Purpose: Partial radial distribution and the cell volume.

Parameters: sys; rmax; bins; type_i; type_j; r; g; count; cap; n_out; volume

Units: rmax in angstrom; volume in cubic angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_count

Purpose: Number of primitive rings.

Parameters: rings

Units: dimensionless count

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_vertices

Purpose: Atom indices of one ring.

Parameters: rings; index; atoms; cap; n_out

Units: dimensionless indices

Example: `docs/examples/minimal_fortran.f90`

### seams_running_cn

Purpose: Running coordination number along the partial RDF.

Parameters: sys; rmax; bins; type_i; type_j; r; cn; cap; n_out

Units: rmax in angstrom; coordination is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_signature_cages

Purpose: Closed cages whose ring census matches a signature.

Parameters: sys; cutoff; typ; graph; spec; n_cages; n_atoms

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_site_table

Purpose: Per-atom site kind from a site-table spec.

Parameters: sys; spec; kinds; len

Units: dimensionless site codes

Example: `docs/examples/minimal_fortran.f90`

### seams_steinhardt

Purpose: Local q_l and the Lechner-Dellago average.

Parameters: sys; cutoff; typ; order_l; ql; ql_len; qlbar; qlbar_len

Units: cutoff in angstrom; q_l is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_system_count

Purpose: Number of points in the frame.

Parameters: sys

Units: dimensionless count

Example: `docs/examples/minimal_fortran.f90`

### seams_system_read

Purpose: Read one frame from a trajectory path.

Parameters: sys; path; frame; type_filter

Units: frame index is 1-based; coordinates in angstrom

Example: `docs/examples/minimal_fortran.f90`

### signature_count

Purpose: Closed-cage count for one signature on the cutoff graph.

Parameters: sys; spec; n_cages; n_atoms

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### site_kinds

Purpose: Per-atom site kind for a Fortran spec string.

Parameters: sys; spec; kinds; n

Units: dimensionless site codes

Example: `docs/examples/minimal_fortran.f90`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/fortran/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug fortran` writes the 2400 by 1800 PNG.
```

