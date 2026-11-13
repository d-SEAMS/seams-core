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
population counts. `rule` is `SEAMS_CHILL_PUBLISHED` (four closest neighbours in
q_lm, Nguyen and Molinero equation 1) or `SEAMS_CHILL_OVITO` (every neighbour
inside the cutoff). `seams_label_name` maps a label to cubic, hexagonal, water,
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

## Frames from a file

`seams_path_nframes` reports how many frames a path holds.
`seams_system_read` copies frame 1, or a later 1-based index, into a
handle. A negative type filter keeps every atom. `seams_system_set_atoms`
is the array entry with `size_t` lengths. `seams_system_set_frame` remains
the int-length entry.

## Order parameters, cages, symmetry, and the front

`seams_ql` writes local q3, q4, or q6, and the Lechner-Dellago average
when the average pointer is present. `seams_q6bar` summarises that
average. The nucleation cut used in the quickstart is 0.36.
`seams_chill` is the original CHILL rule. `SEAMS_CHILL_PUBLISHED` keeps
inclusive edges and `SEAMS_CHILL_OVITO` keeps the open edges.
`seams_chill_plus` is unchanged.
`seams_cages_find` scores hexagonal and double-diamond membership.
`SEAMS_GRAPH_SEEDED` is the mutual-seed graph. `seams_cages_signature`
counts closed cages of one census, such as `512`.
`seams_system_symmetry` writes the shell label of each atom. A
non-positive residual selects the library ceiling.
`seams_brine_analyse` measures a planar ice front from a caller-supplied
ice assignment.
`seams_steinhardt` is the cutoff-shell q_l, including q3, with the
Lechner-Dellago average when the second buffer is present.
`seams_neighbours_knn` builds a k-nearest list, mutual or one-way.
`seams_rdf`, `seams_running_cn`, and `seams_coordination` are the
partial distribution and its integral.
`seams_hbonds` counts hydrogen bonds. `seams_density_axis` is the
number density along one edge. `seams_fingerprint` names the bonded
graph. `seams_prisms` and `seams_monolayer` are the prism and sheet
censuses. `seams_site_table`, `seams_contact_pairs`, and `seams_domains`
read a site spec such as `1=waterO`. `seams_cages` and
`seams_signature_cages` are the cage searches.

`seams_handle_error` is the detail stored on the handle that failed.
`seams_last_error` remains the thread-local detail.
`seams_neighbours_copy` fills caller buffers. A null buffer queries the
length. `seams_neighbours_csr` still returns pointers owned by the handle
until `seams_neighbours_destroy`.

One handle is used by one thread. Distinct handles may run together.
`seams_set_num_threads` is process-global.

`examples/c` builds the quickstart against the installed library.
`pkg-config --cflags --libs seams` and `find_package(seams)` both name
that install. The design of the ABI is
`docs/orgmode/explanation/c-api.org`.

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
:lines: 102-466
```

```{literalinclude} ../../../examples/c/quickstart.c
:language: c
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

### SEAMS_CAGE_BOTH

Purpose: Atom that belongs to both cage types.

Parameters: none

Units: dimensionless code 3

Example: `docs/examples/minimal_c.c`

### SEAMS_CAGE_CUBIC

Purpose: Atom that belongs only to a double-diamond cage.

Parameters: none

Units: dimensionless code 2

Example: `docs/examples/minimal_c.c`

### SEAMS_CAGE_HEX

Purpose: Atom that belongs only to a hexagonal cage.

Parameters: none

Units: dimensionless code 1

Example: `docs/examples/minimal_c.c`

### SEAMS_CAGE_WATER

Purpose: Atom that belongs to neither cage.

Parameters: none

Units: dimensionless code 0

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

### SEAMS_DEPRECATED

Purpose: Mark a shipped name and name its replacement.

Parameters: message

Units: attribute text

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

### SEAMS_GRAPH_SEEDED

Purpose: Cage graph: mutual seeds and union completions.

Parameters: none

Units: dimensionless code 2

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

### SEAMS_SITE_ACCEPTOR

Purpose: Site kind of an acceptor atom.

Parameters: none

Units: dimensionless code 5

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_ANION

Purpose: Site kind of an anion.

Parameters: none

Units: dimensionless code 2

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_APOLAR

Purpose: Site kind of an apolar atom.

Parameters: none

Units: dimensionless code 7

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_CATION

Purpose: Site kind of a cation head.

Parameters: none

Units: dimensionless code 1

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_DONOR_H

Purpose: Site kind of a donor hydrogen.

Parameters: none

Units: dimensionless code 4

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_POLAR

Purpose: Site kind of a polar atom.

Parameters: none

Units: dimensionless code 6

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_SOLVENT

Purpose: Site kind of a solvent atom.

Parameters: none

Units: dimensionless code 10

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_TAIL

Purpose: Site kind of a tail atom.

Parameters: none

Units: dimensionless code 3

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_UNSPECIFIED

Purpose: Site kind that the table did not name.

Parameters: none

Units: dimensionless code 0

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_WATER_H

Purpose: Site kind of a water hydrogen.

Parameters: none

Units: dimensionless code 9

Example: `docs/examples/minimal_c.c`

### SEAMS_SITE_WATER_O

Purpose: Site kind of a water oxygen.

Parameters: none

Units: dimensionless code 8

Example: `docs/examples/minimal_c.c`

### SEAMS_VERSION_MAJOR

Purpose: Major component of this release.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### SEAMS_VERSION_MINOR

Purpose: Minor component of this release.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### SEAMS_VERSION_NUMBER

Purpose: Packed version, (major << 16) | (minor << 8) | patch.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### SEAMS_VERSION_PATCH

Purpose: Patch component of this release.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### SEAMS_VERSION_STRING

Purpose: Release text written by the build.

Parameters: none

Units: text

Example: `docs/examples/minimal_c.c`

### seams_abi

Purpose: Return the ABI number of this libseams.

Parameters: none

Units: dimensionless integer

Example: `docs/examples/minimal_c.c`

### seams_brine

Purpose: Opaque handle of one ice-front measurement.

Parameters: none

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_brine_analyse

Purpose: Measure a planar ice front from an ice assignment.

Parameters: sys; ice; n_ice; water_type; ion_type; n_bins; axis; out

Units: lengths in angstrom; ice flags are dimensionless

Example: `docs/examples/minimal_c.c`

### seams_brine_destroy

Purpose: Free a brine-front handle.

Parameters: front

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_brine_fronts

Purpose: Report the two planar crossings of the ice slab.

Parameters: front; z_enter; z_exit; thickness; found

Units: lengths in angstrom

Example: `docs/examples/minimal_c.c`

### seams_brine_pockets

Purpose: Count brine pockets and the liquid reservoir.

Parameters: front; n_pockets; reservoir_waters

Units: dimensionless counts

Example: `docs/examples/minimal_c.c`

### seams_brine_profile

Purpose: Copy the concentration profile along the front axis.

Parameters: front; z_mid; molality; ice_fraction; cap; n_out

Units: z in angstrom; molality in mol/kg; ice fraction is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_brine_rejection

Purpose: Report the molal rejection coefficient.

Parameters: front; c_ice; c_liquid; k; defined

Units: molality in mol/kg; k is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_cage_label

Purpose: Enumerator of per-atom cage membership.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_c.c`

### seams_cage_set

Purpose: Opaque handle of hexagonal and double-diamond membership.

Parameters: none

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_cages

Purpose: Build hexagonal and double-diamond cage membership.

Parameters: sys; cutoff; type; graph; k; complete; out

Units: cutoff in angstrom

Example: `docs/examples/minimal_c.c`

### seams_cages_counts

Purpose: Disjoint cage population counts.

Parameters: cages; n_hex; n_cubic; n_water; n_both

Units: dimensionless counts

Example: `docs/examples/minimal_c.c`

### seams_cages_destroy

Purpose: Free a cage-membership handle.

Parameters: cages

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_cages_find

Purpose: Score hexagonal and double-diamond cages.

Parameters: sys; cutoff; type; graph; k; complete; out

Units: cutoff in angstrom; k is a dimensionless neighbour count

Example: `docs/examples/minimal_c.c`

### seams_cages_labels

Purpose: Write one cage code per atom.

Parameters: cages; per_atom; len

Units: dimensionless codes

Example: `docs/examples/minimal_c.c`

### seams_cages_signature

Purpose: Count closed cages of one ring census.

Parameters: sys; cutoff; type; graph; spec; n_cages; n_atoms

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_c.c`

### seams_chill

Purpose: Write one original CHILL label per atom.

Parameters: sys; cutoff; type; graph; rule; labels; labels_len; counts9

Units: cutoff in angstrom; labels are dimensionless

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

### seams_contact_pairs

Purpose: Count mutual nearest unlike contacts of a site table.

Parameters: sys; site_spec; n_pairs; n_cation; n_anion

Units: dimensionless counts

Example: `docs/examples/minimal_c.c`

### seams_coordination

Purpose: Coordination number integrated to rmax.

Parameters: sys; rmax; bins; type_i; type_j; cn

Units: rmax in angstrom; cn is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_density_axis

Purpose: Number density along one cell axis.

Parameters: sys; type; bins; axis; z; rho; cap; n_out

Units: z in angstrom; rho is number density

Example: `docs/examples/minimal_c.c`

### seams_domains

Purpose: Largest polar or apolar domain.

Parameters: sys; cutoff; site_spec; subset; n_subset; largest; percolation

Units: cutoff in angstrom; percolation is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_fingerprint

Purpose: Topological fingerprint of the bonded graph.

Parameters: sys; cutoff; type; k; graph; hops; colour_types; out

Units: cutoff in angstrom; hops is a dimensionless depth

Example: `docs/examples/minimal_c.c`

### seams_fingerprint_atom

Purpose: Local key of one atom.

Parameters: fp; index; buf; cap; n_out

Units: text

Example: `docs/examples/minimal_c.c`

### seams_fingerprint_destroy

Purpose: Free a fingerprint handle.

Parameters: fp

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_fingerprint_key

Purpose: Frame key of a fingerprint.

Parameters: fp; buf; cap; n_out

Units: text

Example: `docs/examples/minimal_c.c`

### seams_fingerprint_rings

Purpose: Primitive-ring census stored on a fingerprint.

Parameters: fp; counts; counts_len

Units: dimensionless ring counts

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

### seams_handle_error

Purpose: Return the detail string stored on a handle.

Parameters: handle

Units: text

Example: `docs/examples/minimal_c.c`

### seams_hbonds

Purpose: Count hydrogen bonds on one frame.

Parameters: sys; cutoff; oxygen_type; hydrogen_type; dist_cutoff; angle_cutoff; n_bonds

Units: distances in angstrom; angle in degrees; n_bonds is dimensionless

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

### seams_monolayer

Purpose: Primitive-ring histogram of a monolayer.

Parameters: sys; cutoff; type; max_depth; sheet_area; counts; counts_len

Units: cutoff in angstrom; sheet_area in square angstrom

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

### seams_neighbours_copy

Purpose: Copy the neighbour list into caller buffers.

Parameters: nb; offsets; offsets_cap; n_offsets; cols; cols_cap; nnz

Units: dimensionless indices; a null buffer queries the length

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

### seams_neighbours_knn

Purpose: Build a k-nearest neighbour list.

Parameters: sys; k; cutoff; type; mutual; out

Units: cutoff in angstrom; k is a dimensionless neighbour count

Example: `docs/examples/minimal_c.c`

### seams_path_nframes

Purpose: Count frames in a trajectory path.

Parameters: path; nframes

Units: dimensionless count

Example: `docs/examples/minimal_c.c`

### seams_prisms

Purpose: Perfect and imperfect prism counts by ring size.

Parameters: sys; cutoff; type; perfect; imperfect; len

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_c.c`

### seams_q4q6

Purpose: Write local q4 and q6 on the cutoff shell.

Parameters: sys; cutoff; type; q4; q6

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_c.c`

### seams_q6bar

Purpose: Write Lechner-Dellago q6-bar and the ice count above a cut.

Parameters: sys; cutoff; type; ice_cut; q6bar; q6bar_len; mean; median; n_ice

Units: cutoff in angstrom; q6-bar is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_ql

Purpose: Write local q_l and, on request, its neighbour average.

Parameters: sys; cutoff; type; order_l; graph; ql; ql_len; qlbar; qlbar_len

Units: cutoff in angstrom; q_l is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_rdf

Purpose: Partial radial distribution.

Parameters: sys; rmax; bins; type_i; type_j; r; g; count; cap; n_out; volume

Units: r and rmax in angstrom; g is dimensionless; volume in cubic angstrom

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

### seams_rings_vertices

Purpose: Copy the atom indices of one ring.

Parameters: rings; index; atoms; cap; n_out

Units: dimensionless indices

Example: `docs/examples/minimal_c.c`

### seams_running_cn

Purpose: Running coordination number along the partial RDF.

Parameters: sys; rmax; bins; type_i; type_j; r; cn; cap; n_out

Units: r in angstrom; cn is dimensionless

Example: `docs/examples/minimal_c.c`

### seams_set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Units: dimensionless thread count; 0 restores the library default

Example: `docs/examples/minimal_c.c`

### seams_signature_cages

Purpose: Count closed cages of one ring census.

Parameters: sys; cutoff; type; graph; spec; n_cages; n_atoms

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_c.c`

### seams_site_kind

Purpose: Enumerator of site-table kinds.

Parameters: none

Units: dimensionless

Example: `docs/examples/minimal_c.c`

### seams_site_table

Purpose: Per-atom site kind from a site-table spec.

Parameters: sys; spec; kinds; len

Units: dimensionless kind codes

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

### seams_steinhardt

Purpose: Local Steinhardt q_l on the cutoff shell, with the average on request.

Parameters: sys; cutoff; type; order_l; ql; ql_len; qlbar; qlbar_len

Units: cutoff in angstrom; q_l is dimensionless

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

### seams_system_read

Purpose: Read one frame from a trajectory path.

Parameters: sys; path; frame; type_filter

Units: frame index is 1-based; coordinates in angstrom

Example: `docs/examples/minimal_c.c`

### seams_system_set_atoms

Purpose: Store points and a cell, with size_t lengths.

Parameters: sys; xyz; n; box9; box_len; types; ids

Units: coordinates and box lengths in angstrom

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

### seams_system_symmetry

Purpose: Write one ice-shell label per atom.

Parameters: sys; shell_rmsd; labels; labels_len

Units: shell_rmsd is a dimensionless residual; labels are dimensionless codes

Example: `docs/examples/minimal_c.c`

### seams_topo

Purpose: Opaque handle of one topological fingerprint.

Parameters: none

Units: handle

Example: `docs/examples/minimal_c.c`

### seams_version

Purpose: Return the packed release number.

Parameters: none

Units: dimensionless integer, (major << 16) | (minor << 8) | patch

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
