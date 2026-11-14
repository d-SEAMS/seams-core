# C

`seams.h` is ABI 1 of libseams. A handle is opaque. Coordinates are Angstrom.
A cell is nine numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`. A zero tilt is orthorhombic.
`lx, ly, lz` are edge lengths.

`docs/examples/c_api_tour.c` calls every exported function on the in-tree fixtures.
`meson test docs_c_api_tour` compiles that file, runs it, and requires each printed
`capi` line below to match the run.

## Handle lifecycle

`seams_system_create` allocates an empty frame. The caller owns it and passes it to
`seams_system_destroy`. `NULL` is a valid argument to every destroy function.
`seams_system_set_atoms` stores points and a nine-number cell. Lengths are `size_t`.
`seams_system_set_frame` is the same store with an `int` count.
`seams_system_set_ortho` stores an orthorhombic cell. `seams_system_read` copies one
1-based frame from a trajectory. A negative type filter keeps every atom. A positive
filter keeps that species. `seams_system_count` is the number of points, or 0 for a
null handle. Neighbour, ring, cage, brine, and fingerprint handles are separate
objects. Destroying the frame does not free them.

One thread uses a handle. Another thread may use a different handle at the same time.
`seams_ql`, `seams_q4q6`, `seams_q6bar`, and `seams_steinhardt` share one
process-global scratch, so those four calls do not overlap.

## Status codes

`seams_status` is zero on success.

| code | value | `seams_status_string` |
| --- | --- | --- |
| `SEAMS_OK` | 0 | ok |
| `SEAMS_ERR_NULL` | 1 | null pointer |
| `SEAMS_ERR_RANGE` | 2 | value out of range |
| `SEAMS_ERR_BOX` | 3 | box is not a positive cell |
| `SEAMS_ERR_ARG` | 4 | bad argument |
| `SEAMS_ERR_STATE` | 5 | system has no frame |
| `SEAMS_ERR_NOMEM` | 6 | out of memory |
| `SEAMS_ERR_INTERNAL` | 7 | internal error |

## seams_last_error and seams_handle_error

`seams_last_error` is the detail of the last failing call on this thread. A later
successful call on the same thread clears it. The storage is thread-local.
`seams_handle_error` is the copy stored on the handle that failed. That copy stays
until the next call on that handle, or until destroy. Pass `NULL` and
`seams_handle_error` returns the thread-local string. The caller does not free either
pointer.

On an empty frame, `seams_chill_plus` returns `SEAMS_ERR_STATE`. The tour prints:

```
capi seams_last_error system has no frame
capi seams_handle_error no-frame
```

## Threads

`seams_set_num_threads(1)` selects one OpenMP thread. `seams_set_num_threads(0)`
restores the library default. A negative width returns `SEAMS_ERR_RANGE`.
`seams_get_num_threads` reports the width the next call will use.
`seams_set_num_threads` is process-global. Do not call it while another thread is
inside `seams_get_num_threads`.

## What the tour runs

The eight-atom cell is the diamond in `docs/examples/ice.c`. Every atom is cubic.
`input/traj/ice1h.lammpstrj` is the 128-oxygen ice Ih cell. `input/traj/genice_sI.lammpstrj`
is the 46-oxygen sI cage. `input/traj/exampleTraj.lammpstrj` is the mixed TIP4P dump
used for hydrogen bonds, sites, and domains. The cutoff is 3.5 angstroms unless a
line below says otherwise.

`docs/examples/core_counts.c` reads `input/traj/ice1h.lammpstrj` and
`input/traj/genice_sII.lammpstrj` and checks the C core: 128 hexagonal
waters and 256 six-rings, then 136 clathrate waters, 144 five-rings and
16 six-rings. `meson test docs_core_counts` compiles and runs it.

```{literalinclude} ../../../docs/examples/core_counts.c
:language: c
```

## Functions


### seams_abi

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_abi 1
```

### seams_version

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_version 134656
```

### seams_version_major

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_version_major 2
```

### seams_version_minor

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_version_minor 14
```

### seams_version_patch

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_version_patch 0
```

### seams_version_string

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_version_string 2.14.0
```

### seams_status_string

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_status_string ok
```

### seams_label_name

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_label_name hexagonal
```

### seams_ring_symmetry_name

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_ring_symmetry_name hexagon
```

### seams_set_num_threads

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_set_num_threads 1
```

### seams_get_num_threads

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_get_num_threads 1
```

### seams_system_count

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_count 0
```

### seams_system_create

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_create 1
```

### seams_last_error

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_last_error system has no frame
```

### seams_handle_error

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_handle_error no-frame
```

### seams_system_destroy

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_destroy ok
```

### seams_neighbours_destroy

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_neighbours_destroy ok
```

### seams_rings_destroy

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_destroy ok
```

### seams_cages_destroy

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_cages_destroy ok
```

### seams_brine_destroy

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_brine_destroy ok
```

### seams_fingerprint_destroy

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_fingerprint_destroy ok
```

### seams_system_set_atoms

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_set_atoms 8
```

### seams_system_set_frame

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_set_frame 8
```

### seams_system_set_ortho

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_set_ortho 8
```

### seams_chill_plus

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_chill_plus cubic 8
```

### seams_path_nframes

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_path_nframes 1
```

### seams_system_read

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_read 128
```

### seams_neighbours_build

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_neighbours_build 128
```

### seams_neighbours_count

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_neighbours_count 128
```

### seams_neighbours_copy

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_neighbours_copy 129 512
```

### seams_neighbours_csr

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_neighbours_csr 129 512
```

### seams_neighbours_knn

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_neighbours_knn 128
```

### seams_chill

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_chill hexagonal 128
```

```
capi seams_chill_plus hexagonal 128
```

### seams_q4q6

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_q4q6 0.625855
```

### seams_ql

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_ql 0.562449
```

### seams_q6bar

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_q6bar 0.562449 128
```

### seams_steinhardt

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_steinhardt 0.562449
```

### seams_rings_find

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_find 256
```

### seams_rings_count

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_count 256
```

### seams_rings_histogram

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_histogram 256
```

### seams_rings_labels

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_labels 6
```

### seams_rings_vertices

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_vertices 6
```

### seams_rings_fraction

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_fraction 1.000000
```

### seams_cages_find

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_cages_find 128 0 0 0
```

### seams_cages_counts

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_cages_counts 128
```

### seams_cages_labels

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_cages_labels 0
```

### seams_cages

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_cages 0
```

```
capi seams_cages_destroy ok
```

### seams_cages_signature

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_cages_signature 2 40
```

### seams_signature_cages

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_signature_cages 2 40
```

### seams_system_symmetry

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_symmetry 128
```

### seams_brine_analyse

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_brine_analyse ok
```

### seams_brine_fronts

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_brine_fronts 1
```

### seams_brine_rejection

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_brine_rejection 0
```

### seams_brine_profile

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_brine_profile 20
```

### seams_brine_pockets

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_brine_pockets 1 48
```

### seams_rdf

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rdf 20 4158.665
```

### seams_running_cn

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_running_cn 26.205
```

### seams_coordination

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_coordination 4.031
```

### seams_hbonds

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_hbonds 460
```

### seams_density_axis

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_density_axis 4
```

### seams_fingerprint

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_fingerprint ok
```

### seams_fingerprint_key

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_fingerprint_key 16
```

### seams_fingerprint_atom

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_fingerprint_atom 16
```

### seams_fingerprint_rings

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_fingerprint_rings 256
```

### seams_contact_pairs

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_contact_pairs 0 250
```

### seams_domains

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_domains 750 750
```

### seams_prisms

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_prisms 0 0
```

### seams_monolayer

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_monolayer 256
```

### seams_site_table

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_site_table 250
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### SEAMS_ABI

Purpose: ABI number compiled into this header. Declared in include/seams.h.

Parameters: none

Returns: The macro value (1).

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_BOTH

Purpose: Both cages. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 3.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_CUBIC

Purpose: Double-diamond cage only. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 2.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_HEX

Purpose: Hexagonal cage only. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_WATER

Purpose: In neither cage. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 0.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_CHILL_OVITO

Purpose: Published edges are inclusive. OVITO edges are open. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_CHILL_PUBLISHED

Purpose: Published edges are inclusive. OVITO edges are open. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 0.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_DEPRECATED

Purpose: Mark a shipped declaration and name its replacement. Declared in include/seams.h.

Parameters: message

Returns: The macro value (attribute).

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_ARG

Purpose: A caller argument is the wrong shape. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 4.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_BOX

Purpose: The cell is missing or has no volume. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 3.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_INTERNAL

Purpose: The library hit an internal failure. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 7.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_NOMEM

Purpose: An allocation failed. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 6.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_NULL

Purpose: A required pointer is null. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_RANGE

Purpose: A number is outside the accepted range. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 2.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_STATE

Purpose: The handle has no frame yet. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 5.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_GRAPH_CUTOFF

Purpose: Neighbour graph. Seeded is a cage graph: mutual seeds, union completions. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 0.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_GRAPH_MUTUAL

Purpose: Neighbour graph. Seeded is a cage graph: mutual seeds, union completions. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_GRAPH_SEEDED

Purpose: Neighbour graph. Seeded is a cage graph: mutual seeds, union completions. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 2.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_CLATHRATE

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 4.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_CUBIC

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 0.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_HEXAGONAL

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_INTERCLATHRATE

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 5.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_INTERFACIAL

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 3.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_RECUBIC

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 7.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_REHEX

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 8.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_UNCLASSIFIED

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 6.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_WATER

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 2.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_OK

Purpose: The call succeeded. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 0.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_BOTH_BASAL

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 4.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_BOTH_PRISMATIC

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 5.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_DDC

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_DEFORMED_PRISM

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 7.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_HC_BASAL

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 2.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_HC_PRISMATIC

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 3.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_HEXAGON

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 10.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_MIXED_PRISM

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 8.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_OTHER

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 11.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_PENTAGON

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 9.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_PRISM

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 6.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_UNCLASSIFIED

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 0.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_ACCEPTOR

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 5.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_ANION

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 2.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_APOLAR

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 7.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_CATION

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_DONOR_H

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 4.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_POLAR

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 6.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_SOLVENT

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 10.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_TAIL

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 3.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_UNSPECIFIED

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 0.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_WATER_H

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 9.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_WATER_O

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumerator value 8.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_MAJOR

Purpose: Major component of this release. Declared in include/seams.h.

Parameters: none

Returns: The macro value (integer).

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_MINOR

Purpose: Minor component of this release. Declared in include/seams.h.

Parameters: none

Returns: The macro value (integer).

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_NUMBER

Purpose: Packed version, (major << 16) | (minor << 8) | patch. Declared in include/seams.h.

Parameters: none

Returns: The macro value (integer).

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_PATCH

Purpose: Patch component of this release. Declared in include/seams.h.

Parameters: none

Returns: The macro value (integer).

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_STRING

Purpose: Release text written by the build. Declared in include/seams.h.

Parameters: none

Returns: The macro value (text).

Ownership: No allocation.

Errors: None.

Units: text

Example: `docs/examples/c_api_tour.c`

### seams_abi

Purpose: ABI number compiled into this library. Declared in include/seams.h.

Parameters: none

Returns: SEAMS_ABI, which is 1.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_brine

Purpose: One planar ice/brine front measurement. Declared in include/seams.h.

Parameters: none

Returns: The opaque type.

Ownership: The caller creates the object and frees it with the matching destroy function.

Errors: None. This is a type.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_brine_analyse

Purpose: Planar ice/brine front from an ice assignment. Declared in include/seams.h.

Parameters: sys; ice; n_ice; water_type; ion_type; n_bins; axis; out

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller frees *out with seams_brine_destroy().

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_brine_destroy

Purpose: Free a brine-front handle. Declared in include/seams.h.

Parameters: front

Returns: Nothing.

Ownership: Frees the caller-owned handle. delete of a null pointer is a no-op.

Errors: None.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_brine_fronts

Purpose: The two planar crossings of the ice slab. Declared in include/seams.h.

Parameters: front; z_enter; z_exit; thickness; found

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_brine_pockets

Purpose: Brine pockets trapped in the ice, excluding the liquid reservoir. Declared in include/seams.h.

Parameters: front; n_pockets; reservoir_waters

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_brine_profile

Purpose: Concentration profile along the front axis. Declared in include/seams.h.

Parameters: front; z_mid; molality; ice_fraction; cap; n_out

Returns: SEAMS_OK when every non-NULL column was filled, or when every column pointer was NULL and n_out was written.

Ownership: The caller owns the columns.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_brine_rejection

Purpose: Molal rejection coefficient k = c_ice / c_liquid. Declared in include/seams.h.

Parameters: front; c_ice; c_liquid; k; defined

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_cage_label

Purpose: Per-atom cage membership from the hexagonal and double-diamond search. The four codes are disjoint. Declared in include/seams.h.

Parameters: none

Returns: The enumeration type.

Ownership: No allocation.

Errors: None. This is a type.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_cage_set

Purpose: Hexagonal and double-diamond membership for one frame. Declared in include/seams.h.

Parameters: none

Returns: The opaque type.

Ownership: The caller creates the object and frees it with the matching destroy function.

Errors: None. This is a type.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_cages

Purpose: Hexagonal and double-diamond cages. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; k; complete; out

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: On success the caller owns the new seams_cage_set and frees it with the matching destroy function. Buffers the caller passed stay with the caller.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_cages_counts

Purpose: Disjoint cage population counts. Declared in include/seams.h.

Parameters: cages; n_hex; n_cubic; n_water; n_both

Returns: SEAMS_OK, or SEAMS_ERR_NULL when cages is NULL.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_cages_destroy

Purpose: Free a cage-membership handle. Declared in include/seams.h.

Parameters: cages

Returns: Nothing.

Ownership: Frees the caller-owned handle. delete of a null pointer is a no-op.

Errors: None.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_cages_find

Purpose: Hexagonal and double-diamond cage membership. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; k; complete; out

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller frees *out with seams_cages_destroy().

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_cages_labels

Purpose: Per-atom cage code. Declared in include/seams.h.

Parameters: cages; per_atom; len

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: The caller owns per_atom.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_cages_signature

Purpose: Count closed cages whose ring census matches a signature. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; spec; n_cages; n_atoms

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_chill

Purpose: Original CHILL labels, published or OVITO edges. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; rule; labels; labels_len; counts9

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns labels and counts9.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_chill_plus

Purpose: CHILL+ labels. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; rule; labels; counts9

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns labels and counts9.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_chill_rule

Purpose: Published edges are inclusive. OVITO edges are open. Declared in include/seams.h.

Parameters: none

Returns: The enumeration type.

Ownership: No allocation.

Errors: None. This is a type.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_contact_pairs

Purpose: Mutual nearest unlike contacts of a site table. Declared in include/seams.h.

Parameters: sys; site_spec; n_pairs; n_cation; n_anion

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_coordination

Purpose: Coordination number integrated to rmax. Declared in include/seams.h.

Parameters: sys; rmax; bins; type_i; type_j; cn

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_density_axis

Purpose: Number density along one cell axis. Declared in include/seams.h.

Parameters: sys; type; bins; axis; z; rho; cap; n_out

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_domains

Purpose: Largest domain of a polar or apolar site subset. Declared in include/seams.h.

Parameters: sys; cutoff; site_spec; subset; n_subset; largest; percolation

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint

Purpose: Topological fingerprint of the bonded graph. Declared in include/seams.h.

Parameters: sys; cutoff; type; k; graph; hops; colour_types; out

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: On success the caller owns the new seams_topo and frees it with the matching destroy function. Buffers the caller passed stay with the caller.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_atom

Purpose: Local key of one atom. Declared in include/seams.h.

Parameters: fp; index; buf; cap; n_out

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_destroy

Purpose: Free a fingerprint. Declared in include/seams.h.

Parameters: fp

Returns: Nothing.

Ownership: Frees the caller-owned handle. delete of a null pointer is a no-op.

Errors: None.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_key

Purpose: Frame key text. Declared in include/seams.h.

Parameters: fp; buf; cap; n_out

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_rings

Purpose: Primitive-ring census stored on the fingerprint. Declared in include/seams.h.

Parameters: fp; counts; counts_len

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_get_num_threads

Purpose: Team width the next analysis call will request. Declared in include/seams.h.

Parameters: none

Returns: A positive count.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_graph

Purpose: Neighbour graph. Seeded is a cage graph: mutual seeds, union completions. Declared in include/seams.h.

Parameters: none

Returns: The enumeration type.

Ownership: No allocation.

Errors: None. This is a type.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_handle_error

Purpose: Detail string stored on a handle. Declared in include/seams.h.

Parameters: handle

Returns: The handle's string, or the thread-local string when handle is NULL.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/c_api_tour.c`

### seams_hbonds

Purpose: Hydrogen-bond count on one frame. Declared in include/seams.h.

Parameters: sys; cutoff; oxygen_type; hydrogen_type; dist_cutoff; angle_cutoff; n_bonds

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_label

Purpose: CHILL and CHILL+ atom labels. Matches molSys::atom_state_type. Declared in include/seams.h.

Parameters: none

Returns: The enumeration type.

Ownership: No allocation.

Errors: None. This is a type.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_label_name

Purpose: Name of a CHILL label. Declared in include/seams.h.

Parameters: label

Returns: Static string.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/c_api_tour.c`

### seams_lammps_bounds

Purpose: Fill ox, oy, oz, lx, ly, lz, xy, xz, yz from LAMMPS box bounds. Declared in include/seams_lammps_box.h.

Parameters: lo; hi; tilt; box9

Returns: Nothing. The nine cell numbers are written to box9.

Ownership: The caller owns lo, hi, tilt, and box9. The function does not allocate.

Errors: None. The caller supplies three-vectors and a nine-number buffer.

Units: lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_last_error

Purpose: Thread-local detail of the last failing call on this thread. Declared in include/seams.h.

Parameters: none

Returns: Pointer to a thread-local string.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/c_api_tour.c`

### seams_monolayer

Purpose: Primitive-ring histogram of a monolayer. Declared in include/seams.h.

Parameters: sys; cutoff; type; max_depth; sheet_area; counts; counts_len

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_neighbours

Purpose: One neighbour list owned by the caller, built from a frame. Declared in include/seams.h.

Parameters: none

Returns: The opaque type.

Ownership: The caller creates the object and frees it with the matching destroy function.

Errors: None. This is a type.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_build

Purpose: Build a cutoff or mutual neighbour list. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; out

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller frees *out with seams_neighbours_destroy().

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_copy

Purpose: Copy the neighbour CSR into caller buffers. Declared in include/seams.h.

Parameters: nb; offsets; offsets_cap; n_offsets; cols; cols_cap; nnz

Returns: SEAMS_OK when every non-NULL buffer was filled, or when both buffers were NULL and the lengths were written.

Ownership: The caller owns the buffers.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_count

Purpose: Number of neighbour rows. Declared in include/seams.h.

Parameters: nb

Returns: The row count, or 0 when nb is NULL.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_csr

Purpose: Borrow the neighbour list as CSR. Declared in include/seams.h.

Parameters: nb; offsets; n_offsets; cols; nnz

Returns: SEAMS_OK, or SEAMS_ERR_NULL.

Ownership: The pointers are owned by nb and remain valid until seams_neighbours_destroy(nb).

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_destroy

Purpose: Free a neighbour list. Declared in include/seams.h.

Parameters: nb

Returns: Nothing.

Ownership: Frees the caller-owned handle. delete of a null pointer is a no-op.

Errors: None.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_knn

Purpose: K-nearest neighbour list. Declared in include/seams.h.

Parameters: sys; k; cutoff; type; mutual; out

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: On success the caller owns the new seams_neighbours and frees it with the matching destroy function. Buffers the caller passed stay with the caller.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_path_nframes

Purpose: Number of frames in a trajectory path. Declared in include/seams.h.

Parameters: path; nframes

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_INTERNAL.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_prisms

Purpose: Perfect and imperfect prism counts by ring size. Declared in include/seams.h.

Parameters: sys; cutoff; type; perfect; imperfect; len

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_q4q6

Purpose: Local q4 and q6 on the cutoff shell. Declared in include/seams.h.

Parameters: sys; cutoff; type; q4; q6

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns q4 and q6.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_q6bar

Purpose: Lechner-Dellago q6-bar and the ice count above a cut. Declared in include/seams.h.

Parameters: sys; cutoff; type; ice_cut; q6bar; q6bar_len; mean; median; n_ice

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns q6bar.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_ql

Purpose: Local Steinhardt q_l, and the Lechner-Dellago average on request. Declared in include/seams.h.

Parameters: sys; cutoff; type; order_l; graph; ql; ql_len; qlbar; qlbar_len

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns both arrays.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_rdf

Purpose: Partial radial distribution. Declared in include/seams.h.

Parameters: sys; rmax; bins; type_i; type_j; r; g; count; cap; n_out; volume

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer. A NULL buffer with a length pointer queries the required count and returns SEAMS_OK. A short buffer returns SEAMS_ERR_RANGE and writes nothing.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_ring_symmetry

Purpose: Ring topology. Size classes are used when a six-ring is not a cage face. Declared in include/seams.h.

Parameters: none

Returns: The enumeration type.

Ownership: No allocation.

Errors: None. This is a type.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_ring_symmetry_name

Purpose: Name of a ring-symmetry code. Declared in include/seams.h.

Parameters: symmetry

Returns: Static string.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/c_api_tour.c`

### seams_rings

Purpose: One primitive-ring set owned by the caller. Declared in include/seams.h.

Parameters: none

Returns: The opaque type.

Ownership: The caller creates the object and frees it with the matching destroy function.

Errors: None. This is a type.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_rings_count

Purpose: Number of primitive rings. Declared in include/seams.h.

Parameters: rings

Returns: The ring count, or 0 when rings is NULL.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_rings_destroy

Purpose: Free a ring set. Declared in include/seams.h.

Parameters: rings

Returns: Nothing.

Ownership: Frees the caller-owned handle. delete of a null pointer is a no-op.

Errors: None.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_rings_find

Purpose: Primitive rings (Franzblau) up to a maximum size. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; max_depth; out

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller frees *out with seams_rings_destroy().

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_rings_fraction

Purpose: Fraction of atoms that sit on a ring. Declared in include/seams.h.

Parameters: rings; size; fraction

Returns: SEAMS_OK, or SEAMS_ERR_NULL.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_rings_histogram

Purpose: Count primitive rings by size. Declared in include/seams.h.

Parameters: rings; counts; counts_len

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: The caller owns counts.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_rings_labels

Purpose: Size and symmetry of each ring. Declared in include/seams.h.

Parameters: rings; sizes; symmetry

Returns: SEAMS_OK, or SEAMS_ERR_NULL.

Ownership: The caller owns both arrays.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_rings_vertices

Purpose: Atom indices of one ring. Declared in include/seams.h.

Parameters: rings; index; atoms; cap; n_out

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: The caller owns atoms.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_running_cn

Purpose: Running coordination number along the partial RDF. Declared in include/seams.h.

Parameters: sys; rmax; bins; type_i; type_j; r; cn; cap; n_out

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_set_num_threads

Purpose: Set the OpenMP team used by later analysis calls. Declared in include/seams.h.

Parameters: n

Returns: SEAMS_OK, or SEAMS_ERR_RANGE when n is negative.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_signature_cages

Purpose: Closed cages of one ring census. Declared in include/seams.h.

Parameters: sys; cutoff; type; graph; spec; n_cages; n_atoms

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_site_kind

Purpose: Site kind. Values match the engine's site table. Declared in include/seams.h.

Parameters: none

Returns: The enumeration type.

Ownership: No allocation.

Errors: None. This is a type.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_site_table

Purpose: Per-atom site kind from a site-table spec. Declared in include/seams.h.

Parameters: sys; spec; kinds; len

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_status

Purpose: Status codes. Zero is success. Declared in include/seams.h.

Parameters: none

Returns: The enumeration type.

Ownership: No allocation.

Errors: None. This is a type.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_status_string

Purpose: Short name of a status code. Declared in include/seams.h.

Parameters: status

Returns: Pointer to a static string.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/c_api_tour.c`

### seams_steinhardt

Purpose: Local Steinhardt q_l on the cutoff shell, and the average on request. Declared in include/seams.h.

Parameters: sys; cutoff; type; order_l; ql; ql_len; qlbar; qlbar_len

Returns: Seams_status. SEAMS_OK (0) is success.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: cutoff and lengths in angstrom

Example: `docs/examples/c_api_tour.c`

### seams_system

Purpose: One frame. The caller creates and destroys it. Declared in include/seams.h.

Parameters: none

Returns: The opaque type.

Ownership: The caller creates the object and frees it with the matching destroy function.

Errors: None. This is a type.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_system_count

Purpose: Number of points in the frame. Declared in include/seams.h.

Parameters: sys

Returns: The atom count, or 0 when sys is NULL or holds no frame.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_system_create

Purpose: Allocate an empty frame. Declared in include/seams.h.

Parameters: none

Returns: A new handle, or NULL when allocation fails.

Ownership: The caller owns the handle and frees it with the matching destroy function.

Errors: NULL. seams_last_error() then reads 'out of memory'.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_system_destroy

Purpose: Free a frame. Declared in include/seams.h.

Parameters: sys

Returns: Nothing.

Ownership: Frees the caller-owned handle. delete of a null pointer is a no-op.

Errors: None.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_system_read

Purpose: Read one frame from a trajectory file. Declared in include/seams.h.

Parameters: sys; path; frame; type_filter

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The library copies the frame out of the file.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_system_set_atoms

Purpose: Store points and a general cell. Declared in include/seams.h.

Parameters: sys; xyz; n; box9; box_len; types; ids

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The library copies the arrays.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_system_set_frame

Purpose: Store points and a general LAMMPS cell. Declared in include/seams.h.

Parameters: sys; xyz; n; box9; types; ids

Returns: Same codes as seams_system_set_atoms().

Ownership: The library copies the arrays.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_system_set_ortho

Purpose: Store points and an orthorhombic cell. Declared in include/seams.h.

Parameters: sys; xyz; n; ox; oy; oz; lx; ly; lz; types; ids

Returns: Same codes as seams_system_set_frame().

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_system_symmetry

Purpose: Per-atom ice shell label. Declared in include/seams.h.

Parameters: sys; shell_rmsd; labels; labels_len

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_INTERNAL, SEAMS_ERR_NOMEM.

Ownership: The caller owns labels.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Units: as named by the parameters

Example: `docs/examples/c_api_tour.c`

### seams_topo

Purpose: Topological fingerprint of one bonded frame. Declared in include/seams.h.

Parameters: none

Returns: The opaque type.

Ownership: The caller creates the object and frees it with the matching destroy function.

Errors: None. This is a type.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_version

Purpose: Packed release number. Declared in include/seams.h.

Parameters: none

Returns: SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_version_major

Purpose: Major version. Declared in include/seams.h.

Parameters: none

Returns: The integer result.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_version_minor

Purpose: Minor version. Declared in include/seams.h.

Parameters: none

Returns: The integer result.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_version_patch

Purpose: Patch version. Declared in include/seams.h.

Parameters: none

Returns: The integer result.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/c_api_tour.c`

### seams_version_string

Purpose: Version text, for example "2.14.0". Declared in include/seams.h.

Parameters: none

Returns: Pointer to a static string.

Ownership: The library owns the string. The caller does not free it.

Errors: None.

Units: text

Example: `docs/examples/c_api_tour.c`
<!-- API_SYMBOLS_END -->

The program those lines come from:

```{literalinclude} ../../../docs/examples/c_api_tour.c
```


## Structure still

```{figure} ../tutorials/molecules/c/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug c` writes the 2400 by 1800 PNG.
```

