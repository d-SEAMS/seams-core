# C API

`seams.h` is ABI 1 of libseams. A handle is opaque. Coordinates are Angstrom.
A cell is nine numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`. A zero tilt is orthorhombic.
`lx, ly, lz` are edge lengths.

`docs/examples/minimal_c.c` is the short program. `docs/examples/c_api_tour.c` calls every exported function on the in-tree fixtures.
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

No C API call aborts or exits the process. Every C++ exception is caught inside the
library and returned as one of these codes.

## Length-checked entry points

The ABI 1 calls `seams_chill_plus`, `seams_q4q6`, and `seams_rings_labels` write into
buffers without a length, and `seams_rings_histogram` takes an `int` length. New code
calls `seams_chill_plus_n`, `seams_q4q6_n`, `seams_rings_labels_n`, and
`seams_rings_histogram_n`, which take a `size_t` length for every buffer and return
`SEAMS_ERR_RANGE` when one is short. `seams_system_natoms`, `seams_neighbours_nrows`,
and `seams_rings_nrings` return their count through a `size_t` pointer and report a
null or empty handle as a status. The older calls stay exported with the same
behaviour.

## Linking

libseams is compiled with hidden visibility and a linker export list, so the only
exported symbols are the `seams_*` functions in `seams.h`. `SEAMS_API` marks them.
The library build defines `SEAMS_BUILDING`, which selects `dllexport` on Windows;
consumers get `dllimport`. Define `SEAMS_STATIC` before including `seams.h` to link
a static archive. `tests/capi_symbols.txt` is the exported list, and the
`capi_symbols` test compares it with the header and the built library.

## seams_last_error and seams_handle_error

`seams_last_error_message` and its ABI 1 spelling `seams_last_error` return the same
string. `seams_last_error` is the detail of the last failing call on this thread. A later
successful call on the same thread clears it. The storage is thread-local.
`seams_handle_error` is the copy stored on the handle that failed. A call that takes
the handle non-const clears that copy on entry. A call that takes it `const` writes the
copy only when it fails, so a successful `const` call leaves the last failure in place.
The copy stays until the next failing call or non-const call on that handle, or until
destroy. Pass `NULL` and
`seams_handle_error` returns the thread-local string. The caller does not free either
pointer.

On an empty frame, `seams_chill_plus` returns `SEAMS_ERR_STATE`. The tour prints:

```
capi seams_last_error system has no frame
capi seams_handle_error no-frame
```

## Const handles

A function that takes a `const` handle is a read-only accessor. On success it writes
nothing to the handle; on failure it writes only the error copy, under a lock. Several
threads may make `const` calls on one handle at the same time, as long as no non-const
call on that handle overlaps them. Read the message with `seams_last_error_message`
when threads share a handle. `tests/capi_const_threads.c` runs eight threads this way.

The frame accessors are `seams_system_count`, `seams_system_natoms`, `seams_system_box`,
`seams_system_positions`, `seams_system_types`, and `seams_system_ids`. Every query on
a neighbour, ring, cage, brine, or fingerprint handle is also `const`. The ring label
cache that `seams_rings_labels` fills on the first symmetry request is built once,
under a lock.

The analyses on a `seams_system` stay non-const, and each declaration in `seams.h`
says why under "Handle access":

- the loads replace the stored frame;
- `seams_chill`, `seams_chill_plus`, and `seams_chill_plus_n` write CHILL+ labels into
  the stored frame;
- `seams_ql`, `seams_q4q6`, `seams_q4q6_n`, `seams_q6bar`, and `seams_steinhardt` use
  process-global q_lm scratch, so two of them must not overlap even on different
  handles;
- the others keep no scratch on the handle, and are non-const only because they clear
  the error copy on entry, which is what lets `seams_handle_error` report that call.

## Reloading a frame

`seams_system_set_atoms` reuses the storage of the frame it replaces. When the ids
are `NULL`, or 1 to n in order, and the atom count is unchanged, it also keeps the id
index. `seams_system_set_positions` overwrites only the coordinates, and optionally
the cell, for a frame of the same atom count; types and ids stay. A trajectory loop
calls `seams_system_set_atoms` once and `seams_system_set_positions` for each later
frame. There is no borrowed-pointer load. Every analysis reads the library's per-atom
records, so a view onto caller memory would be copied on every analysis call instead
of once per frame. A failed load or reload leaves the earlier frame in place.

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

`seams_chill_plus` writes one label per atom and, when `counts9` is not null, nine
population counts. `rule` is `SEAMS_CHILL_PUBLISHED` (four closest neighbours in
q_lm, Nguyen and Molinero equation 1) or `SEAMS_CHILL_OVITO` (every neighbour
inside the cutoff). `seams_label_name` maps a label to cubic, hexagonal, water,
interfacial, clathrate, interClathrate, unclassified, reCubic, or reHex.

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

### seams_last_error_message

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_last_error_message system has no frame
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

### seams_system_natoms

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_natoms 128
```

### seams_system_set_positions

The call is in `docs/examples/c_api_tour.c`. It moves the eight-atom cell and
keeps its types and ids. The run prints:

```
capi seams_system_set_positions 8
```

### seams_system_box

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_box 0.000 0.000 0.000 6.351 6.351 6.351
```

### seams_system_positions

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_positions 0.250 0.250 0.250
```

### seams_system_types

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_types 8 of type 1
```

### seams_system_ids

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_system_ids 1 8
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

### seams_neighbours_nrows

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_neighbours_nrows 128
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

### seams_chill_plus_n

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_chill_plus_n hexagonal 128
```

### seams_q4q6

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_q4q6 0.625855
```

### seams_q4q6_n

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_q4q6_n 0.625855
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

### seams_rings_nrings

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_nrings 256
```

### seams_rings_histogram

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_histogram 256
```

### seams_rings_histogram_n

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_histogram_n 256
```

### seams_rings_labels

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_labels 6
```

### seams_rings_labels_n

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_labels_n 6
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

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_BUILDING

Purpose: Defined only while libseams itself is compiled; selects dllexport on Windows.

Parameters: none

Units: none

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_BOTH

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_CUBIC

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_HEX

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_CAGE_WATER

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_CHILL_OVITO

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_CHILL_PUBLISHED

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_DEPRECATED

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_ARG

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_BOX

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_INTERNAL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_NOMEM

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_NULL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_RANGE

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_ERR_STATE

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_GRAPH_CUTOFF

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_GRAPH_MUTUAL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_GRAPH_SEEDED

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_CLATHRATE

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_CUBIC

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_HEXAGONAL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_INTERCLATHRATE

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_INTERFACIAL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_RECUBIC

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_REHEX

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_UNCLASSIFIED

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_LABEL_WATER

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_OK

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_BOTH_BASAL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_BOTH_PRISMATIC

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_DDC

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_DEFORMED_PRISM

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_HC_BASAL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_HC_PRISMATIC

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_HEXAGON

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_MIXED_PRISM

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_OTHER

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_PENTAGON

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_PRISM

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_RING_UNCLASSIFIED

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_ACCEPTOR

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_ANION

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_APOLAR

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_CATION

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_DONOR_H

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_POLAR

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_SOLVENT

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_TAIL

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_UNSPECIFIED

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_WATER_H

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_SITE_WATER_O

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_STATIC

Purpose: Define before including seams.h when linking a static libseams.

Parameters: none

Units: none

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_MAJOR

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_MINOR

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_NUMBER

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_PATCH

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### SEAMS_VERSION_STRING

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_abi

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_brine

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_brine_analyse

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_brine_destroy

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_brine_fronts

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_brine_pockets

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_brine_profile

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_brine_rejection

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cage_label

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cage_set

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cages

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cages_counts

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cages_destroy

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cages_find

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cages_labels

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_cages_signature

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_chill

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_chill_plus

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_chill_plus_n

Purpose: CHILL+ labels into a caller buffer of checked length.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_chill_rule

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_contact_pairs

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_coordination

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_density_axis

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_domains

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_atom

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_destroy

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_key

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_rings

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_get_num_threads

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_graph

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_handle_error

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_hbonds

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_label

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_label_name

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_last_error

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_last_error_message

Purpose: Thread-local detail of the last failing call; same string as seams_last_error.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_monolayer

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_build

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_copy

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_count

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_csr

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_destroy

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_knn

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_nrows

Purpose: Row count of a neighbour list, returned through a status.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_path_nframes

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_prisms

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_q4q6

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_q4q6_n

Purpose: Local q4 and q6 into caller buffers of checked length.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_q6bar

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_ql

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rdf

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_ring_symmetry

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_ring_symmetry_name

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_count

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_destroy

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_find

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_fraction

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_histogram

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_histogram_n

Purpose: Ring-size histogram into a size_t-length buffer.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_labels

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_labels_n

Purpose: Size and symmetry of each ring into buffers of checked length.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_nrings

Purpose: Primitive ring count, returned through a status.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_rings_vertices

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_running_cn

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_set_num_threads

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_signature_cages

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_site_kind

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_site_table

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_status

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_status_string

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_steinhardt

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_box

Purpose: Copy the stored cell out as the nine LAMMPS numbers; read-only on a const handle.

Parameters: see the declaration in include/seams.h.

Units: Angstrom

Example: `docs/examples/c_api_tour.c`

### seams_system_count

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_create

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_destroy

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_ids

Purpose: Copy the stored atom ids out; read-only on a const handle.

Parameters: see the declaration in include/seams.h.

Units: none

Example: `docs/examples/c_api_tour.c`

### seams_system_natoms

Purpose: Atom count of a frame, returned through a status.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_positions

Purpose: Copy the stored coordinates out; read-only on a const handle.

Parameters: see the declaration in include/seams.h.

Units: Angstrom

Example: `docs/examples/c_api_tour.c`

### seams_system_read

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_atoms

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_frame

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_ortho

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_positions

Purpose: Overwrite the coordinates, and optionally the cell, of the stored frame in place.

Parameters: see the declaration in include/seams.h.

Units: Angstrom

Example: `docs/examples/c_api_tour.c`

### seams_system_symmetry

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_system_types

Purpose: Copy the stored atom types out; read-only on a const handle.

Parameters: see the declaration in include/seams.h.

Units: none

Example: `docs/examples/c_api_tour.c`

### seams_topo

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_version

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_version_major

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_version_minor

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_version_patch

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

Example: `docs/examples/c_api_tour.c`

### seams_version_string

Purpose: Public name in include/seams.h.

Parameters: see the declaration in include/seams.h.

Units: as documented on that declaration.

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

