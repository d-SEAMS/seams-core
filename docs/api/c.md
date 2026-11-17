# C

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

One thread uses a handle. Another thread may use a different handle at the same time,
including for `seams_ql`, `seams_q4q6`, `seams_q6bar`, and `seams_steinhardt`: their
q_lm scratch belongs to the calling thread. `tests/capi_q_threads.c` runs them, and
CHILL+, on eight handles in eight threads.

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
behavior.

## Linking

libseams is compiled with hidden visibility and a linker export list, so the only
exported symbols are the `seams_*` functions in `seams.h`. `SEAMS_API` marks them.
The library build defines `SEAMS_BUILDING`, which selects `dllexport` on Windows;
consumers get `dllimport`. Define `SEAMS_STATIC` before including `seams.h` to link
a static archive. `tests/capi_symbols.txt` is the exported list, and the
`capi_symbols` test compares it with the header and the built library.
`capi_abi` compares the exported names with `scripts/abi_baseline.txt`. When
libabigail's `abidiff` is on `PATH` at configure time, it also builds libseams once
more with debug info and compares its types with the ABI 1 corpus in
`scripts/abi_baseline.abi`: added functions pass, and anything else fails. The
`abi` pixi environment carries conda-forge's libabigail, so
`pixi run -e abi meson setup <builddir>` finds it, and the `capi-abi` CI job runs
that way. Distribution packages work too: `abigail-tools` on Debian and Ubuntu,
`libabigail` on Fedora. The corpus is format 2.4, which libabigail 2.9 and 3.0
read. Without
`abidiff`, or with an older one, `capi_abi` still checks the exported names, says
the type step was skipped, and passes.

The 2.13.0 release had no libseams. Its C interface was `seams_lang.h` in
libseams_lang, which the R and Julia packages still use. `capi_abi_2_13` builds
libseams_lang with debug info and compares its `seams_*` functions with the 2.13.0
corpus in `scripts/abi_2.13.0_seams_lang.abi`. Functions added since 2.13.0 pass,
and a removed function or a changed signature fails. `capi_abi_2_13_selftest`
edits copies of the corpus to prove the comparison catches both. Without
`abidiff`, both tests are skipped.

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
- the others, including `seams_ql`, `seams_q4q6`, `seams_q4q6_n`, `seams_q6bar`, and
  `seams_steinhardt`, keep no scratch on the handle, and are non-const only because
  they clear the error copy on entry, which is what lets `seams_handle_error` report
  that call. The q_lm scratch of the q calls belongs to the calling thread.

## Reloading a frame

`seams_system_set_atoms` reuses the storage of the frame it replaces. When the ids
are `NULL`, or 1 to n in order, and the atom count is unchanged, it also keeps the id
index. `seams_system_set_positions` overwrites only the coordinates, and optionally
the cell, for a frame of the same atom count; types and ids stay. A trajectory loop
calls `seams_system_set_atoms` once and `seams_system_set_positions` for each later
frame. There is no borrowed-pointer load. Every analysis reads the library's per-atom
records, so a view onto caller memory would be copied on every analysis call instead
of once per frame. A failed load or reload leaves the earlier frame in place.
`seams_system_set_positions` keeps that guarantee by checking every coordinate
before it writes any, so it reads the input twice and needs no extra memory. The
check scans from the end of the input, so the write pass finds its start still in
cache. In a trajectory loop a reload then costs 2 to 5% more than writing the
coordinates into the C++ frame directly at 4096 and 110592 molecules, and about 12%
more at 884736, where the second read of the 21 MB input comes from L3.
`seams_system_stage_positions(sys, 1)` makes it read the input once instead: it
writes into a second set of per-atom records, checked in the same pass, and swaps
that set in at the end. That reload stays within 5% of the direct write at each of
those sizes. The set costs about 72 bytes per atom, 15 to 18% of the resident
memory of a process that has loaded an ice frame of 0.1 to 0.9 million molecules
and run CHILL+ on it.

For large frames staged mode is the zero-overhead path. Once the input no longer
fits in L2, the default check costs a full second read of it, and staging is the
only way to read it once. The times below are per reload in a trajectory loop: the
median of 21 batches of 100 calls on one thread, from `bench_capi_overhead`. The
direct write is the C++ loop that writes x, y, z and the label into a built frame.

| molecules | direct write | default | staged |
| --- | --- | --- | --- |
| 4096 | 0.0068 ms | 0.0070 ms (+2.3%) | 0.0072 ms (+4.4%) |
| 110592 | 0.3882 ms | 0.4067 ms (+4.8%) | 0.3664 ms (-5.6%) |
| 884736 | 3.1568 ms | 3.5415 ms (+12.2%) | 3.1939 ms (+1.2%) |

Each interquartile range is at most 2.2% of its median. A loop over frames of a
million atoms or more should call `seams_system_stage_positions(sys, 1)` once after
the first `seams_system_set_atoms`, if the extra 72 bytes per atom fit.

## Threads

`seams_set_num_threads(1)` selects one OpenMP thread. `seams_set_num_threads(0)`
restores the library default. A negative width returns `SEAMS_ERR_RANGE`.
`seams_get_num_threads` reports the width the next call will use.
`seams_set_num_threads` is process-global. Do not call it while another thread is
inside `seams_get_num_threads`.

## Examples by family

Each program below covers one family and checks its own results on
`input/traj/ice1h.lammpstrj`. `meson test capi_example_<family>` compiles it as C11
with `-Wall -Wextra -Wpedantic -Werror` and runs it. The argument is the repository root.
`capi_installed_examples` installs the build into a staging directory and builds all
six out of tree against that install, once with the flags from `pkg-config seams` and
once with `find_package(seams)` and `docs/examples/c/CMakeLists.txt`. Each binary must
load the installed `libseams.so` and pass.

Neighbours, `docs/examples/c/neighbours.c`:

```{literalinclude} ../../../docs/examples/c/neighbours.c
:language: c
```

Rings, `docs/examples/c/rings.c`:

```{literalinclude} ../../../docs/examples/c/rings.c
:language: c
```

CHILL+, `docs/examples/c/chill_plus.c`:

```{literalinclude} ../../../docs/examples/c/chill_plus.c
:language: c
```

q4 and q6, `docs/examples/c/q4q6.c`:

```{literalinclude} ../../../docs/examples/c/q4q6.c
:language: c
```

Cages, `docs/examples/c/cages.c`:

```{literalinclude} ../../../docs/examples/c/cages.c
:language: c
```

Errors, `docs/examples/c/errors.c`, needs no input file:

```{literalinclude} ../../../docs/examples/c/errors.c
:language: c
```

## Documentation check

`capi_docs_complete` checks that every function in `tests/capi_symbols.txt` has a
doc comment in `seams.h` and that its section below carries the same text.
`python3 docs/capi_docs.py --write include/seams.h tests/capi_symbols.txt docs/api/c.md`
rewrites those sections from the header.

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

### seams_system_stage_positions

The call is in `docs/examples/c_api_tour.c`. It reloads the cell once through the
second set of records and then frees that set. The run prints:

```
capi seams_system_stage_positions 8
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

### seams_rings_chordless

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_rings_chordless 256
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

### seams_signature_cages_edge

The call is in `docs/examples/c_api_tour.c`. The run prints:

```
capi seams_signature_cages_edge 2 40
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

### SEAMS_BUILDING

Purpose: Defined only while libseams itself is compiled; selects dllexport on Windows.

Parameters: none

Returns: The result declared for this name in include/seams.h.

Ownership: The caller owns any buffer it allocated. Library text is not freed by the caller.

Errors: None.

Units: none

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

### SEAMS_STATIC

Purpose: Define before including seams.h when linking a static libseams.

Parameters: none

Returns: The result declared for this name in include/seams.h.

Ownership: The caller owns any buffer it allocated. Library text is not freed by the caller.

Errors: None.

Units: none

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

Purpose: ABI number compiled into this library.

```c
int seams_abi(void);
```

Parameters: none.

Returns: SEAMS_ABI, which is 1.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
if (seams_abi() != SEAMS_ABI) { return 1; }
```

The tour call is in `docs/examples/c_api_tour.c`.

Ownership: No allocation.

Errors: None.

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

Purpose: Planar ice/brine front from an ice assignment.

ice[i] is non-zero when atom i is ice. water_type 0 keeps every atom whose type is not ion_type. ion_type 0 means the frame has no ions. n_bins 0 selects 20 bins. axis is 0, 1, or 2. The ice fraction cut is 0.5, and the bulk ice and liquid cuts are 0.85 and 0.15.

```c
seams_status seams_brine_analyse(seams_system *sys, const int *ice, size_t n_ice, int water_type, int ion_type, int n_bins, int axis, seams_brine **out);
```

Parameters: `sys`, `ice`, `n_ice`, `water_type`, `ion_type`, `n_bins`, `axis`, `out`.

- `sys`: Frame handle.
- `ice`: Per-atom ice flags, length n_ice.
- `n_ice`: Must equal the frame count.
- `water_type`: Species kept as water, or 0.
- `ion_type`: Species kept as ions, or 0 for none.
- `n_bins`: Bin count along axis, or 0 for 20.
- `axis`: 0, 1, or 2.
- `out`: Receives a new handle. On failure *out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller frees *out with seams_brine_destroy().

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_brine *front = NULL;
if (seams_brine_analyse(sys, ice, n, 1, 2, 0, 2, &front) == SEAMS_OK) {
  seams_brine_destroy(front);
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_brine_destroy

Purpose: Free a brine-front handle.

```c
void seams_brine_destroy(seams_brine *front);
```

Parameters: `front`.

- `front`: Handle from seams_brine_analyse(), or NULL. NULL is ignored.

Returns: nothing.

Ownership: Releases the handle and its profile.

Thread safety: The caller must not use front after this returns.

Units: as given for each parameter above.

Example:

```c
seams_brine_destroy(front);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_brine_fronts

Purpose: The two planar crossings of the ice slab.

Lengths are Angstrom, in the frame's cell.

```c
seams_status seams_brine_fronts(const seams_brine *front, double *z_enter, double *z_exit, double *thickness, int *found);
```

Parameters: `front`, `z_enter`, `z_exit`, `thickness`, `found`.

- `front`: Brine handle.
- `z_enter`: Where the ice fraction crosses the cut on entering the ice along the positive axis, or NULL.
- `z_exit`: The crossing on leaving the ice, or NULL.
- `thickness`: Ice run along the axis, or NULL.
- `found`: Receives 1 when both crossings exist, otherwise 0, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL or every output pointer is NULL.

Ownership: The caller owns front; the call borrows it. z_enter, z_exit, thickness and found are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on front, and not alongside a non-const call on front.

Units: as given for each parameter above.

Example:

```c
double thick = 0;
int found = 0;
seams_brine_fronts(front, NULL, NULL, &thick, &found);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_brine_pockets

Purpose: Brine pockets trapped in the ice, excluding the liquid reservoir.

```c
seams_status seams_brine_pockets(const seams_brine *front, int *n_pockets, int *reservoir_waters);
```

Parameters: `front`, `n_pockets`, `reservoir_waters`.

- `front`: Brine handle.
- `n_pockets`: Receives the pocket count, or NULL.
- `reservoir_waters`: Receives the water count of the liquid reservoir, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL or both output pointers are NULL.

Ownership: The caller owns front; the call borrows it. n_pockets and reservoir_waters are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on front, and not alongside a non-const call on front.

Units: as given for each parameter above.

Example:

```c
int pockets = 0;
seams_brine_pockets(front, &pockets, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_brine_profile

Purpose: Concentration profile along the front axis.

z_mid is Angstrom. molality is mol/kg. ice_fraction is dimensionless.

```c
seams_status seams_brine_profile(const seams_brine *front, double *z_mid, double *molality, double *ice_fraction, size_t cap, size_t *n_out);
```

Parameters: `front`, `z_mid`, `molality`, `ice_fraction`, `cap`, `n_out`.

- `front`: Brine handle.
- `z_mid`: Bin centers, or NULL to skip the column.
- `molality`: Bin molalities, or NULL to skip the column.
- `ice_fraction`: Bin ice fractions, or NULL to skip the column.
- `cap`: Capacity of each non-NULL column.
- `n_out`: Receives the bin count.

Returns: SEAMS_OK when every non-NULL column was filled, or when every column pointer was NULL and n_out was written. SEAMS_ERR_RANGE when a non-NULL column is short. SEAMS_ERR_NULL when front or n_out is NULL.

Ownership: The caller owns the columns.

Thread safety: Read-only: may run alongside other const calls on front, and not alongside a non-const call on front.

Units: as given for each parameter above.

Example:

```c
size_t bins = 0;
seams_brine_profile(front, NULL, NULL, NULL, 0, &bins);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_brine_rejection

Purpose: Molal rejection coefficient k = c_ice / c_liquid.

Molality is formula units per kilogram of water.

```c
seams_status seams_brine_rejection(const seams_brine *front, double *c_ice, double *c_liquid, double *k, int *defined);
```

Parameters: `front`, `c_ice`, `c_liquid`, `k`, `defined`.

- `front`: Brine handle.
- `c_ice`: Bulk ice molality, or NULL.
- `c_liquid`: Bulk liquid molality, or NULL.
- `k`: The ratio, or NULL. Unchanged when the ratio is undefined.
- `defined`: Receives 1 when the liquid bulk contains salt or the ice bulk does not, otherwise 0, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL or every output pointer is NULL.

Ownership: The caller owns front; the call borrows it. c_ice, c_liquid, k and defined are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on front, and not alongside a non-const call on front.

Units: as given for each parameter above.

Example:

```c
double k = 0;
int defined = 0;
seams_brine_rejection(front, NULL, NULL, &k, &defined);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Hexagonal and double-diamond cages. Same search as seams_cages_find().

```c
seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `k`, `complete`, `out`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `graph`: SEAMS_GRAPH_CUTOFF, SEAMS_GRAPH_MUTUAL, or SEAMS_GRAPH_SEEDED.
- `k`: Neighbour count for the seeded graph. Zero selects 4.
- `complete`: Non-zero runs ring-adjacent completion on seeded labels.
- `out`: Receives a ::seams_cage_set. On failure *out is NULL.

Returns: Same codes as seams_cages_find().

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller frees *out with seams_cages_destroy().

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_cage_set *cages = NULL;
seams_cages(sys, 3.5, 0, SEAMS_GRAPH_SEEDED, 4, 1, &cages);
seams_cages_destroy(cages);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_cages_counts

Purpose: Disjoint cage population counts.

The four counts sum to the number of atoms the search kept.

```c
seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex, int *n_cubic, int *n_water, int *n_both);
```

Parameters: `cages`, `n_hex`, `n_cubic`, `n_water`, `n_both`.

- `cages`: Cage handle.
- `n_hex`: Hexagonal only, or NULL.
- `n_cubic`: Double-diamond only, or NULL.
- `n_water`: Neither cage, or NULL.
- `n_both`: Both cages, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL when cages is NULL or every count pointer is NULL.

Ownership: The caller owns cages; the call borrows it. n_hex, n_cubic, n_water and n_both are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on cages, and not alongside a non-const call on cages.

Units: as given for each parameter above.

Example:

```c
int hex = 0, cubic = 0;
seams_cages_counts(cages, &hex, &cubic, NULL, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_cages_destroy

Purpose: Free a cage-membership handle.

```c
void seams_cages_destroy(seams_cage_set *cages);
```

Parameters: `cages`.

- `cages`: Handle from seams_cages_find() or seams_cages(), or NULL. NULL is ignored.

Returns: nothing.

Ownership: Releases the handle and its per-atom labels.

Thread safety: The caller must not use cages after this returns.

Units: as given for each parameter above.

Example:

```c
seams_cages_destroy(cages);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_cages_find

Purpose: Hexagonal and double-diamond cage membership.

```c
seams_status seams_cages_find(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `k`, `complete`, `out`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom. The seeded graph uses cutoff + 1.5 as the candidate shell.
- `type`: Species filter.
- `graph`: SEAMS_GRAPH_CUTOFF, SEAMS_GRAPH_MUTUAL, or SEAMS_GRAPH_SEEDED.
- `k`: Neighbour count for the seeded graph. Zero selects 4.
- `complete`: Non-zero runs ring-adjacent completion on the seeded labels. Ignored for the cutoff and mutual graphs.
- `out`: Receives a new handle. On failure *out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller frees *out with seams_cages_destroy().

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_cage_set *cages = NULL;
if (seams_cages_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 0, 0, &cages) == SEAMS_OK) {
  seams_cages_destroy(cages);
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_cages_labels

Purpose: Per-atom cage code.

```c
seams_status seams_cages_labels(const seams_cage_set *cages, int *per_atom, size_t len);
```

Parameters: `cages`, `per_atom`, `len`.

- `cages`: Cage handle.
- `per_atom`: Caller buffer of ::seams_cage_label values, one per atom in the frame. Atoms the species filter dropped are SEAMS_CAGE_WATER.
- `len`: Length of per_atom. Must be at least the frame count stored on the handle.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: The caller owns per_atom.

Thread safety: Read-only: may run alongside other const calls on cages, and not alongside a non-const call on cages.

Units: as given for each parameter above.

Example:

```c
seams_cages_labels(cages, per_atom, n);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_cages_signature

Purpose: Count closed cages whose ring census matches a signature.

spec is a census such as "512", "51262", or "51264", or the names "hc" and "ddc". n_cages is the number of closed cages. n_atoms is the number of distinct vertices in those cages.

```c
seams_status seams_cages_signature(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `spec`, `n_cages`, `n_atoms`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `spec`: Signature text. Must be non-NULL.
- `n_cages`: Receives the cage count.
- `n_atoms`: Receives the vertex count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns sys; the call borrows it. spec is read during the call only. n_cages and n_atoms are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
size_t n_cages = 0, n_atoms = 0;
seams_cages_signature(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", &n_cages, &n_atoms);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_chill

Purpose: Original CHILL labels, published or OVITO edges.

Published CHILL is Moore's four-neighbour rule with inclusive c_ij edges. OVITO mode uses the same thresholds and the open edges. Labels match ::seams_label.

```c
seams_status seams_chill(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, size_t labels_len, int counts9[9]);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `rule`, `labels`, `labels_len`, `counts9`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter. Zero and negative keep every atom.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `rule`: SEAMS_CHILL_PUBLISHED or SEAMS_CHILL_OVITO.
- `labels`: Caller buffer of labels_len ints.
- `labels_len`: Must be at least the frame count.
- `counts9`: Nine ints that receive the population of each label, or NULL to skip the histogram.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const: writes each atom's CHILL+ label into the stored frame and clears the handle's error slot on entry. Coordinates, cell, types, and ids are not changed.

Ownership: The caller owns labels and counts9.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
int counts[9];
seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n, counts);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_chill_plus

Purpose: CHILL+ labels, without a buffer length.

ABI 1 entry. labels must hold one int per atom in the frame. New callers use seams_chill_plus_n(), which checks the length.

```c
seams_status seams_chill_plus(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, int counts9[9]);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `rule`, `labels`, `counts9`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter. Zero and negative keep every atom.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `rule`: SEAMS_CHILL_PUBLISHED or SEAMS_CHILL_OVITO.
- `labels`: Caller buffer of one label per atom.
- `counts9`: Nine ints, or NULL.

Returns: Same codes as seams_chill_plus_n().

Handle access: Non-const: writes each atom's CHILL+ label into the stored frame and clears the handle's error slot on entry. Coordinates, cell, types, and ids are not changed.

Ownership: The caller owns labels and counts9.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_chill_plus_n

Purpose: CHILL+ labels into a buffer of known length.

rule selects the q_lm shell. SEAMS_CHILL_PUBLISHED, the default, is Nguyen and Molinero equation 1: the four closest neighbours. SEAMS_CHILL_OVITO is every neighbour inside the cutoff, which is the OVITO Chill+ modifier. Labels match ::seams_label.

```c
seams_status seams_chill_plus_n(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, size_t labels_len, int counts9[9]);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `rule`, `labels`, `labels_len`, `counts9`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter. Zero and negative keep every atom.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `rule`: SEAMS_CHILL_PUBLISHED or SEAMS_CHILL_OVITO.
- `labels`: Caller buffer of labels_len ints.
- `labels_len`: Must be at least the frame count.
- `counts9`: Nine ints that receive the population of each label, or NULL.

Returns: SEAMS_OK; SEAMS_ERR_NULL for a NULL sys or labels; SEAMS_ERR_STATE when sys holds no frame; SEAMS_ERR_RANGE for a short labels buffer or a non-positive or non-finite cutoff; SEAMS_ERR_ARG for an unknown rule or graph; SEAMS_ERR_NOMEM; SEAMS_ERR_INTERNAL.

Handle access: Non-const: writes each atom's CHILL+ label into the stored frame and clears the handle's error slot on entry. Coordinates, cell, types, and ids are not changed.

Ownership: The caller owns labels and counts9.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
int counts[9];
seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                   labels, n, counts);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Mutual nearest unlike contacts of a site table.

Cation and anion counts are molecules in the ion cloud.

```c
seams_status seams_contact_pairs(seams_system *sys, const char *site_spec, size_t *n_pairs, int *n_cation, int *n_anion);
```

Parameters: `sys`, `site_spec`, `n_pairs`, `n_cation`, `n_anion`.

- `sys`: Frame handle.
- `site_spec`: Text such as "1=cationHead,2=anion". Must be non-NULL.
- `n_pairs`: Receives the contact count, or NULL.
- `n_cation`: Receives the cation molecule count, or NULL.
- `n_anion`: Receives the anion molecule count, or NULL. At least one of the three outputs must be non-NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns sys; the call borrows it. site_spec is read during the call only. n_pairs, n_cation and n_anion are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
size_t pairs = 0;
seams_contact_pairs(sys, "1=cationHead,2=anion", &pairs, NULL, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_coordination

Purpose: Coordination number integrated to rmax.

```c
seams_status seams_coordination(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *cn);
```

Parameters: `sys`, `rmax`, `bins`, `type_i`, `type_j`, `cn`.

- `sys`: Frame handle.
- `rmax`: Positive finite outer radius, Angstrom.
- `bins`: Bin count. Zero selects about rmax / 0.1.
- `type_i`: Central species, as in seams_rdf().
- `type_j`: Counted species, as in seams_rdf().
- `cn`: Receives the coordination number.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns sys; the call borrows it. cn is caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
double cn = 0;
seams_coordination(sys, 3.5, 0, 1, 1, &cn);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_density_axis

Purpose: Number density along one cell axis.

z is Angstrom. rho is number density, atoms per cubic Angstrom.

```c
seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis, double *z, double *rho, size_t cap, size_t *n_out);
```

Parameters: `sys`, `type`, `bins`, `axis`, `z`, `rho`, `cap`, `n_out`.

- `sys`: Frame handle.
- `type`: Species filter. Zero keeps every atom.
- `bins`: Bin count. Zero selects about one bin per 0.1 Angstrom.
- `axis`: 0, 1, or 2.
- `z`: Bin centers, or NULL to query.
- `rho`: Density per bin, or NULL to query.
- `cap`: Capacity of each non-NULL column.
- `n_out`: Receives the bin count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns the columns.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
double z[50], rho[50];
size_t bins = 0;
seams_density_axis(sys, 0, 50, 2, z, rho, 50, &bins);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_domains

Purpose: Largest domain of a polar or apolar site subset.

```c
seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec, int subset, int *n_subset, int *largest, double *percolation);
```

Parameters: `sys`, `cutoff`, `site_spec`, `subset`, `n_subset`, `largest`, `percolation`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cluster cutoff, Angstrom.
- `site_spec`: Site table text. Must be non-NULL.
- `subset`: 0 is apolar, 1 is polar.
- `n_subset`: Receives the number of sites in the subset, or NULL.
- `largest`: Receives the size of the largest domain, or NULL.
- `percolation`: Receives the largest domain size over the subset size, or NULL. At least one output must be non-NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns sys; the call borrows it. site_spec is read during the call only. n_subset, largest and percolation are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
int largest = 0;
seams_domains(sys, 5.0, "1=cationHead,2=anion,3=tail", 0, NULL, &largest, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint

Purpose: Topological fingerprint of the bonded graph.

```c
seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k, int graph, int hops, int color_types, seams_topo **out);
```

Parameters: `sys`, `cutoff`, `type`, `k`, `graph`, `hops`, `color_types`, `out`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `k`: Neighbour count for the mutual graph. Zero selects 4.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `hops`: Neighbourhood depth. Zero selects 2.
- `color_types`: Non-zero colors vertices by LAMMPS type.
- `out`: Receives a ::seams_topo. On failure *out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller frees *out with seams_fingerprint_destroy().

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_topo *fp = NULL;
seams_fingerprint(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 0, &fp);
seams_fingerprint_destroy(fp);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_atom

Purpose: Local key of one atom. Same buffer rule as seams_fingerprint_key().

```c
seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf, size_t cap, size_t *n_out);
```

Parameters: `fp`, `index`, `buf`, `cap`, `n_out`.

- `fp`: Fingerprint handle.
- `index`: Atom index in the frame.
- `buf`: Caller buffer of cap bytes, or NULL to query.
- `cap`: Capacity of buf, including the terminating NUL.
- `n_out`: Receives the key length without the terminating NUL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE for a bad index or a short buffer.

Ownership: The caller owns buf.

Thread safety: Read-only: may run alongside other const calls on fp, and not alongside a non-const call on fp.

Units: as given for each parameter above.

Example:

```c
char key[256];
size_t len = 0;
seams_fingerprint_atom(fp, 0, key, sizeof key, &len);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_destroy

Purpose: Free a fingerprint.

```c
void seams_fingerprint_destroy(seams_topo *fp);
```

Parameters: `fp`.

- `fp`: Handle from seams_fingerprint(), or NULL. NULL is ignored.

Returns: nothing.

Ownership: Releases the handle and its key strings.

Thread safety: The caller must not use fp after this returns.

Units: as given for each parameter above.

Example:

```c
seams_fingerprint_destroy(fp);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_key

Purpose: Frame key text.

```c
seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap, size_t *n_out);
```

Parameters: `fp`, `buf`, `cap`, `n_out`.

- `fp`: Fingerprint handle.
- `buf`: Caller buffer of cap bytes, or NULL to query.
- `cap`: Capacity of buf, including the terminating NUL.
- `n_out`: Receives the key length without the terminating NUL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when buf is short.

Ownership: The caller owns buf. The copy is NUL-terminated.

Thread safety: Read-only: may run alongside other const calls on fp, and not alongside a non-const call on fp.

Units: as given for each parameter above.

Example:

```c
size_t len = 0;
seams_fingerprint_key(fp, NULL, 0, &len);
char *key = malloc(len + 1);
seams_fingerprint_key(fp, key, len + 1, &len);
free(key);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_fingerprint_rings

Purpose: Primitive-ring census stored on the fingerprint.

```c
seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts, size_t counts_len);
```

Parameters: `fp`, `counts`, `counts_len`.

- `fp`: Fingerprint handle.
- `counts`: Caller buffer. counts[k] receives the number of rings of size k. Sizes past the buffer are dropped.
- `counts_len`: Length of counts. Must be positive.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when counts_len is 0.

Ownership: The caller owns counts.

Thread safety: Read-only: may run alongside other const calls on fp, and not alongside a non-const call on fp.

Units: as given for each parameter above.

Example:

```c
int census[13];
seams_fingerprint_rings(fp, census, 13);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_get_num_threads

Purpose: Team width the next analysis call will request.

```c
int seams_get_num_threads(void);
```

Parameters: none.

Returns: A positive count. After seams_set_num_threads(0) this is the saved default, which is omp_get_max_threads() when OpenMP is linked and 1 otherwise.

Thread safety: Observes the process-global width. Do not call it concurrently with seams_set_num_threads().

Units: as given for each parameter above.

Example:

```c
int width = seams_get_num_threads();
```

The tour call is in `docs/examples/c_api_tour.c`.

Ownership: No allocation.

Errors: None.

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

Purpose: Detail string stored on a handle.

```c
const char *seams_handle_error(const void *handle);
```

Parameters: `handle`.

- `handle`: A ::seams_system, ::seams_neighbours, ::seams_rings, ::seams_cage_set, ::seams_topo, or ::seams_brine, or NULL.

Returns: The handle's string, or the thread-local string when handle is NULL. Never NULL. A call taking the handle non-const clears the string on entry. A call taking it const writes the string only when it fails, so after a successful const call the string still holds the last failure.

Ownership: Owned by the handle. Valid until the next failing call on that handle, the next non-const call on it, or destroy.

Thread safety: Reading a handle requires that no other thread is inside a call on that same handle.

Units: as given for each parameter above.

Example:

```c
if (seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                       labels, n, NULL) != SEAMS_OK) {
  fprintf(stderr, "%s\n", seams_handle_error(sys));
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_hbonds

Purpose: Hydrogen-bond count on one frame.

```c
seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int *n_bonds);
```

Parameters: `sys`, `cutoff`, `oxygen_type`, `hydrogen_type`, `dist_cutoff`, `angle_cutoff`, `n_bonds`.

- `sys`: Frame handle.
- `cutoff`: Positive finite donor-acceptor cutoff, Angstrom.
- `oxygen_type`: Heavy-atom species. Zero keeps every atom that is not hydrogen_type.
- `hydrogen_type`: Positive hydrogen species.
- `dist_cutoff`: Acceptor-H distance, Angstrom.
- `angle_cutoff`: O-O-H angle, degrees.
- `n_bonds`: Receives the number of hydrogen bonds.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns sys; the call borrows it. n_bonds is caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
int n_bonds = 0;
seams_hbonds(sys, 3.5, 1, 2, 2.4, 30.0, &n_bonds);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Name of a CHILL label.

```c
const char *seams_label_name(int label);
```

Parameters: `label`.

- `label`: A ::seams_label value.

Returns: Static string. Unknown codes return "unknown".

Ownership: The caller does not free it.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
printf("%s\n", seams_label_name(SEAMS_LABEL_HEXAGONAL));
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_lammps_bounds

Purpose: Fill ox, oy, oz, lx, ly, lz, xy, xz, yz from LAMMPS box bounds. Declared in include/seams_lammps_box.h.

Parameters: lo; hi; tilt; box9

Returns: Nothing. The nine cell numbers are written to box9.

Ownership: The caller owns lo, hi, tilt, and box9. The function does not allocate.

Errors: None. The caller supplies three-vectors and a nine-number buffer.

Units: lengths in angstrom

Example: `docs/examples/c_api_tour.c`

- `label`: A ::seams_label value.
Thread safety: Pure.
```c
printf("%s\n", seams_label_name(SEAMS_LABEL_HEXAGONAL));
```
The tour call is in `docs/examples/c_api_tour.c`.

### seams_last_error

Purpose: Thread-local detail of the last failing call on this thread.

ABI 1 spelling of seams_last_error_message(). Both return the same string.

```c
const char *seams_last_error(void);
```

Parameters: none.

Returns: Pointer to a thread-local string. Never NULL.

Ownership: Owned by the library. Valid until the next call on this thread that sets or clears the error.

Thread safety: The storage is thread-local.

Units: as given for each parameter above.

Example:

```c
if (seams_system_set_atoms(sys, xyz, n, box9, 9, NULL, NULL) != SEAMS_OK) {
  fprintf(stderr, "%s\n", seams_last_error());
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_last_error_message

Purpose: Thread-local detail of the last failing call on this thread.

The string is empty after a successful call on this thread. Prefer seams_handle_error() when the failing call had a handle: that copy survives a later call on another handle. seams_last_error() is the ABI 1 spelling and returns the same string.

```c
const char *seams_last_error_message(void);
```

Parameters: none.

Returns: Pointer to a thread-local string. Never NULL.

Ownership: Owned by the library. Valid until the next call on this thread that sets or clears the error.

Thread safety: The storage is thread-local.

Units: as given for each parameter above.

Example:

```c
seams_status rc = seams_system_natoms(NULL, NULL);
fprintf(stderr, "%s: %s\n", seams_status_string(rc), seams_last_error_message());
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_monolayer

Purpose: Primitive-ring histogram of a monolayer.

```c
seams_status seams_monolayer(seams_system *sys, double cutoff, int type, int max_depth, double sheet_area, int *counts, size_t counts_len);
```

Parameters: `sys`, `cutoff`, `type`, `max_depth`, `sheet_area`, `counts`, `counts_len`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `max_depth`: Largest ring size, at least 3.
- `sheet_area`: Sheet area, square Angstrom. Must be positive.
- `counts`: counts[k] is the number of rings of size k.
- `counts_len`: Must be at least max_depth + 1.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns counts.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
int counts[7];
seams_monolayer(sys, 3.5, 0, 6, 900.0, counts, 7);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Build a cutoff or mutual neighbour list.

```c
seams_status seams_neighbours_build(seams_system *sys, double cutoff, int type, int graph, seams_neighbours **out);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `out`.

- `sys`: Frame handle that already holds a frame.
- `cutoff`: Positive finite cutoff, Angstrom. The mutual graph keeps the four nearest neighbours inside cutoff + 1.5.
- `type`: Species to keep. Zero and negative keep every atom.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL. SEAMS_GRAPH_SEEDED is rejected here.
- `out`: Receives a new handle. On failure *out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller frees *out with seams_neighbours_destroy(). The list is a copy and does not borrow sys.

Thread safety: Exclusive use of sys. The new handle is private to the caller.

Units: as given for each parameter above.

Example:

```c
seams_neighbours *nb = NULL;
if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK) {
  seams_neighbours_destroy(nb);
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_copy

Purpose: Copy the neighbour CSR into caller buffers.

```c
seams_status seams_neighbours_copy(const seams_neighbours *nb, int *offsets, size_t offsets_cap, size_t *n_offsets, int *cols, size_t cols_cap, size_t *nnz);
```

Parameters: `nb`, `offsets`, `offsets_cap`, `n_offsets`, `cols`, `cols_cap`, `nnz`.

- `nb`: Neighbour handle.
- `offsets`: Caller buffer of offsets_cap ints, or NULL to query.
- `offsets_cap`: Capacity of offsets. Ignored when offsets is NULL.
- `n_offsets`: Receives the required offset count (n + 1).
- `cols`: Caller buffer of cols_cap ints, or NULL to query.
- `cols_cap`: Capacity of cols. Ignored when cols is NULL.
- `nnz`: Receives the required column count.

Returns: SEAMS_OK when every non-NULL buffer was filled, or when both buffers were NULL and the lengths were written. SEAMS_ERR_RANGE when a non-NULL buffer is short. SEAMS_ERR_NULL when nb, n_offsets, or nnz is NULL.

Ownership: The caller owns the buffers. The handle is not modified.

Thread safety: Read-only: may run alongside other const calls on nb, and not alongside a non-const call on nb.

Units: as given for each parameter above.

Example:

```c
size_t n_off = 0, nnz = 0;
seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz);
int *off = malloc(n_off * sizeof *off), *cols = malloc(nnz * sizeof *cols);
seams_neighbours_copy(nb, off, n_off, &n_off, cols, nnz, &nnz);
free(off);
free(cols);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_count

Purpose: Number of neighbour rows, as an int.

ABI 1 query. seams_neighbours_nrows() reports errors.

```c
int seams_neighbours_count(const seams_neighbours *nb);
```

Parameters: `nb`.

- `nb`: Neighbour handle, or NULL.

Returns: The row count, or 0 when nb is NULL.

Ownership: The caller owns nb; the call borrows it. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on nb, and not alongside a non-const call on nb.

Units: as given for each parameter above.

Example:

```c
int rows = seams_neighbours_count(nb);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_csr

Purpose: Borrow the neighbour list as CSR.

Offsets has n + 1 entries. Columns omit the self entry.

```c
seams_status seams_neighbours_csr(const seams_neighbours *nb, const int **offsets, int *n_offsets, const int **cols, int *nnz);
```

Parameters: `nb`, `offsets`, `n_offsets`, `cols`, `nnz`.

- `nb`: Neighbour handle.
- `offsets`: Receives the offset pointer.
- `n_offsets`: Receives the offset count, n + 1.
- `cols`: Receives the column pointer. May be NULL when the list has no pairs; *nnz is then 0.
- `nnz`: Receives the column count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL.

Ownership: The pointers are owned by nb and remain valid until seams_neighbours_destroy(nb). Use seams_neighbours_copy() for caller buffers.

Thread safety: Exclusive use of nb. The returned pointers must not be used concurrently with a call that frees nb.

Units: as given for each parameter above.

Example:

```c
const int *off = NULL, *cols = NULL;
int n_off = 0, nnz = 0;
seams_neighbours_csr(nb, &off, &n_off, &cols, &nnz);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_destroy

Purpose: Free a neighbour list.

```c
void seams_neighbours_destroy(seams_neighbours *nb);
```

Parameters: `nb`.

- `nb`: Handle from seams_neighbours_build() or seams_neighbours_knn(), or NULL. NULL is ignored.

Returns: nothing.

Ownership: Releases the handle and its CSR arrays.

Thread safety: The caller must not use nb after this returns.

Units: as given for each parameter above.

Example:

```c
seams_neighbours_destroy(nb);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_knn

Purpose: k-nearest neighbour list.

```c
seams_status seams_neighbours_knn(seams_system *sys, int k, double cutoff, int type, int mutual, seams_neighbours **out);
```

Parameters: `sys`, `k`, `cutoff`, `type`, `mutual`, `out`.

- `sys`: Frame handle.
- `k`: Positive neighbour count.
- `cutoff`: Nominal positive cutoff, Angstrom. Candidates are taken inside cutoff + 1.5.
- `type`: Species filter. Zero and negative keep every atom.
- `mutual`: Non-zero keeps the mutual graph. Zero keeps the one-way k-nearest graph.
- `out`: Receives a neighbour handle. On failure *out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller frees *out with seams_neighbours_destroy().

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_neighbours *nb = NULL;
seams_neighbours_knn(sys, 4, 3.5, 0, 1, &nb);
seams_neighbours_destroy(nb);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_neighbours_nrows

Purpose: Number of neighbour rows.

```c
seams_status seams_neighbours_nrows(const seams_neighbours *nb, size_t *n);
```

Parameters: `nb`, `n`.

- `nb`: Neighbour handle.
- `n`: Receives the row count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL when nb or n is NULL.

Ownership: The caller owns nb; the call borrows it. n is caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on nb, and not alongside a non-const call on nb.

Units: as given for each parameter above.

Example:

```c
size_t rows = 0;
seams_neighbours_nrows(nb, &rows);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_path_nframes

Purpose: Number of frames in a trajectory path.

```c
seams_status seams_path_nframes(const char *path, size_t *nframes);
```

Parameters: `path`, `nframes`.

- `path`: Filesystem path.
- `nframes`: Receives the count. Must be non-NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_INTERNAL. *nframes is unchanged on failure.

Ownership: path is read during the call only. nframes is caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Safe for distinct paths.

Units: as given for each parameter above.

Example:

```c
size_t nframes = 0;
seams_status rc = seams_path_nframes("ice.lammpstrj", &nframes);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_prisms

Purpose: Perfect and imperfect prism counts by ring size.

Rings come from the mutual four-nearest graph inside cutoff + 1.5.

```c
seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect, int *imperfect, size_t len);
```

Parameters: `sys`, `cutoff`, `type`, `perfect`, `imperfect`, `len`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `perfect`: Caller buffer. perfect[k] receives the perfect k-prisms.
- `imperfect`: Caller buffer. imperfect[k] receives the deformed k-prisms.
- `len`: Length of each buffer. At least 7, so that sizes 3 through 6 are filled.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns both buffers.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
int perfect[7], imperfect[7];
seams_prisms(sys, 3.5, 0, perfect, imperfect, 7);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_q4q6

Purpose: Local q4 and q6 on the cutoff shell, without buffer lengths.

ABI 1 pair. Each non-NULL array has one double per atom. New callers use seams_q4q6_n(), which checks the lengths.

```c
seams_status seams_q4q6(seams_system *sys, double cutoff, int type, double *q4, double *q6);
```

Parameters: `sys`, `cutoff`, `type`, `q4`, `q6`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `q4`: Local q4, or NULL.
- `q6`: Local q6, or NULL. q4 and q6 may not both be NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time.

Ownership: The caller owns q4 and q6.

Thread safety: Exclusive use of sys. Any q call on another handle may run beside it.

Units: as given for each parameter above.

Example:

```c
seams_q4q6(sys, 3.5, 0, q4, q6);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_q4q6_n

Purpose: Local q4 and q6 on the cutoff shell into buffers of known length.

```c
seams_status seams_q4q6_n(seams_system *sys, double cutoff, int type, double *q4, size_t q4_len, double *q6, size_t q6_len);
```

Parameters: `sys`, `cutoff`, `type`, `q4`, `q4_len`, `q6`, `q6_len`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `q4`: Local q4, or NULL.
- `q4_len`: Length of q4. Ignored when q4 is NULL.
- `q6`: Local q6, or NULL. q4 and q6 may not both be NULL.
- `q6_len`: Length of q6. Ignored when q6 is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE for a short buffer or a bad cutoff, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time.

Ownership: The caller owns q4 and q6.

Thread safety: Exclusive use of sys. Any q call on another handle may run beside it.

Units: as given for each parameter above.

Example:

```c
seams_q4q6_n(sys, 3.5, 0, q4, n, q6, n);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_q6bar

Purpose: Lechner-Dellago q6-bar and the ice count above a cut.

The average is the same neighbour mean of q_6m as seams_ql() of order 6 on the cutoff graph. n_ice counts finite values strictly above ice_cut. The nucleation census uses 0.36. mean and median ignore non-finite values. An empty finite set stores 0.

```c
seams_status seams_q6bar(seams_system *sys, double cutoff, int type, double ice_cut, double *q6bar, size_t q6bar_len, double *mean, double *median, int *n_ice);
```

Parameters: `sys`, `cutoff`, `type`, `ice_cut`, `q6bar`, `q6bar_len`, `mean`, `median`, `n_ice`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `ice_cut`: Finite threshold on q6-bar.
- `q6bar`: Per-atom average, or NULL to keep only the summary.
- `q6bar_len`: Length of q6bar. Ignored when q6bar is NULL.
- `mean`: Receives the mean, or NULL.
- `median`: Receives the median, or NULL.
- `n_ice`: Receives the count above ice_cut, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. At least one of q6bar, mean, median, and n_ice must be non-NULL.

Handle access: Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time.

Ownership: The caller owns q6bar.

Thread safety: Exclusive use of sys. Any q call on another handle may run beside it.

Units: as given for each parameter above.

Example:

```c
double mean = 0, median = 0;
int n_ice = 0;
seams_q6bar(sys, 3.5, 0, 0.36, NULL, 0, &mean, &median, &n_ice);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_ql

Purpose: Local Steinhardt q_l, and the Lechner-Dellago average on request.

order_l is 3, 4, 6, 8, or 12. Either array may be NULL. A non-NULL array has one double per atom and its length must be at least the frame count. The average is the neighbour mean of q_lm. The cutoff graph reads the cutoff shell. The mutual graph reads the four nearest neighbours inside cutoff + 1.5.

```c
seams_status seams_ql(seams_system *sys, double cutoff, int type, int order_l, int graph, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len);
```

Parameters: `sys`, `cutoff`, `type`, `order_l`, `graph`, `ql`, `ql_len`, `qlbar`, `qlbar_len`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter. Zero and negative keep every atom.
- `order_l`: Spherical-harmonic degree.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `ql`: Local q_l, or NULL.
- `ql_len`: Length of ql. Ignored when ql is NULL.
- `qlbar`: Averaged q_l, or NULL to skip the average.
- `qlbar_len`: Length of qlbar. Ignored when qlbar is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time.

Ownership: The caller owns both arrays.

Thread safety: Exclusive use of sys. Any q call on another handle may run beside it.

Units: as given for each parameter above.

Example:

```c
seams_ql(sys, 3.5, 0, 6, SEAMS_GRAPH_CUTOFF, q6, n, q6bar, n);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rdf

Purpose: Partial radial distribution.

r is Angstrom. g is dimensionless. volume is cubic Angstrom.

```c
seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *g, int *count, size_t cap, size_t *n_out, double *volume);
```

Parameters: `sys`, `rmax`, `bins`, `type_i`, `type_j`, `r`, `g`, `count`, `cap`, `n_out`, `volume`.

- `sys`: Frame handle.
- `rmax`: Positive finite outer radius, Angstrom.
- `bins`: Bin count. Zero selects about rmax / 0.1.
- `type_i`: First species. Zero and negative keep every atom of a single-species frame, or every atom as one species otherwise.
- `type_j`: Second species. Zero and negative reuse the type_i species.
- `r`: Bin centers, or NULL to query.
- `g`: g(r) per bin, or NULL to query.
- `count`: Pair count per bin, or NULL to query.
- `cap`: Capacity of each non-NULL column.
- `n_out`: Receives the bin count.
- `volume`: Receives the cell volume, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns the columns.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
size_t bins = 0;
seams_rdf(sys, 8.0, 80, 0, 0, NULL, NULL, NULL, 0, &bins, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Name of a ring-symmetry code.

```c
const char *seams_ring_symmetry_name(int symmetry);
```

Parameters: `symmetry`.

- `symmetry`: A ::seams_ring_symmetry value.

Returns: Static string. Unknown codes return "unknown".

Ownership: The caller does not free it.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
printf("%s\n", seams_ring_symmetry_name(SEAMS_RING_DDC));
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

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

Purpose: Number of primitive rings, as an int.

ABI 1 query. seams_rings_nrings() reports errors.

```c
int seams_rings_count(const seams_rings *rings);
```

Parameters: `rings`.

- `rings`: Ring handle, or NULL.

Returns: The ring count, or 0 when rings is NULL.

Ownership: The caller owns rings; the call borrows it. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings.

Units: as given for each parameter above.

Example:

```c
int nr = seams_rings_count(rings);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_rings_destroy

Purpose: Free a ring set.

```c
void seams_rings_destroy(seams_rings *rings);
```

Parameters: `rings`.

- `rings`: Handle from seams_rings_find(), or NULL. NULL is ignored.

Returns: nothing.

Ownership: Releases the handle and every ring it holds.

Thread safety: The caller must not use rings after this returns.

Units: as given for each parameter above.

Example:

```c
seams_rings_destroy(rings);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_rings_find

Purpose: Primitive rings (Franzblau) up to a maximum size.

```c
seams_status seams_rings_find(seams_system *sys, double cutoff, int type, int graph, int max_depth, seams_rings **out);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `max_depth`, `out`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `max_depth`: Largest ring size, from 3 through 12.
- `out`: Receives a new handle. On failure *out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller frees *out with seams_rings_destroy().

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_rings *rings = NULL;
if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK) {
  seams_rings_destroy(rings);
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_fraction

Purpose: Fraction of atoms that sit on a ring.

```c
seams_status seams_rings_fraction(const seams_rings *rings, int size, double *fraction);
```

Parameters: `rings`, `size`, `fraction`.

- `rings`: Ring handle.
- `size`: Ring size to count. Zero counts atoms that sit on any ring.
- `fraction`: Receives a value in [0, 1].

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_NOMEM.

Ownership: The caller owns rings; the call borrows it. fraction is caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings.

Units: as given for each parameter above.

Example:

```c
double f6 = 0;
seams_rings_fraction(rings, 6, &f6);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_histogram

Purpose: Count primitive rings by size, with an int length.

ABI 1 entry. New callers use seams_rings_histogram_n().

```c
seams_status seams_rings_histogram(const seams_rings *rings, int *counts, int counts_len);
```

Parameters: `rings`, `counts`, `counts_len`.

- `rings`: Ring handle.
- `counts`: Caller buffer. counts[k] receives the number of rings of size k. Sizes outside the buffer are ignored.
- `counts_len`: Length of counts. Must be positive.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: The caller owns counts.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings.

Units: as given for each parameter above.

Example:

```c
int hist[13];
seams_rings_histogram(rings, hist, 13);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_histogram_n

Purpose: Count primitive rings by size into a buffer of known length.

```c
seams_status seams_rings_histogram_n(const seams_rings *rings, int *counts, size_t counts_len);
```

Parameters: `rings`, `counts`, `counts_len`.

- `rings`: Ring handle.
- `counts`: Caller buffer. counts[k] receives the number of rings of size k. Sizes outside the buffer are ignored.
- `counts_len`: Length of counts. Must be positive.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when counts_len is 0.

Ownership: The caller owns counts.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings.

Units: as given for each parameter above.

Example:

```c
int hist[13];
seams_rings_histogram_n(rings, hist, 13);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_labels

Purpose: Size and symmetry of each ring, without buffer lengths.

ABI 1 entry. Each non-NULL buffer has length seams_rings_count(rings). New callers use seams_rings_labels_n(), which checks the lengths.

```c
seams_status seams_rings_labels(const seams_rings *rings, int *sizes, int *symmetry);
```

Parameters: `rings`, `sizes`, `symmetry`.

- `rings`: Ring handle.
- `sizes`: One int per ring, or NULL.
- `symmetry`: One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns both arrays.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. The first symmetry request fills a cache on the handle under a lock; later requests only read it.

Units: as given for each parameter above.

Example:

```c
seams_rings_labels(rings, sizes, symmetry);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_labels_n

Purpose: Size and symmetry of each ring into buffers of known length.

```c
seams_status seams_rings_labels_n(const seams_rings *rings, int *sizes, size_t sizes_len, int *symmetry, size_t symmetry_len);
```

Parameters: `rings`, `sizes`, `sizes_len`, `symmetry`, `symmetry_len`.

- `rings`: Ring handle.
- `sizes`: One int per ring, or NULL.
- `sizes_len`: Length of sizes. Ignored when sizes is NULL.
- `symmetry`: One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.
- `symmetry_len`: Length of symmetry. Ignored when symmetry is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when a buffer is shorter than the ring count, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: The caller owns both arrays.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings.

Units: as given for each parameter above.

Example:

```c
size_t nr = 0;
seams_rings_nrings(rings, &nr);
int *sym = malloc(nr * sizeof *sym);
seams_rings_labels_n(rings, NULL, 0, sym, nr);
free(sym);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_nrings

Purpose: Number of primitive rings.

```c
seams_status seams_rings_nrings(const seams_rings *rings, size_t *n);
```

Parameters: `rings`, `n`.

- `rings`: Ring handle.
- `n`: Receives the ring count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL when rings or n is NULL.

Ownership: The caller owns rings; the call borrows it. n is caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings.

Units: as given for each parameter above.

Example:

```c
size_t nr = 0;
seams_rings_nrings(rings, &nr);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_vertices

Purpose: Atom indices of one ring.

```c
seams_status seams_rings_vertices(const seams_rings *rings, int index, int *atoms, size_t cap, size_t *n_out);
```

Parameters: `rings`, `index`, `atoms`, `cap`, `n_out`.

- `rings`: Ring handle.
- `index`: Ring index, from 0 to seams_rings_count(rings) - 1.
- `atoms`: Caller buffer, or NULL to query.
- `cap`: Capacity of atoms. Ignored when atoms is NULL.
- `n_out`: Receives the vertex count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: The caller owns atoms. Indices are positions in the frame.

Thread safety: Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings.

Units: as given for each parameter above.

Example:

```c
int atoms[12];
size_t len = 0;
seams_rings_vertices(rings, 0, atoms, 12, &len);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_running_cn

Purpose: Running coordination number along the partial RDF.

```c
seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *cn, size_t cap, size_t *n_out);
```

Parameters: `sys`, `rmax`, `bins`, `type_i`, `type_j`, `r`, `cn`, `cap`, `n_out`.

- `sys`: Frame handle.
- `rmax`: Positive finite outer radius, Angstrom.
- `bins`: Bin count. Zero selects about rmax / 0.1.
- `type_i`: Central species, as in seams_rdf().
- `type_j`: Counted species, as in seams_rdf().
- `r`: Bin centers, Angstrom, or NULL to query.
- `cn`: Running coordination number, or NULL to query.
- `cap`: Capacity of each non-NULL column.
- `n_out`: Receives the bin count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns the columns.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
double r[80], cn[80];
size_t bins = 0;
seams_running_cn(sys, 8.0, 80, 0, 0, r, cn, 80, &bins);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_set_num_threads

Purpose: Set the OpenMP team used by later analysis calls.

```c
seams_status seams_set_num_threads(int n);
```

Parameters: `n`.

- `n`: Positive width, or 0 to restore the width captured at first use. A negative width is rejected.

Returns: SEAMS_OK, or SEAMS_ERR_RANGE when n is negative, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Thread safety: Process-global. Do not call it concurrently with analysis. The stored width is then read by every handle.

Units: as given for each parameter above.

Example:

```c
seams_set_num_threads(1);
```

The tour call is in `docs/examples/c_api_tour.c`.

Ownership: The caller owns every buffer it passes. This call does not allocate a handle the caller must free.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_signature_cages

Purpose: Closed cages of one ring census. Same search as seams_cages_signature().

```c
seams_status seams_signature_cages(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `spec`, `n_cages`, `n_atoms`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `spec`: Signature text such as "512" or "hc". Must be non-NULL.
- `n_cages`: Receives the cage count.
- `n_atoms`: Receives the vertex count.

Returns: Same codes as seams_cages_signature().

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns sys; the call borrows it. spec is read during the call only. n_cages and n_atoms are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
size_t n_cages = 0, n_atoms = 0;
seams_signature_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "hc", &n_cages, &n_atoms);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Per-atom site kind from a site-table spec.

```c
seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds, size_t len);
```

Parameters: `sys`, `spec`, `kinds`, `len`.

- `sys`: Frame handle.
- `spec`: Site table text such as "1=cationHead,2=anion". Must be non-NULL.
- `kinds`: One ::seams_site_kind per atom.
- `len`: Length of kinds. Must be at least the frame count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns kinds.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_site_table(sys, "1=cationHead,2=anion", kinds, n);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Short name of a status code.

```c
const char *seams_status_string(seams_status status);
```

Parameters: `status`.

- `status`: A ::seams_status value.

Returns: Pointer to a static string. Unknown codes return "unknown status".

Ownership: The caller does not free it.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
fprintf(stderr, "%s\n", seams_status_string(SEAMS_ERR_BOX));
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_steinhardt

Purpose: Local Steinhardt q_l on the cutoff shell, and the average on request.

This is seams_ql() with the cutoff graph.

```c
seams_status seams_steinhardt(seams_system *sys, double cutoff, int type, int order_l, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len);
```

Parameters: `sys`, `cutoff`, `type`, `order_l`, `ql`, `ql_len`, `qlbar`, `qlbar_len`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `order_l`: 3, 4, 6, 8, or 12.
- `ql`: Local q_l, or NULL.
- `ql_len`: Length of ql. Ignored when ql is NULL.
- `qlbar`: Averaged q_l, or NULL to skip the Lechner-Dellago average.
- `qlbar_len`: Length of qlbar. Ignored when qlbar is NULL.

Returns: Same codes as seams_ql().

Handle access: Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time.

Ownership: The caller owns both arrays.

Thread safety: Exclusive use of sys. Any q call on another handle may run beside it.

Units: as given for each parameter above.

Example:

```c
seams_steinhardt(sys, 3.5, 0, 6, q6, n, NULL, 0);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system

Purpose: One frame. The caller creates and destroys it. Declared in include/seams.h.

Parameters: none

Returns: The opaque type.

Ownership: The caller creates the object and frees it with the matching destroy function.

Errors: None. This is a type.

Units: handle

Example: `docs/examples/c_api_tour.c`

### seams_system_box

Purpose: Copy the stored cell out as nine LAMMPS numbers.

```c
seams_status seams_system_box(const seams_system *sys, double *box9, size_t box_len);
```

Parameters: `sys`, `box9`, `box_len`.

- `sys`: Frame handle.
- `box9`: Caller buffer: ox, oy, oz, lx, ly, lz, xy, xz, yz. The values match the cell that was stored, to rounding.
- `box_len`: Must be 9.

Returns: SEAMS_OK; SEAMS_ERR_NULL; SEAMS_ERR_RANGE when box_len != 9; SEAMS_ERR_STATE when no frame is stored.

Ownership: The caller owns box9.

Thread safety: Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys.

Units: as given for each parameter above.

Example:

```c
double cell[9];
seams_system_box(sys, cell, 9);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_count

Purpose: Number of points in the frame, as an int.

ABI 1 query. It cannot report an error; seams_system_natoms() can.

```c
int seams_system_count(const seams_system *sys);
```

Parameters: `sys`.

- `sys`: Frame handle, or NULL.

Returns: The atom count, or 0 when sys is NULL or holds no frame.

Ownership: The caller owns sys; the call borrows it. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys.

Units: as given for each parameter above.

Example:

```c
int n = seams_system_count(sys);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_system_create

Purpose: Allocate an empty frame.

```c
seams_system *seams_system_create(void);
```

Parameters: none.

Returns: A new handle, or NULL on allocation failure. seams_last_error_message() then explains the failure.

Ownership: The caller owns the handle and passes it to seams_system_destroy().

Thread safety: Safe. The new handle is not shared yet.

Units: as given for each parameter above.

Example:

```c
seams_system *sys = seams_system_create();
if (sys == NULL) { return SEAMS_ERR_NOMEM; }
seams_system_destroy(sys);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: NULL. seams_last_error() then reads 'out of memory'.

Example: `docs/examples/c_api_tour.c`

### seams_system_destroy

Purpose: Free a frame.

```c
void seams_system_destroy(seams_system *sys);
```

Parameters: `sys`.

- `sys`: Handle from seams_system_create(), or NULL. NULL is ignored.

Returns: nothing.

Ownership: The handle and the frame it stored are released. Neighbour, ring, cage, and brine handles built from it stay valid: they hold their own copies.

Thread safety: The caller must not use sys after this returns, and must not call it concurrently with any use of sys.

Units: as given for each parameter above.

Example:

```c
seams_system_destroy(sys);
sys = NULL;
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_system_ids

Purpose: Copy the stored atom ids out.

```c
seams_status seams_system_ids(const seams_system *sys, int *ids, size_t len);
```

Parameters: `sys`, `ids`, `len`.

- `sys`: Frame handle.
- `ids`: Caller buffer for n ints, in stored order.
- `len`: Length of ids. Must be at least n.

Returns: SEAMS_OK; SEAMS_ERR_NULL; SEAMS_ERR_STATE when no frame is stored; SEAMS_ERR_RANGE for a short buffer.

Ownership: The caller owns ids.

Thread safety: Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys.

Units: as given for each parameter above.

Example:

```c
seams_system_ids(sys, ids, n);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_natoms

Purpose: Number of points in the frame.

```c
seams_status seams_system_natoms(const seams_system *sys, size_t *n);
```

Parameters: `sys`, `n`.

- `sys`: Frame handle.
- `n`: Receives the atom count.

Returns: SEAMS_OK, SEAMS_ERR_NULL when sys or n is NULL, or SEAMS_ERR_STATE when sys holds no frame. *n is unchanged on failure.

Ownership: The caller owns sys; the call borrows it. n is caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys.

Units: as given for each parameter above.

Example:

```c
size_t n = 0;
if (seams_system_natoms(sys, &n) == SEAMS_OK) {
  labels = malloc(n * sizeof *labels);
  free(labels);
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_positions

Purpose: Copy the stored coordinates out.

```c
seams_status seams_system_positions(const seams_system *sys, double *xyz, size_t len);
```

Parameters: `sys`, `xyz`, `len`.

- `sys`: Frame handle.
- `xyz`: Caller buffer for 3 n doubles, Angstrom, in stored order.
- `len`: Length of xyz in doubles. Must be at least 3 n.

Returns: SEAMS_OK; SEAMS_ERR_NULL; SEAMS_ERR_STATE when no frame is stored; SEAMS_ERR_RANGE for a short buffer.

Ownership: The caller owns xyz.

Thread safety: Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys.

Units: as given for each parameter above.

Example:

```c
size_t n = 0;
seams_system_natoms(sys, &n);
double *xyz = malloc(3 * n * sizeof *xyz);
seams_system_positions(sys, xyz, 3 * n);
free(xyz);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_read

Purpose: Read one frame from a trajectory file.

```c
seams_status seams_system_read(seams_system *sys, const char *path, int frame, int type_filter);
```

Parameters: `sys`, `path`, `frame`, `type_filter`.

- `sys`: Frame handle. On failure the previously stored frame is kept.
- `path`: Filesystem path. CON, LAMMPS dump, and the formats the linked readers accept.
- `frame`: Frame index, 1-based, matching the file order.
- `type_filter`: Negative keeps every atom. Zero guesses the oxygen species. A positive value keeps that species number.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE (frame index or empty frame), SEAMS_ERR_ARG (unreadable file), SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Handle access: Non-const: replaces the stored frame. On failure the old frame stays.

Ownership: The library copies the frame out of the file.

Thread safety: Exclusive use of sys. Concurrent reads of distinct paths on distinct handles are safe.

Units: as given for each parameter above.

Example:

```c
seams_status rc = seams_system_read(sys, "ice.lammpstrj", 1, -1);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_atoms

Purpose: Store points and a general cell. Lengths are size_t.

On failure the previously stored frame, if any, is kept. A reload on the same handle reuses its per-atom storage. When `ids` is NULL or 1..n and `n` matches the previous load, the id table is kept as well, so reloading a trajectory frame by frame costs one pass over the input after the finiteness check. seams_system_set_positions() also skips the type and id pass.

```c
seams_status seams_system_set_atoms(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len, const int *types, const int *ids);
```

Parameters: `sys`, `xyz`, `n`, `box9`, `box_len`, `types`, `ids`.

- `sys`: Frame handle.
- `xyz`: n times 3 Cartesian coordinates, Angstrom. Every value must be finite.
- `n`: Atom count. Zero is rejected. A count that does not fit in the frame's integer width is rejected.
- `box9`: Cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.
- `box_len`: Must be 9.
- `types`: Atom types, or NULL to store type 1 for every atom.
- `ids`: Atom ids, or NULL to store 1..n. Ids must be unique.

Returns: SEAMS_OK; SEAMS_ERR_NULL for a NULL sys, xyz, or box9; SEAMS_ERR_RANGE for n == 0 or box_len != 9; SEAMS_ERR_BOX for a non-finite cell or a non-positive edge; SEAMS_ERR_ARG for a non-finite coordinate or repeated ids; SEAMS_ERR_NOMEM; SEAMS_ERR_INTERNAL.

Handle access: Non-const: replaces the stored frame. On failure the old frame stays.

Ownership: The library copies the arrays. The caller keeps them.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
const double xyz[6] = {0, 0, 0, 1.5, 0, 0};
const double box9[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
seams_status rc = seams_system_set_atoms(sys, xyz, 2, box9, 9, NULL, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_frame

Purpose: Store points and a general LAMMPS cell, with an int count.

ABI 1 entry. New callers use seams_system_set_atoms().

```c
seams_status seams_system_set_frame(seams_system *sys, const double *xyz, int n, const double *box9, const int *types, const int *ids);
```

Parameters: `sys`, `xyz`, `n`, `box9`, `types`, `ids`.

- `sys`: Frame handle.
- `xyz`: n times 3 coordinates, Angstrom.
- `n`: Atom count. A non-positive count is rejected.
- `box9`: Nine cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.
- `types`: Types, or NULL for type 1.
- `ids`: Ids, or NULL for 1..n.

Returns: Same codes as seams_system_set_atoms().

Handle access: Non-const: replaces the stored frame. On failure the old frame stays.

Ownership: The library copies the arrays.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_status rc = seams_system_set_frame(sys, xyz, 2, box9, NULL, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_ortho

Purpose: Store points and an orthorhombic cell.

```c
seams_status seams_system_set_ortho(seams_system *sys, const double *xyz, int n, double ox, double oy, double oz, double lx, double ly, double lz, const int *types, const int *ids);
```

Parameters: `sys`, `xyz`, `n`, `ox`, `oy`, `oz`, `lx`, `ly`, `lz`, `types`, `ids`.

- `sys`: Frame handle.
- `xyz`: n times 3 coordinates, Angstrom.
- `n`: Atom count.
- `ox`: Origin x, Angstrom.
- `oy`: Origin y, Angstrom.
- `oz`: Origin z, Angstrom.
- `lx`: Edge length x, Angstrom.
- `ly`: Edge length y, Angstrom.
- `lz`: Edge length z, Angstrom.
- `types`: Types, or NULL.
- `ids`: Ids, or NULL.

Returns: Same codes as seams_system_set_frame(). The tilt is zero.

Handle access: Non-const: replaces the stored frame. On failure the old frame stays.

Ownership: The library copies the arrays.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_status rc = seams_system_set_ortho(sys, xyz, 2, 0, 0, 0, 10, 10, 10, NULL, NULL);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_set_positions

Purpose: Overwrite the coordinates, and optionally the cell, of the stored frame.

Types and ids stay as stored. Per-atom labels written by an earlier seams_chill() or seams_chill_plus() are reset. On failure the frame is unchanged. By default the call checks every coordinate first and then writes them in place, which reads xyz twice and needs no extra memory. After seams_system_stage_positions(sys, 1) it reads xyz once instead (see there).

```c
seams_status seams_system_set_positions(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len);
```

Parameters: `sys`, `xyz`, `n`, `box9`, `box_len`.

- `sys`: Frame handle with a stored frame.
- `xyz`: n times 3 Cartesian coordinates, Angstrom. Every value must be finite.
- `n`: Atom count. Must equal the stored count.
- `box9`: New cell (ox, oy, oz, lx, ly, lz, xy, xz, yz), or NULL to keep the stored cell.
- `box_len`: 9 when box9 is not NULL; ignored otherwise.

Returns: SEAMS_OK; SEAMS_ERR_NULL for a NULL sys or xyz; SEAMS_ERR_STATE when no frame is stored; SEAMS_ERR_RANGE when n differs from the stored count or box_len != 9; SEAMS_ERR_BOX for a non-finite cell or a non-positive edge; SEAMS_ERR_ARG for a non-finite coordinate; SEAMS_ERR_NOMEM; SEAMS_ERR_INTERNAL.

Handle access: Non-const: replaces the stored coordinates and cell.

Ownership: The library copies xyz and box9. The caller keeps them and may reuse the buffers as soon as the call returns.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_system_set_atoms(sys, frame0, n, box9, 9, types, NULL);
for (size_t f = 1; f < nframes; ++f) {
  seams_system_set_positions(sys, frames[f], n, NULL, 0);
  seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED,
                     labels, n, NULL);
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_stage_positions

Purpose: Trade memory for a single pass in seams_system_set_positions().

With enable nonzero the handle keeps a second set of per-atom records, about 72 bytes per atom. seams_system_set_positions() then checks and writes the coordinates into that set in one pass and swaps it with the stored frame when every value was finite, so a failed call still leaves the frame unchanged. The first reload after each load fills the second set from the frame. With enable 0, the default, the second set is freed. The setting survives loads.

```c
seams_status seams_system_stage_positions(seams_system *sys, int enable);
```

Parameters: `sys`, `enable`.

- `sys`: Frame handle, with or without a stored frame.
- `enable`: Nonzero to keep the second set, 0 to free it.

Returns: SEAMS_OK; SEAMS_ERR_NULL for a NULL sys.

Handle access: Non-const: changes how later reloads use memory. The stored frame is not touched.

Ownership: The library owns the second set and frees it with enable 0 or seams_system_destroy().

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_system_stage_positions(sys, 1);
for (size_t f = 1; f < nframes; ++f) {
  seams_system_set_positions(sys, frames[f], n, NULL, 0);
}
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_symmetry

Purpose: Per-atom ice shell label.

labels[i] is 0 for ice Ih, 1 for ice Ic, 2 for stacking disorder, or -1 when the shell does not match. A non-positive shell_rmsd selects the library ceiling (the scale-free second-shell residual, 0.15).

```c
seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd, int *labels, size_t labels_len);
```

Parameters: `sys`, `shell_rmsd`, `labels`, `labels_len`.

- `sys`: Frame handle.
- `shell_rmsd`: Positive scale-free residual ceiling, or a non-positive value for the library ceiling.
- `labels`: Caller buffer.
- `labels_len`: Must be at least the frame count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_INTERNAL, SEAMS_ERR_NOMEM.

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns labels.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
seams_system_symmetry(sys, 0.0, labels, n);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_system_types

Purpose: Copy the stored atom types out.

```c
seams_status seams_system_types(const seams_system *sys, int *types, size_t len);
```

Parameters: `sys`, `types`, `len`.

- `sys`: Frame handle.
- `types`: Caller buffer for n ints, in stored order.
- `len`: Length of types. Must be at least n.

Returns: SEAMS_OK; SEAMS_ERR_NULL; SEAMS_ERR_STATE when no frame is stored; SEAMS_ERR_RANGE for a short buffer.

Ownership: The caller owns types.

Thread safety: Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys.

Units: as given for each parameter above.

Example:

```c
seams_system_types(sys, types, n);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

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

Purpose: Packed release number.

```c
int seams_version(void);
```

Parameters: none.

Returns: SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
if (seams_version() < ((2 << 16) | (14 << 8))) { return 1; }
```

The tour call is in `docs/examples/c_api_tour.c`.

Ownership: No allocation.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_version_major

Purpose: Major version of the library at run time.

```c
int seams_version_major(void);
```

Parameters: none.

Returns: SEAMS_VERSION_MAJOR of the build.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
int major = seams_version_major();
```

The tour call is in `docs/examples/c_api_tour.c`.

Ownership: No allocation.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_version_minor

Purpose: Minor version of the library at run time.

```c
int seams_version_minor(void);
```

Parameters: none.

Returns: SEAMS_VERSION_MINOR of the build.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
int minor = seams_version_minor();
```

The tour call is in `docs/examples/c_api_tour.c`.

Ownership: No allocation.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_version_patch

Purpose: Patch version of the library at run time.

```c
int seams_version_patch(void);
```

Parameters: none.

Returns: SEAMS_VERSION_PATCH of the build.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
int patch = seams_version_patch();
```

The tour call is in `docs/examples/c_api_tour.c`.

Ownership: No allocation.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_version_string

Purpose: Version text, for example "2.14.0".

```c
const char *seams_version_string(void);
```

Parameters: none.

Returns: Pointer to a static string. Never NULL.

Ownership: The caller does not free it.

Thread safety: Pure.

Units: as given for each parameter above.

Example:

```c
printf("d-SEAMS %s\n", seams_version_string());
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: None.

Example: `docs/examples/c_api_tour.c`

### seams_signature_cages_edge

Purpose: Closed cages of one ring census with an explicit face-edge window.

Same search as seams_signature_cages(). edge NULL or "" keeps the default window 2.65,3.15 in the coordinate units of the file. "none" leaves the cage bound to the face census.

```c
seams_status seams_signature_cages_edge(seams_system *sys, double cutoff, int type, int graph, const char *spec, const char *edge, size_t *n_cages, size_t *n_atoms);
```

Parameters: `sys`, `cutoff`, `type`, `graph`, `spec`, `edge`, `n_cages`, `n_atoms`.

- `sys`: Frame handle.
- `cutoff`: Positive finite cutoff, Angstrom.
- `type`: Species filter.
- `graph`: SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.
- `spec`: Signature text such as "512" or "hc". Must be non-NULL.
- `edge`: Face-edge window, or NULL. "none" turns the window off.
- `n_cages`: Receives the cage count.
- `n_atoms`: Receives the vertex count.

Returns: Same codes as seams_cages_signature().

Handle access: Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle.

Ownership: The caller owns sys; the call borrows it. spec and edge are read during the call only. n_cages and n_atoms are caller storage that the call fills. The library keeps no pointer to any argument.

Thread safety: Exclusive use of sys.

Units: as given for each parameter above.

Example:

```c
size_t n_cages = 0, n_atoms = 0;
seams_signature_cages_edge(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, "512", "none", &n_cages, &n_atoms);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`

### seams_rings_chordless

Purpose: Chordless cycles of one size.

A cycle is chordless when no edge of the graph joins two of its vertices that are not next to each other on the cycle. Every primitive ring is chordless; a chordless cycle need not be primitive. Each cycle is written once as size atom indices: it starts at its lowest index and runs towards the smaller of that atom's two cycle neighbours. Cycles are sorted lexicographically, so the list does not depend on the thread count.

```c
seams_status seams_rings_chordless(const seams_neighbours *nb, int size, int *rings, size_t rings_len, size_t *n_rings, int *per_atom, size_t per_atom_len);
```

Parameters: `nb`, `size`, `rings`, `rings_len`, `n_rings`, `per_atom`, `per_atom_len`.

- `nb`: Graph from seams_neighbours_build().
- `size`: Cycle length, from 3 through 12.
- `rings`: Caller buffer of size times the cycle count ints, or NULL to skip the list.
- `rings_len`: Length of rings in ints. Ignored when rings is NULL.
- `n_rings`: Receives the cycle count. Required when rings is not NULL.
- `per_atom`: Caller buffer of one int per atom, or NULL. Receives the number of cycles through each atom.
- `per_atom_len`: Length of per_atom. Ignored when per_atom is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. rings, n_rings and per_atom may not all be NULL. A short rings buffer returns SEAMS_ERR_RANGE after n_rings and per_atom are written.

Ownership: The caller owns rings and per_atom. Indices are positions in the frame.

Thread safety: Exclusive use of nb.

Units: as given for each parameter above.

Example:

```c
size_t n_six = 0;
seams_rings_chordless(nb, 6, NULL, 0, &n_six, NULL, 0);
```

The tour call is in `docs/examples/c_api_tour.c`.

Errors: SEAMS_ERR_NULL, SEAMS_ERR_RANGE, SEAMS_ERR_BOX, SEAMS_ERR_ARG, SEAMS_ERR_STATE, SEAMS_ERR_NOMEM, or SEAMS_ERR_INTERNAL. seams_last_error() holds the detail, and seams_handle_error(handle) when the call took a handle.

Example: `docs/examples/c_api_tour.c`
<!-- API_SYMBOLS_END -->

The program those lines come from:

```{literalinclude} ../../../docs/examples/c_api_tour.c
```


## Structure still

```{figure} ../user/tutorials/molecules/c/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih colored by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug c` writes the 2400 by 1800 PNG.
```

