# C API

`seams.h` is ABI 1 of libseams. Handles are opaque. Coordinates are Angstrom.
A cell is nine numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`. A zero tilt is orthorhombic.
`lx, ly, lz` are edge lengths.

## Version and errors

`seams_abi` returns 1. `seams_version_major`, `seams_version_minor`, `seams_version_patch`
and `seams_version_string` identify the release. The build writes those
macros from the project version (2.14.0).

`seams_status` is `SEAMS_OK`, `SEAMS_ERR_NULL`, `SEAMS_ERR_RANGE`, `SEAMS_ERR_BOX`,
`SEAMS_ERR_ARG`, `SEAMS_ERR_STATE`, `SEAMS_ERR_NOMEM`, or `SEAMS_ERR_INTERNAL`.
`seams_status_string` names a code. `seams_last_error` is the thread-local detail.

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
