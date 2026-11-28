<!-- https://docs.dseams.info/capi/rings.md -->

# Rings {#rings}

Primitive rings and their histogram.

The stable C interface provides 120 exported functions through `include/seams.h`.

```c
// Copyright (c) d-SEAMS developers
/* Ice Ih has 256 primitive six-rings at 3.5 angstroms: 64 basal, 192 prismatic. */
#include "seams.h"

#include <stdio.h>
#include <stdlib.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_rings *rings = NULL;
  int hist[13] = {0};
  int sized[13] = {0};
  size_t nrings = 0;
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) != SEAMS_OK ||
      seams_rings_histogram(rings, hist, 13) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  if (seams_rings_nrings(rings, &nrings) != SEAMS_OK || nrings != 256 ||
      seams_rings_histogram_n(rings, sized, 13) != SEAMS_OK || sized[6] != 256) {
    fprintf(stderr, "nrings %zu hist6 %d\n", nrings, sized[6]);
    return 1;
  }
  int *sizes = calloc(nrings, sizeof(int));
  int *symmetry = calloc(nrings, sizeof(int));
  if (sizes == NULL || symmetry == NULL) {
    return 1;
  }
  if (seams_rings_labels_n(rings, sizes, nrings, symmetry, nrings) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error_message());
    return 1;
  }
  int basal = 0;
  int prismatic = 0;
  int hexagon = 0;
  for (size_t i = 0; i < nrings; ++i) {
    if (sizes[i] != 6) {
      fprintf(stderr, "ring %zu size %d\n", i, sizes[i]);
      return 1;
    }
    if (symmetry[i] == SEAMS_RING_HC_BASAL) {
      ++basal;
    } else if (symmetry[i] == SEAMS_RING_HC_PRISMATIC) {
      ++prismatic;
    } else if (symmetry[i] == SEAMS_RING_HEXAGON) {
      ++hexagon;
    } else {
      fprintf(stderr, "ring %zu symmetry %d\n", i, symmetry[i]);
      return 1;
    }
  }
  if (basal != 64 || prismatic != 192 || hexagon != 0) {
    fprintf(stderr, "symmetry %d %d %d\n", basal, prismatic, hexagon);
    return 1;
  }
  printf("rings6 %d\n", hist[6]);
  free(sizes);
  free(symmetry);
  seams_rings_destroy(rings);
  seams_system_destroy(sys);
  return hist[6] == 256 ? 0 : 1;
}
```

## seams_ring_symmetry_name {#seams-ring-symmetry-name}

```c
const char *seams_ring_symmetry_name(int symmetry)
```

Name of a ring-symmetry code.

**`symmetry`**

A ::seams_ring_symmetry value.

Returns: Static string. Unknown codes return “unknown”.

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code printf(“%s\\n”, seams_ring_symmetry_name(SEAMS_RING_DDC)); @endcode

## seams_rings_find {#seams-rings-find}

```c
seams_status seams_rings_find(seams_system *sys, double cutoff, int type, int graph, int max_depth, seams_rings **out)
```

Primitive rings (Franzblau) up to a maximum size.

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom.

**`type`**

Species filter.

**`graph`**

::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

**`max_depth`**

Largest ring size, from 3 through 12.

**`out`**

Receives a new handle. On failure *out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees *out with seams_rings_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_rings *rings = NULL; if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK) { seams_rings_destroy(rings); } @endcode

## seams_rings_from_csr {#seams-rings-from-csr}

```c
seams_status seams_rings_from_csr( const int *offsets, size_t offsets_len, const int *columns, size_t columns_len, int max_depth, const seams_system *geometry, seams_rings **out)
```

Primitive rings from an explicit caller-supplied CSR graph.

**`offsets`**

CSR offsets.

**`offsets_len`**

Number of offsets; at least one.

**`columns`**

CSR columns, or NULL when columns_len is zero.

**`columns_len`**

Number of columns.

**`max_depth`**

Largest ring size, from 3 through 12.

**`geometry`**

Optional frame supplying coordinates and periodic cell.

**`out`**

Receives a new ring handle.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

offsets has n_vertices + 1 entries and starts at zero. columns contains zero-based vertex indices. Directed input is symmetrized by edge union; self edges and duplicates are ignored. The Franzblau search is the same implementation used by seams_rings_find(). geometry may be NULL. When supplied, its atom count must equal n_vertices and winding queries are enabled.

@par Handle access Read-only access to geometry. @par Ownership The caller owns all arrays and frees *out with seams_rings_destroy(). @par Thread safety Read-only calls may share geometry. The returned handle is independent.

## seams_rings_destroy {#seams-rings-destroy}

```c
void seams_rings_destroy(seams_rings *rings)
```

Free a ring set.

**`rings`**

Handle from seams_rings_find(), or NULL. NULL is ignored.

@par Ownership Releases the handle and every ring it holds. @par Thread safety The caller must not use rings after this returns. @par Example @code seams_rings_destroy(rings); @endcode

## seams_rings_count {#seams-rings-count}

```c
int seams_rings_count(const seams_rings *rings)
```

Number of primitive rings, as an int.

**`rings`**

Ring handle, or NULL.

Returns: The ring count, or 0 when rings is NULL.

ABI 1 query. seams_rings_nrings() reports errors.

@par Ownership The caller owns rings; the call borrows it. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int nr = seams_rings_count(rings); @endcode

## seams_rings_histogram {#seams-rings-histogram}

```c
seams_status seams_rings_histogram(const seams_rings *rings, int *counts, int counts_len)
```

Count primitive rings by size, with an int length.

**`rings`**

Ring handle.

**`counts`**

Caller buffer. counts\[k\] receives the number of rings of size k. Sizes outside the buffer are ignored.

**`counts_len`**

Length of counts. Must be positive.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

ABI 1 entry. New callers use seams_rings_histogram_n().

@par Ownership The caller owns counts. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int hist\[13\]; seams_rings_histogram(rings, hist, 13); @endcode

## seams_rings_labels {#seams-rings-labels}

```c
seams_status seams_rings_labels(const seams_rings *rings, int *sizes, int *symmetry)
```

Size and symmetry of each ring, without buffer lengths.

**`rings`**

Ring handle.

**`sizes`**

One int per ring, or NULL.

**`symmetry`**

One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

ABI 1 entry. Each non-NULL buffer has length seams_rings_count(rings). New callers use seams_rings_labels_n(), which checks the lengths.

@par Ownership The caller owns both arrays. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. The first symmetry request fills a cache on the handle under a lock; later requests only read it. @par Example @code seams_rings_labels(rings, sizes, symmetry); @endcode

## seams_rings_vertices {#seams-rings-vertices}

```c
seams_status seams_rings_vertices(const seams_rings *rings, int index, int *atoms, size_t cap, size_t *n_out)
```

Atom indices of one ring.

**`rings`**

Ring handle.

**`index`**

Ring index, from 0 to seams_rings_count(rings) - 1.

**`atoms`**

Caller buffer, or NULL to query.

**`cap`**

Capacity of atoms. Ignored when atoms is NULL.

**`n_out`**

Receives the vertex count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

@par Ownership The caller owns atoms. Indices are positions in the frame. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int atoms\[12\]; size_t len = 0; seams_rings_vertices(rings, 0, atoms, 12, &len); @endcode

## seams_rings_winding {#seams-rings-winding}

```c
seams_status seams_rings_winding( const seams_rings *rings, int axis, int *axis_winding, size_t axis_winding_len, int *contractible, size_t contractible_len)
```

Periodic-axis winding and contractibility of every ring.

**`rings`**

Ring handle.

**`axis`**

Periodic cell axis, 0, 1, or 2.

**`axis_winding`**

Signed winding output, or NULL.

**`axis_winding_len`**

Length of axis_winding.

**`contractible`**

One for a zero three-component winding, or NULL.

**`contractible_len`**

Length of contractible.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE.

axis_winding receives the oriented integer winding along axis (0, 1, or 2): the minimum-image bond steps summed in the order seams_rings_vertices() lists the ring, in units of that cell vector. A ring that crosses the cell once along +x is +1 on axis 0; the reverse traversal is -1. contractible receives 1 exactly when all three winding components are zero. Either output may be NULL, but not both. Geometry must have been supplied to seams_rings_from_csr(), or carried by seams_rings_find(). @par Ownership The caller owns both output buffers. @par Thread safety Read-only: may run alongside other const calls on rings.

## seams_rings_chordless {#seams-rings-chordless}

```c
seams_status seams_rings_chordless(const seams_neighbours *nb, int size, int *rings, size_t rings_len, size_t *n_rings, int *per_atom, size_t per_atom_len)
```

Chordless cycles of one size.

**`nb`**

Graph from seams_neighbours_build().

**`size`**

Cycle length, from 3 through 12.

**`rings`**

Caller buffer of size times the cycle count ints, or NULL to skip the list.

**`rings_len`**

Length of rings in ints. Ignored when rings is NULL.

**`n_rings`**

Receives the cycle count. Required when rings is not NULL.

**`per_atom`**

Caller buffer of one int per atom, or NULL. Receives the number of cycles through each atom.

**`per_atom_len`**

Length of per_atom. Ignored when per_atom is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL. rings, n_rings and per_atom may not all be NULL. A short rings buffer returns ::SEAMS_ERR_RANGE after n_rings and per_atom are written.

A cycle is chordless when no edge of the graph joins two of its vertices that are not next to each other on the cycle. Every primitive ring is chordless; a chordless cycle need not be primitive. Each cycle is written once as size atom indices: it starts at its lowest index and runs towards the smaller of that atom’s two cycle neighbours. Cycles are sorted lexicographically, so the list does not depend on the thread count. @par Ownership The caller owns rings and per_atom. Indices are positions in the frame. @par Thread safety Exclusive use of nb. @par Example @code size_t n_six = 0; seams_rings_chordless(nb, 6, NULL, 0, &n_six, NULL, 0); @endcode

## seams_rings_fraction {#seams-rings-fraction}

```c
seams_status seams_rings_fraction(const seams_rings *rings, int size, double *fraction)
```

Fraction of atoms that sit on a ring.

**`rings`**

Ring handle.

**`size`**

Ring size to count. Zero counts atoms that sit on any ring.

**`fraction`**

Receives a value in \[0, 1\].

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_NOMEM.

@par Ownership The caller owns rings; the call borrows it. fraction is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code double f6 = 0; seams_rings_fraction(rings, 6, &f6); @endcode

## seams_rings_nrings {#seams-rings-nrings}

```c
seams_status seams_rings_nrings(const seams_rings *rings, size_t *n)
```

Number of primitive rings.

**`rings`**

Ring handle.

**`n`**

Receives the ring count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL when rings or n is NULL.

@par Ownership The caller owns rings; the call borrows it. n is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code size_t nr = 0; seams_rings_nrings(rings, &nr); @endcode

## seams_rings_labels_n {#seams-rings-labels-n}

```c
seams_status seams_rings_labels_n(const seams_rings *rings, int *sizes, size_t sizes_len, int *symmetry, size_t symmetry_len)
```

Size and symmetry of each ring into buffers of known length.

**`rings`**

Ring handle.

**`sizes`**

One int per ring, or NULL.

**`sizes_len`**

Length of sizes. Ignored when sizes is NULL.

**`symmetry`**

One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.

**`symmetry_len`**

Length of symmetry. Ignored when symmetry is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when a buffer is shorter than the ring count, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Ownership The caller owns both arrays. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code size_t nr = 0; seams_rings_nrings(rings, &nr); int *sym = malloc(nr * sizeof *sym); seams_rings_labels_n(rings, NULL, 0, sym, nr); free(sym); @endcode

## seams_rings_histogram_n {#seams-rings-histogram-n}

```c
seams_status seams_rings_histogram_n(const seams_rings *rings, int *counts, size_t counts_len)
```

Count primitive rings by size into a buffer of known length.

**`rings`**

Ring handle.

**`counts`**

Caller buffer. counts\[k\] receives the number of rings of size k. Sizes outside the buffer are ignored.

**`counts_len`**

Length of counts. Must be positive.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when counts_len is 0.

@par Ownership The caller owns counts. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int hist\[13\]; seams_rings_histogram_n(rings, hist, 13); @endcode
