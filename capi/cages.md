<!-- https://docs.dseams.info/capi/cages.md -->

# Cages {#cages}

Hexagonal, cubic, and signature cages.

The stable C interface provides 120 exported functions through `include/seams.h`.

```c
// Copyright (c) d-SEAMS developers
/* Every Ice Ih oxygen is in a hexagonal cage and not a cubic one. */
#include "seams.h"

#include <stdio.h>

int main(void) {
  seams_system *sys = seams_system_create();
  seams_cage_set *cages = NULL;
  int hex = 0;
  int cubic = 0;
  int water = 0;
  int both = 0;
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    return 1;
  }
  if (seams_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 4, 0, &cages) != SEAMS_OK ||
      seams_cages_counts(cages, &hex, &cubic, &water, &both) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  printf("cages %d %d %d %d\n", hex, cubic, water, both);
  seams_cages_destroy(cages);
  seams_system_destroy(sys);
  return hex == 128 && cubic == 0 && water == 0 && both == 0 ? 0 : 1;
}
```

## seams_cages_find {#seams-cages-find}

```c
seams_status seams_cages_find(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)
```

Hexagonal and double-diamond cage membership.

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom. The seeded graph uses cutoff + 1.5 as the candidate shell.

**`type`**

Species filter.

**`graph`**

::SEAMS_GRAPH_CUTOFF, ::SEAMS_GRAPH_MUTUAL, or ::SEAMS_GRAPH_SEEDED.

**`k`**

Neighbour count for the seeded graph. Zero selects 4.

**`complete`**

Non-zero runs ring-adjacent completion on the seeded labels. Ignored for the cutoff and mutual graphs.

**`out`**

Receives a new handle. On failure *out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees *out with seams_cages_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_cage_set *cages = NULL; if (seams_cages_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 0, 0, &cages) == SEAMS_OK) { seams_cages_destroy(cages); } @endcode

## seams_cages_destroy {#seams-cages-destroy}

```c
void seams_cages_destroy(seams_cage_set *cages)
```

Free a cage-membership handle.

**`cages`**

Handle from seams_cages_find() or seams_cages(), or NULL. NULL is ignored.

@par Ownership Releases the handle and its per-atom labels. @par Thread safety The caller must not use cages after this returns. @par Example @code seams_cages_destroy(cages); @endcode

## seams_cages_counts {#seams-cages-counts}

```c
seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex, int *n_cubic, int *n_water, int *n_both)
```

Disjoint cage population counts.

**`cages`**

Cage handle.

**`n_hex`**

Hexagonal only, or NULL.

**`n_cubic`**

Double-diamond only, or NULL.

**`n_water`**

Neither cage, or NULL.

**`n_both`**

Both cages, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL when cages is NULL or every count pointer is NULL.

The four counts sum to the number of atoms the search kept.

@par Ownership The caller owns cages; the call borrows it. n_hex, n_cubic, n_water and n_both are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on cages, and not alongside a non-const call on cages. @par Example @code int hex = 0, cubic = 0; seams_cages_counts(cages, &hex, &cubic, NULL, NULL); @endcode

## seams_cages_labels {#seams-cages-labels}

```c
seams_status seams_cages_labels(const seams_cage_set *cages, int *per_atom, size_t len)
```

Per-atom cage code.

**`cages`**

Cage handle.

**`per_atom`**

Caller buffer of ::seams_cage_label values, one per atom in the frame. Atoms the species filter dropped are ::SEAMS_CAGE_WATER.

**`len`**

Length of per_atom. Must be at least the frame count stored on the handle.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

@par Ownership The caller owns per_atom. @par Thread safety Read-only: may run alongside other const calls on cages, and not alongside a non-const call on cages. @par Example @code seams_cages_labels(cages, per_atom, n); @endcode

## seams_cages_signature {#seams-cages-signature}

```c
seams_status seams_cages_signature(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)
```

Count closed cages whose ring census matches a signature.

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom.

**`type`**

Species filter.

**`graph`**

::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

**`spec`**

Signature text. Must be non-NULL.

**`n_cages`**

Receives the cage count.

**`n_atoms`**

Receives the vertex count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

spec is a census such as “512”, “51262”, or “51264”, or the names “hc” and “ddc”. n_cages is the number of closed cages. n_atoms is the number of distinct vertices in those cages.

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys; the call borrows it. spec is read during the call only. n_cages and n_atoms are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Exclusive use of sys. @par Example @code size_t n_cages = 0, n_atoms = 0; seams_cages_signature(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, “512”, &n_cages, &n_atoms); @endcode

## seams_system_symmetry {#seams-system-symmetry}

```c
seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd, int *labels, size_t labels_len)
```

Per-atom ice shell label.

**`sys`**

Frame handle.

**`shell_rmsd`**

Positive scale-free residual ceiling, or a non-positive value for the library ceiling.

**`labels`**

Caller buffer.

**`labels_len`**

Must be at least the frame count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_INTERNAL, ::SEAMS_ERR_NOMEM.

labels\[i\] is 0 for ice Ih, 1 for ice Ic, 2 for stacking disorder, or -1 when the shell does not match. A non-positive shell_rmsd selects the library ceiling (the scale-free second-shell residual, 0.15).

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns labels. @par Thread safety Exclusive use of sys. @par Example @code seams_system_symmetry(sys, 0.0, labels, n); @endcode

## seams_cages {#seams-cages}

```c
seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)
```

Hexagonal and double-diamond cages. Same search as seams_cages_find().

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom.

**`type`**

Species filter.

**`graph`**

::SEAMS_GRAPH_CUTOFF, ::SEAMS_GRAPH_MUTUAL, or ::SEAMS_GRAPH_SEEDED.

**`k`**

Neighbour count for the seeded graph. Zero selects 4.

**`complete`**

Non-zero runs ring-adjacent completion on seeded labels.

**`out`**

Receives a ::seams_cage_set. On failure *out is NULL.

Returns: Same codes as seams_cages_find().

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees *out with seams_cages_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_cage_set *cages = NULL; seams_cages(sys, 3.5, 0, SEAMS_GRAPH_SEEDED, 4, 1, &cages); seams_cages_destroy(cages); @endcode

## seams_signature_cages {#seams-signature-cages}

```c
seams_status seams_signature_cages(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)
```

Closed cages of one ring census. Same search as seams_cages_signature().

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom.

**`type`**

Species filter.

**`graph`**

::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

**`spec`**

Signature text such as “512” or “hc”. Must be non-NULL.

**`n_cages`**

Receives the cage count.

**`n_atoms`**

Receives the vertex count.

Returns: Same codes as seams_cages_signature().

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys; the call borrows it. spec is read during the call only. n_cages and n_atoms are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Exclusive use of sys. @par Example @code size_t n_cages = 0, n_atoms = 0; seams_signature_cages(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, “hc”, &n_cages, &n_atoms); @endcode

## seams_signature_cages_edge {#seams-signature-cages-edge}

```c
seams_status seams_signature_cages_edge(seams_system *sys, double cutoff, int type, int graph, const char *spec, const char *edge, size_t *n_cages, size_t *n_atoms)
```

Closed cages of one ring census with an explicit face-edge window.

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom.

**`type`**

Species filter.

**`graph`**

::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

**`spec`**

Signature text such as “512” or “hc”. Must be non-NULL.

**`edge`**

Face-edge window, or NULL. “none” turns the window off.

**`n_cages`**

Receives the cage count.

**`n_atoms`**

Receives the vertex count.

Returns: Same codes as seams_cages_signature().

Same search as seams_signature_cages(). edge NULL or “” keeps the default window 2.65,3.15 in the coordinate units of the file. “none” leaves the cage bound to the face census.

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys; the call borrows it. spec and edge are read during the call only. n_cages and n_atoms are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Exclusive use of sys. @par Example @code size_t n_cages = 0, n_atoms = 0; seams_signature_cages_edge(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, “512”, “none”, &n_cages, &n_atoms); @endcode

## seams_prisms {#seams-prisms}

```c
seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect, int *imperfect, size_t len)
```

Perfect and imperfect prism counts by ring size.

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom.

**`type`**

Species filter.

**`perfect`**

Caller buffer. perfect\[k\] receives the perfect k-prisms.

**`imperfect`**

Caller buffer. imperfect\[k\] receives the deformed k-prisms.

**`len`**

Length of each buffer. At least 7, so that sizes 3 through 6 are filled.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

Rings come from the mutual four-nearest graph inside cutoff + 1.5.

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns both buffers. @par Thread safety Exclusive use of sys. @par Example @code int perfect\[7\], imperfect\[7\]; seams_prisms(sys, 3.5, 0, perfect, imperfect, 7); @endcode

## seams_monolayer {#seams-monolayer}

```c
seams_status seams_monolayer(seams_system *sys, double cutoff, int type, int max_depth, double sheet_area, int *counts, size_t counts_len)
```

Primitive-ring histogram of a monolayer.

**`sys`**

Frame handle.

**`cutoff`**

Positive finite cutoff, Angstrom.

**`type`**

Species filter.

**`max_depth`**

Largest ring size, at least 3.

**`sheet_area`**

Sheet area, square Angstrom. Must be positive.

**`counts`**

counts\[k\] is the number of rings of size k.

**`counts_len`**

Must be at least max_depth + 1.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle’s error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns counts. @par Thread safety Exclusive use of sys. @par Example @code int counts\[7\]; seams_monolayer(sys, 3.5, 0, 6, 900.0, counts, 7); @endcode
