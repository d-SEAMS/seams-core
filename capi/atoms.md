<!-- https://docs.dseams.info/capi/atoms.md -->

# Atoms in {#atoms-in}

Read a frame into an opaque handle.

The stable C interface provides 120 exported functions through `include/seams.h`.

```c
// Copyright (c) d-SEAMS developers
/* Read Ice Ih and require 128 atoms. */
#include "seams.h"

#include <stdio.h>

int main(void) {
  seams_system *sys = seams_system_create();
  size_t natoms = 0;
  if (sys == NULL) {
    return 1;
  }
  if (seams_system_read(sys, "input/traj/ice1h.lammpstrj", 1, -1) != SEAMS_OK) {
    fprintf(stderr, "%s\n", seams_last_error());
    return 1;
  }
  int n = seams_system_count(sys);
  if (seams_system_natoms(sys, &natoms) != SEAMS_OK || natoms != 128 || n != 128) {
    fprintf(stderr, "natoms %zu count %d\n", natoms, n);
    return 1;
  }
  printf("atoms %d\n", n);
  seams_system_destroy(sys);
  return 0;
}
```

## seams_system_create {#seams-system-create}

```c
seams_system *seams_system_create(void)
```

Allocate an empty frame.

Returns: A new handle, or NULL on allocation failure. seams_last_error_message() then explains the failure.

@par Ownership The caller owns the handle and passes it to seams_system_destroy(). @par Thread safety Safe. The new handle is not shared yet. @par Example @code seams_system *sys = seams_system_create(); if (sys == NULL) { return SEAMS_ERR_NOMEM; } seams_system_destroy(sys); @endcode

## seams_system_destroy {#seams-system-destroy}

```c
void seams_system_destroy(seams_system *sys)
```

Free a frame.

**`sys`**

Handle from seams_system_create(), or NULL. NULL is ignored.

@par Ownership The handle and the frame it stored are released. Neighbour, ring, cage, and brine handles built from it stay valid: they hold their own copies. @par Thread safety The caller must not use sys after this returns, and must not call it concurrently with any use of sys. @par Example @code seams_system_destroy(sys); sys = NULL; @endcode

## seams_system_set_atoms {#seams-system-set-atoms}

```c
seams_status seams_system_set_atoms(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len, const int *types, const int *ids)
```

Store points and a general cell. Lengths are size_t.

**`sys`**

Frame handle.

**`xyz`**

n times 3 Cartesian coordinates, Angstrom. Every value must be finite.

**`n`**

Atom count. Zero is rejected. A count that does not fit in the frame’s integer width is rejected.

**`box9`**

Cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.

**`box_len`**

Must be 9.

**`types`**

Atom types, or NULL to store type 1 for every atom.

**`ids`**

Atom ids, or NULL to store 1..n. Ids must be unique.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys, xyz, or box9; ::SEAMS_ERR_RANGE for n == 0 or box_len != 9; ::SEAMS_ERR_BOX for a non-finite cell or a non-positive edge; ::SEAMS_ERR_ARG for a non-finite coordinate or repeated ids; ::SEAMS_ERR_NOMEM; ::SEAMS_ERR_INTERNAL.

On failure the previously stored frame, if any, is kept.

A reload on the same handle reuses its per-atom storage. When @p ids is NULL or 1..n and @p n matches the previous load, the id table is kept as well, so reloading a trajectory frame by frame costs one pass over the input after the finiteness check. seams_system_set_positions() also skips the type and id pass.

@par Handle access Non-const: replaces the stored frame. On failure the old frame stays. @par Ownership The library copies the arrays. The caller keeps them. @par Thread safety Exclusive use of sys. @par Example @code const double xyz\[6\] = {0, 0, 0, 1.5, 0, 0}; const double box9\[9\] = {0, 0, 0, 10, 10, 10, 0, 0, 0}; seams_status rc = seams_system_set_atoms(sys, xyz, 2, box9, 9, NULL, NULL); @endcode

## seams_system_set_frame {#seams-system-set-frame}

```c
seams_status seams_system_set_frame(seams_system *sys, const double *xyz, int n, const double *box9, const int *types, const int *ids)
```

Store points and a general LAMMPS cell, with an int count.

**`sys`**

Frame handle.

**`xyz`**

n times 3 coordinates, Angstrom.

**`n`**

Atom count. A non-positive count is rejected.

**`box9`**

Nine cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.

**`types`**

Types, or NULL for type 1.

**`ids`**

Ids, or NULL for 1..n.

Returns: Same codes as seams_system_set_atoms().

ABI 1 entry. New callers use seams_system_set_atoms().

@par Handle access Non-const: replaces the stored frame. On failure the old frame stays. @par Ownership The library copies the arrays. @par Thread safety Exclusive use of sys. @par Example @code seams_status rc = seams_system_set_frame(sys, xyz, 2, box9, NULL, NULL); @endcode

## seams_system_set_ortho {#seams-system-set-ortho}

```c
seams_status seams_system_set_ortho(seams_system *sys, const double *xyz, int n, double ox, double oy, double oz, double lx, double ly, double lz, const int *types, const int *ids)
```

Store points and an orthorhombic cell.

**`sys`**

Frame handle.

**`xyz`**

n times 3 coordinates, Angstrom.

**`n`**

Atom count.

**`ox`**

Origin x, Angstrom.

**`oy`**

Origin y, Angstrom.

**`oz`**

Origin z, Angstrom.

**`lx`**

Edge length x, Angstrom.

**`ly`**

Edge length y, Angstrom.

**`lz`**

Edge length z, Angstrom.

**`types`**

Types, or NULL.

**`ids`**

Ids, or NULL.

Returns: Same codes as seams_system_set_frame(). The tilt is zero.

@par Handle access Non-const: replaces the stored frame. On failure the old frame stays. @par Ownership The library copies the arrays. @par Thread safety Exclusive use of sys. @par Example @code seams_status rc = seams_system_set_ortho(sys, xyz, 2, 0, 0, 0, 10, 10, 10, NULL, NULL); @endcode

## seams_system_read {#seams-system-read}

```c
seams_status seams_system_read(seams_system *sys, const char *path, int frame, int type_filter)
```

Read one frame from a trajectory file.

**`sys`**

Frame handle. On failure the previously stored frame is kept.

**`path`**

Filesystem path. CON, LAMMPS dump, and the formats the linked readers accept.

**`frame`**

Frame index, 1-based, matching the file order.

**`type_filter`**

Negative keeps every atom. Zero guesses the oxygen species. A positive value keeps that species number.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE (frame index or empty frame), ::SEAMS_ERR_ARG (unreadable file), ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const: replaces the stored frame. On failure the old frame stays. @par Ownership The library copies the frame out of the file. @par Thread safety Exclusive use of sys. Concurrent reads of distinct paths on distinct handles are safe. @par Example @code seams_status rc = seams_system_read(sys, “ice.lammpstrj”, 1, -1); @endcode

## seams_path_nframes {#seams-path-nframes}

```c
seams_status seams_path_nframes(const char *path, size_t *nframes)
```

Number of frames in a trajectory path.

**`path`**

Filesystem path.

**`nframes`**

Receives the count. Must be non-NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_INTERNAL. *nframes is unchanged on failure.

@par Ownership path is read during the call only. nframes is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Safe for distinct paths. @par Example @code size_t nframes = 0; seams_status rc = seams_path_nframes(“ice.lammpstrj”, &nframes); @endcode

## seams_system_count {#seams-system-count}

```c
int seams_system_count(const seams_system *sys)
```

Number of points in the frame, as an int.

**`sys`**

Frame handle, or NULL.

Returns: The atom count, or 0 when sys is NULL or holds no frame.

ABI 1 query. It cannot report an error; seams_system_natoms() can.

@par Ownership The caller owns sys; the call borrows it. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys. @par Example @code int n = seams_system_count(sys); @endcode

## seams_system_set_hydrogen_frame {#seams-system-set-hydrogen-frame}

```c
seams_status seams_system_set_hydrogen_frame( seams_system *sys, const double *heavy_xyz, size_t n_heavy, const int *heavy_types, const int *heavy_ids, const int *heavy_mols, const double *hydrogen_xyz, size_t n_hydrogen, const int *hydrogen_types, const int *hydrogen_ids, const int *hydrogen_mols, const double *box9, size_t box_len)
```

Store a heavy-atom cloud and a separate hydrogen cloud.

**`sys`**

Frame handle.

**`heavy_xyz`**

Interleaved heavy coordinates, Angstrom.

**`n_heavy`**

Heavy-atom count. Must be positive.

**`heavy_types`**

One species per heavy atom.

**`heavy_ids`**

One id per heavy atom, or NULL for 1..n_heavy.

**`heavy_mols`**

One molID per heavy atom, or NULL to copy the atom id.

**`hydrogen_xyz`**

Interleaved hydrogen coordinates, Angstrom.

**`n_hydrogen`**

Hydrogen count. Zero allows the hydrogen pointers to be NULL.

**`hydrogen_types`**

One species per hydrogen, required when n_hydrogen is positive.

**`hydrogen_ids`**

One id per hydrogen, or NULL to continue after the heavies.

**`hydrogen_mols`**

One molID per hydrogen, or NULL to copy the atom id.

**`box9`**

ox, oy, oz, lx, ly, lz, xy, xz, yz.

**`box_len`**

Length of box9. Must be 9.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

Heavies occupy frame indices \[0, n_heavy). Hydrogens follow, in the order given, so a hydrogen cloud index h is frame index n_heavy + h. A null id pointer numbers that cloud from 1 across the whole frame. A null mol pointer sets each molID to that atom’s id, which pairs no hydrogen with a heavy atom. Pass molecule ids to attach each hydrogen to a heavy atom. A failed call leaves the previous frame in place. @par Handle access The call replaces the frame on success and clears the handle error on entry. @par Ownership The caller owns sys and every buffer. The call copies the values and keeps no pointer. @par Thread safety Exclusive use of sys. @par Example @code const double heavy\[\] = {0, 0, 0, 2.8, 0, 0}; const double hydrogen\[\] = {1, 0, 0}; const int heavy_types\[\] = {8, 8}; const int heavy_mols\[\] = {1, 2}; const int hydrogen_types\[\] = {2}; const int hydrogen_mols\[\] = {1}; const double box9\[\] = {0, 0, 0, 20, 20, 20, 0, 0, 0}; seams_system_set_hydrogen_frame(sys, heavy, 2, heavy_types, NULL, heavy_mols, hydrogen, 1, hydrogen_types, NULL, hydrogen_mols, box9, 9); @endcode

## seams_system_natoms {#seams-system-natoms}

```c
seams_status seams_system_natoms(const seams_system *sys, size_t *n)
```

Number of points in the frame.

**`sys`**

Frame handle.

**`n`**

Receives the atom count.

Returns: ::SEAMS_OK, ::SEAMS_ERR_NULL when sys or n is NULL, or ::SEAMS_ERR_STATE when sys holds no frame. *n is unchanged on failure.

@par Ownership The caller owns sys; the call borrows it. n is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys. @par Example @code size_t n = 0; if (seams_system_natoms(sys, &n) == SEAMS_OK) { labels = malloc(n * sizeof *labels); free(labels); } @endcode

## seams_system_set_positions {#seams-system-set-positions}

```c
seams_status seams_system_set_positions(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len)
```

Overwrite the coordinates, and optionally the cell, of the stored frame.

**`sys`**

Frame handle with a stored frame.

**`xyz`**

n times 3 Cartesian coordinates, Angstrom. Every value must be finite.

**`n`**

Atom count. Must equal the stored count.

**`box9`**

New cell (ox, oy, oz, lx, ly, lz, xy, xz, yz), or NULL to keep the stored cell.

**`box_len`**

9 when box9 is not NULL; ignored otherwise.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys or xyz; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE when n differs from the stored count or box_len != 9; ::SEAMS_ERR_BOX for a non-finite cell or a non-positive edge; ::SEAMS_ERR_ARG for a non-finite coordinate; ::SEAMS_ERR_NOMEM; ::SEAMS_ERR_INTERNAL.

Types and ids stay as stored. Per-atom labels written by an earlier seams_chill() or seams_chill_plus() are reset. On failure the frame is unchanged.

By default the call checks every coordinate first and then writes them in place, which reads xyz twice and needs no extra memory. After seams_system_stage_positions(sys, 1) it reads xyz once instead (see there).

@par Handle access Non-const: replaces the stored coordinates and cell. @par Ownership The library copies xyz and box9. The caller keeps them and may reuse the buffers as soon as the call returns. @par Thread safety Exclusive use of sys. @par Example @code seams_system_set_atoms(sys, frame0, n, box9, 9, types, NULL); for (size_t f = 1; f < nframes; ++f) { seams_system_set_positions(sys, frames\[f\], n, NULL, 0); seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n, NULL); } @endcode

## seams_system_stage_positions {#seams-system-stage-positions}

```c
seams_status seams_system_stage_positions(seams_system *sys, int enable)
```

Trade memory for a single pass in seams_system_set_positions().

**`sys`**

Frame handle, with or without a stored frame.

**`enable`**

Nonzero to keep the second set, 0 to free it.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys.

With enable nonzero the handle keeps a second set of per-atom records, about 72 bytes per atom. seams_system_set_positions() then checks and writes the coordinates into that set in one pass and swaps it with the stored frame when every value was finite, so a failed call still leaves the frame unchanged. The first reload after each load fills the second set from the frame. With enable 0, the default, the second set is freed. The setting survives loads.

@par Handle access Non-const: changes how later reloads use memory. The stored frame is not touched. @par Ownership The library owns the second set and frees it with enable 0 or seams_system_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_system_stage_positions(sys, 1); for (size_t f = 1; f < nframes; ++f) { seams_system_set_positions(sys, frames\[f\], n, NULL, 0); } @endcode

## seams_system_box {#seams-system-box}

```c
seams_status seams_system_box(const seams_system *sys, double *box9, size_t box_len)
```

Copy the stored cell out as nine LAMMPS numbers.

**`sys`**

Frame handle.

**`box9`**

Caller buffer: ox, oy, oz, lx, ly, lz, xy, xz, yz. The values match the cell that was stored, to rounding.

**`box_len`**

Must be 9.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE when box_len != 9; ::SEAMS_ERR_STATE when no frame is stored.

@par Ownership The caller owns box9. @par Thread safety Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys. @par Example @code double cell\[9\]; seams_system_box(sys, cell, 9); @endcode

## seams_system_positions {#seams-system-positions}

```c
seams_status seams_system_positions(const seams_system *sys, double *xyz, size_t len)
```

Copy the stored coordinates out.

**`sys`**

Frame handle.

**`xyz`**

Caller buffer for 3 n doubles, Angstrom, in stored order.

**`len`**

Length of xyz in doubles. Must be at least 3 n.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

@par Ownership The caller owns xyz. @par Thread safety Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys. @par Example @code size_t n = 0; seams_system_natoms(sys, &n); double *xyz = malloc(3 * n * sizeof *xyz); seams_system_positions(sys, xyz, 3 * n); free(xyz); @endcode

## seams_system_types {#seams-system-types}

```c
seams_status seams_system_types(const seams_system *sys, int *types, size_t len)
```

Copy the stored atom types out.

**`sys`**

Frame handle.

**`types`**

Caller buffer for n ints, in stored order.

**`len`**

Length of types. Must be at least n.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

@par Ownership The caller owns types. @par Thread safety Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys. @par Example @code seams_system_types(sys, types, n); @endcode

## seams_system_ids {#seams-system-ids}

```c
seams_status seams_system_ids(const seams_system *sys, int *ids, size_t len)
```

Copy the stored atom ids out.

**`sys`**

Frame handle.

**`ids`**

Caller buffer for n ints, in stored order.

**`len`**

Length of ids. Must be at least n.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

@par Ownership The caller owns ids. @par Thread safety Read-only: may run alongside other const calls on sys, and not alongside a non-const call on sys. @par Example @code seams_system_ids(sys, ids, n); @endcode
