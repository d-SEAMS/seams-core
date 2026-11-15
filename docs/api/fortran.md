# Fortran

`seams` is the module a program uses. `seams_c` is the `iso_c_binding`
image of `include/seams.h`. Kinds come from `iso_fortran_env`:
`seams_real` is `real64` and `seams_int` is `int32`, and both are
asserted equal to the C kinds at compile time. Coordinates are
`real(seams_real)` with shape `(3, n)`, x/y/z down the first index,
which is the same memory order as the C array.

Build it with fpm from `bindings/fortran`, or with gfortran. Link
`-lseams`. `FORD` reads `bindings/fortran/ford.md` and writes the
native reference. The book links that HTML at `/fortran/index.html`.
The symbol list below is the same public surface.

```bash
cd bindings/fortran
fpm test --flag "-I$SEAMS_INC" --link-flag "-L$SEAMS_LIB_DIR -Wl,-rpath,$SEAMS_LIB_DIR"
ford ford.md
```

`test_seams.f90` calls `seams_c` and prints `fortran ok`.
`test/test_high.f90` calls `seams` and prints `fortran high ok`.
Both use the 64-atom cubic cell and the same trajectory counts.

## Quickstart

```fortran
program quickstart
  use seams
  implicit none
  integer, parameter :: n = 8
  type(seams_frame) :: frame
  real(seams_real) :: xyz(3, n), lx, ly, lz
  integer(seams_int), allocatable :: labels(:), counts(:)
  ! The eight-site cubic diamond cell, edge 6.35 angstrom.
  lx = 2.75_seams_real * 4 / sqrt(3.0_seams_real)
  ly = lx
  lz = lx
  xyz = lx * reshape(real([0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.5, 0.5, 0.0, &
    0.75, 0.75, 0.25, 0.5, 0.0, 0.5, 0.75, 0.25, 0.75, 0.0, 0.5, 0.5, &
    0.25, 0.75, 0.75], seams_real), [3, n])
  frame = new_frame()
  call frame%set_ortho(xyz, lx, ly, lz)
  call frame%chill_plus(labels, counts)
  call frame%close()
  if (counts(1) /= n) error stop "every site is cubic"
end program quickstart
```

`labels` and `counts` are allocatable. `counts(1)` is cubic ice,
because the C label 0 is the first Fortran element. Clathrate is
`counts(5)`.

## Ownership

A handle is a derived type with a final procedure. Assignment shares
the C object. `:close()` or the last finalizer destroys it. Neighbour
lists, rings, cages, fingerprints, and brine fronts own their own
copies: closing the frame does not free them. Do not use one handle
from two threads.

## Errors

Pass `stat` and `errmsg`, or neither. With neither, a failed call
`error stop`s with the library message. With `stat`, the `seams.h`
code is stored and the program continues. `errmsg` receives the detail
from `seams_handle_error`.

```fortran
program errors
  use seams
  implicit none
  type(seams_frame) :: frame
  real(seams_real) :: xyz(3, 1), lx, ly, lz
  integer(seams_int) :: stat
  character(len=:), allocatable :: err
  xyz = 0
  lx = -1
  ly = 1
  lz = 1
  frame = new_frame()
  call frame%set_ortho(xyz, lx, ly, lz, stat=stat, errmsg=err)
  if (stat /= SEAMS_OK) print *, err
  if (stat == SEAMS_OK) error stop "a negative edge is rejected"
  call frame%close()
end program errors
```

## Cubic ice

`bindings/fortran/example/cubic.f90` is the 64-atom cell, one thread.
This is the output:

```text
version 2.14.0
atoms 64
cubic 64
counts_cubic 64
q6_mean   0.62853936
pairs 256
rings6 128
fraction6   1.00000000
coordination   4.06349206
ql6_first   0.62853936
cages 0 64 0 0
signature512 0 0
rdf_bins 8 2049.208555
hbonds 0
density_bins 4    0.031232
fingerprint 16
contacts 0 0 0
domain 0 0    0.000000
prisms6 0
monolayer6 128
site0 8
symmetry0 1
brine 1    0.000000   12.701706
q6bar_ice 64  0.62853936
```

`cages` is hexagonal, double-diamond, neither, both. All 64 atoms are
double-diamond. `site0` is `SEAMS_SITE_WATER_O` (8). `symmetry0` is
ice Ic (1). `q6bar_ice` counts atoms above 0.36. The cell has no
hydrogens, so `hbonds` is 0. `signature512` is 0 because this cell is
not a clathrate.

## GenIce sI

`seams_set_num_threads` and `seams_get_num_threads` are in the same module.
`bindings/fortran/test_seams.f90` is the program the language-front script
compiles. It calls every `bind(C)` procedure in the module on a 64-atom
cubic cell. `meson test docs_api_fortran` compiles that program with
gfortran and runs it. The run prints `fortran call cubic 64`,
`fortran call q6 0.628539361055`, `fortran call cages 0 64 0 0`,
`fortran call sig 0 0`, `fortran call cn 4.063492063492`,
`fortran call hbond 0`, `fortran call sheet 128`,
`fortran call vert 6 0`, and `fortran call read 128 46`.

The modules are also the fpm package `dseams`. It links an installed
libseams, so tell fpm where that is before `fpm build` or `fpm test`
(`bindings/fortran/test/check_dseams.f90`):

```bash
export FPM_LDFLAGS="$(pkg-config --libs-only-L seams) -Wl,-rpath,$(pkg-config --variable=libdir seams)"
```

fpm takes a git dependency's manifest from the repository root, which
is why `fpm.toml` there points at `bindings/fortran`. A project depends
on it with

```toml
[dependencies]
dseams = { git = "https://github.com/d-SEAMS/seams-core", tag = "v2.14.0" }
```

`bindings/fortran/examples/downstream` is such a project, with a path
dependency.

`seams_c` exports `SEAMS_OK`, `SEAMS_CHILL_PUBLISHED`, `SEAMS_CHILL_OVITO`,
`SEAMS_GRAPH_CUTOFF`, `SEAMS_GRAPH_MUTUAL`, `SEAMS_LABEL_CUBIC`,
`seams_version_string`, `seams_last_error`, `seams_set_num_threads`,
`seams_get_num_threads`, `seams_system_create`, `seams_system_destroy`,
`seams_system_set_ortho`, `seams_system_set_frame`, `seams_chill_plus`,
`seams_q4q6`, `seams_neighbours_build`, `seams_neighbours_destroy`,
`seams_neighbours_csr`, `seams_rings_find`, `seams_rings_destroy`,
`seams_rings_histogram`, and `seams_rings_fraction`.
`seams_c` also binds `seams_handle_error`, `seams_system_count`,
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

`docs/examples/minimal_fortran.f90` reads
`input/traj/genice_sI.lammpstrj`. The program checks 46 atoms, at
least 150 neighbour pairs, clathrate `counts(5) == 46`, 48 five-rings,
54 six-rings, and a q6 mean above 0.5. It prints `fortran example ok`.

```{literalinclude} ../../../docs/examples/minimal_fortran.f90
:language: fortran
```

```{literalinclude} ../../../bindings/fortran/example/cubic.f90
:language: fortran
:lines: 1-40
```

```{literalinclude} ../../../bindings/fortran/src/seams_c.f90
:language: fortran
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### SEAMS_CAGE_BOTH

Purpose: Both cages.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CAGE_CUBIC

Purpose: Double-diamond cage only.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CAGE_HEX

Purpose: Hexagonal cage only.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CAGE_WATER

Purpose: Atom in neither cage.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CHILL_OVITO

Purpose: OVITO CHILL+ rule. Open c_ij edges.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CHILL_PUBLISHED

Purpose: Published CHILL+ rule. Inclusive c_ij edges.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_ARG

Purpose: A caller argument is the wrong shape.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_BOX

Purpose: The cell is missing or has no volume.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_INTERNAL

Purpose: The library hit an internal failure.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_NOMEM

Purpose: An allocation failed.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_NULL

Purpose: A required pointer is null.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_RANGE

Purpose: A number is outside the accepted range.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_STATE

Purpose: The handle has no frame yet.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_GRAPH_CUTOFF

Purpose: Neighbour graph of pairs inside the cutoff.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_GRAPH_MUTUAL

Purpose: Mutual four-nearest neighbour graph.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_GRAPH_SEEDED

Purpose: Cage graph of mutual seeds and union completions.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_CLATHRATE

Purpose: SEAMS_LABEL_CLATHRATE

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_CUBIC

Purpose: Cage graph of mutual seeds and union completions.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless code 0

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_HEXAGONAL

Purpose: SEAMS_LABEL_HEXAGONAL

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_INTERCLATHRATE

Purpose: SEAMS_LABEL_INTERCLATHRATE

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_INTERFACIAL

Purpose: SEAMS_LABEL_INTERFACIAL

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_RECUBIC

Purpose: SEAMS_LABEL_RECUBIC

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_REHEX

Purpose: SEAMS_LABEL_REHEX

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_UNCLASSIFIED

Purpose: SEAMS_LABEL_UNCLASSIFIED

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_WATER

Purpose: SEAMS_LABEL_WATER

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_OK

Purpose: The call succeeded.

Parameters: none

Returns: The integer parameter.

Ownership: No allocation.

Errors: None.

Units: dimensionless status code

Example: `docs/examples/minimal_fortran.f90`

### brine_assign

Purpose: Brine_assign

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### brine_close

Purpose: Brine_close

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### brine_final

Purpose: Brine_final

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### brine_fronts

Purpose: Brine_fronts

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### brine_pockets

Purpose: Brine_pockets

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### brine_profile

Purpose: Brine_profile

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### brine_rejection

Purpose: Brine_rejection

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### cage_assign

Purpose: Cage_assign

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### cage_close

Purpose: Cage_close

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### cage_counts

Purpose: Cage_counts

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### cage_final

Purpose: Cage_final

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### cage_labels

Purpose: Cage_labels

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_assign

Purpose: Frame_assign

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_brine

Purpose: Frame_brine

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_cages

Purpose: Frame_cages

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_chill

Purpose: Frame_chill

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_chill_plus

Purpose: Frame_chill_plus

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_close

Purpose: Frame_close

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_contact_pairs

Purpose: Frame_contact_pairs

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_coordination

Purpose: Coordination number integrated to rmax.

Parameters: optional stat and errmsg; see the procedure

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: rmax in angstrom; coordination is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_count

Purpose: Frame_count

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_density

Purpose: Frame_density

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_domains

Purpose: Frame_domains

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_final

Purpose: Release one owner's reference to a frame.

Parameters: optional stat and errmsg; see the procedure

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless site codes

Example: `docs/examples/minimal_fortran.f90`

### frame_fingerprint

Purpose: Frame_fingerprint

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_hbonds

Purpose: Frame_hbonds

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_knn

Purpose: Frame_knn

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_monolayer

Purpose: Frame_monolayer

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_neighbours

Purpose: Frame_neighbours

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_prisms

Purpose: Frame_prisms

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_q4q6

Purpose: Frame_q4q6

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_q6bar

Purpose: Frame_q6bar

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_ql

Purpose: Frame_ql

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_rdf

Purpose: Frame_rdf

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_read

Purpose: Frame_read

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_rings

Purpose: Frame_rings

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_running_cn

Purpose: Frame_running_cn

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_set_frame

Purpose: Frame_set_frame

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_set_ortho

Purpose: Frame_set_ortho

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_signature

Purpose: Frame_signature

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_site_kinds

Purpose: Frame_site_kinds

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_steinhardt

Purpose: Frame_steinhardt

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### frame_symmetry

Purpose: Frame_symmetry

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### label_name

Purpose: Label_name

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### neighbour_assign

Purpose: Neighbour_assign

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### neighbour_close

Purpose: Neighbour_close

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### neighbour_columns

Purpose: Neighbour_columns

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### neighbour_count

Purpose: Neighbour_count

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### neighbour_final

Purpose: Neighbour_final

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### neighbour_offsets

Purpose: Neighbour_offsets

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### new_frame

Purpose: New_frame

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### num_threads

Purpose: Num_threads

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### path_nframes

Purpose: Path_nframes

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_assign

Purpose: Ring_assign

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_close

Purpose: Ring_close

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_count

Purpose: Ring_count

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_final

Purpose: Ring_final

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_fraction

Purpose: Ring_fraction

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_histogram

Purpose: Ring_histogram

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_sizes

Purpose: Ring_sizes

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_symmetry

Purpose: Ring_symmetry

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### ring_vertices

Purpose: Ring_vertices

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### set_num_threads

Purpose: Set_num_threads

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### status_text

Purpose: Status_text

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### topo_assign

Purpose: Topo_assign

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### topo_atom_key

Purpose: Topo_atom_key

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### topo_close

Purpose: Topo_close

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### topo_final

Purpose: Topo_final

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### topo_key

Purpose: Topo_key

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### topo_rings

Purpose: Topo_rings

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### version_text

Purpose: Version_text

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: as named by the parameters

Example: `docs/examples/minimal_fortran.f90`

### c_chars

Purpose: Copy a Fortran string into a null-terminated C buffer.

Parameters: text

Returns: An allocatable C string, including the trailing null.

Ownership: The caller owns the allocation.

Errors: An allocation failure stops the program.

Units: text

Example: `docs/examples/minimal_fortran.f90`

### domain_size

Purpose: Largest apolar domain for a site-table spec.

Parameters: sys; spec; n_subset; largest; percolation

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: percolation is dimensionless; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### pair_count

Purpose: Mutual nearest unlike contacts for a site-table spec.

Parameters: sys; spec; n_pairs; n_cation; n_anion

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### read_frame

Purpose: Read one trajectory frame into a handle.

Parameters: sys; path; frame; type_filter

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: frame index is 1-based; coordinates in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_cages

Purpose: Hexagonal and double-diamond cage membership.

Parameters: sys; cutoff; typ; graph; k; complete; out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; handle

Example: `docs/examples/minimal_fortran.f90`

### seams_cages_counts

Purpose: Disjoint hexagonal, cubic, water, and both counts.

Parameters: cages; n_hex; n_cubic; n_water; n_both

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_cages_destroy

Purpose: Free a cage-membership handle.

Parameters: cages

Returns: Nothing.

Ownership: Frees the c_ptr the caller owns.

Errors: None.

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_chill

Purpose: Original CHILL labels and nine population counts.

Parameters: sys; cutoff; typ; graph; rule; labels; labels_len; counts

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; codes are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_chill_plus

Purpose: Write one CHILL+ label per atom.

Parameters: sys; cutoff; type; graph; rule; labels; counts

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_contact_pairs

Purpose: Mutual nearest unlike contacts of a site table.

Parameters: sys; spec; n_pairs; n_cation; n_anion

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_coordination

Purpose: Coordination number integrated to rmax.

Parameters: sys; rmax; bins; type_i; type_j; cn

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: rmax in angstrom; coordination is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_density_axis

Purpose: Number density along one cell axis.

Parameters: sys; typ; bins; axis; z; rho; cap; n_out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: density in inverse cubic angstrom; centres in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_domains

Purpose: Largest polar or apolar domain.

Parameters: sys; cutoff; spec; subset; n_subset; largest; percolation

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; percolation is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint

Purpose: Topological fingerprint of the bonded graph.

Parameters: sys; cutoff; typ; k; graph; hops; colour_types; out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; handle

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint_destroy

Purpose: Free a fingerprint handle.

Parameters: fp

Returns: Nothing.

Ownership: Frees the c_ptr the caller owns.

Errors: None.

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint_key

Purpose: Frame fingerprint text.

Parameters: fp; buf; cap; n_out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: text

Example: `docs/examples/minimal_fortran.f90`

### seams_fingerprint_rings

Purpose: Primitive-ring census stored on the fingerprint.

Parameters: fp; counts; counts_len

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_get_num_threads

Purpose: Report the OpenMP width the next call will use.

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless thread count

Example: `docs/examples/minimal_fortran.f90`

### seams_handle_error

Purpose: Detail string stored on a handle.

Parameters: handle

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: text

Example: `docs/examples/minimal_fortran.f90`

### seams_hbonds

Purpose: Hydrogen-bond count.

Parameters: sys; cutoff; oxygen_type; hydrogen_type; dist_cutoff; angle_cutoff; n_bonds

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: distances in angstrom; angle in degrees

Example: `docs/examples/minimal_fortran.f90`

### seams_last_error

Purpose: Return the thread-local error detail.

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: text

Example: `docs/examples/minimal_fortran.f90`

### seams_monolayer

Purpose: Primitive-ring histogram of a monolayer.

Parameters: sys; cutoff; typ; max_depth; sheet_area; counts; counts_len

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; area in square angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_build

Purpose: Build a cutoff or mutual neighbour list.

Parameters: sys; cutoff; type; graph; out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_csr

Purpose: Return the neighbour list as CSR.

Parameters: nb; offsets; n_offsets; cols; nnz

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless indices

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_destroy

Purpose: Free a neighbour-list handle.

Parameters: nb

Returns: Nothing.

Ownership: Frees the c_ptr the caller owns.

Errors: None.

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_neighbours_knn

Purpose: K-nearest neighbour list.

Parameters: sys; k; cutoff; typ; mutual; out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; k is a dimensionless count

Example: `docs/examples/minimal_fortran.f90`

### seams_prisms

Purpose: Perfect and imperfect prism counts by ring size.

Parameters: sys; cutoff; typ; perfect; imperfect; len

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_q4q6

Purpose: Write local q4 and q6.

Parameters: sys; cutoff; type; q4; q6

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_rdf

Purpose: Partial radial distribution and the cell volume.

Parameters: sys; rmax; bins; type_i; type_j; r; g; count; cap; n_out; volume

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: rmax in angstrom; volume in cubic angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_count

Purpose: Number of primitive rings.

Parameters: rings

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless count

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_destroy

Purpose: Free a ring-set handle.

Parameters: rings

Returns: Nothing.

Ownership: Frees the c_ptr the caller owns.

Errors: None.

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_find

Purpose: Find primitive rings up to a maximum size.

Parameters: sys; cutoff; type; graph; max_depth; out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; max_depth is a dimensionless ring size

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_fraction

Purpose: Fraction of atoms on a ring of one size.

Parameters: rings; size; fraction

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless fraction

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_histogram

Purpose: Count primitive rings by size.

Parameters: rings; counts; counts_len

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### seams_rings_vertices

Purpose: Atom indices of one ring.

Parameters: rings; index; atoms; cap; n_out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless indices

Example: `docs/examples/minimal_fortran.f90`

### seams_running_cn

Purpose: Running coordination number along the partial RDF.

Parameters: sys; rmax; bins; type_i; type_j; r; cn; cap; n_out

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: rmax in angstrom; coordination is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless thread count

Example: `docs/examples/minimal_fortran.f90`

### seams_signature_cages

Purpose: Closed cages whose ring census matches a signature.

Parameters: sys; cutoff; typ; graph; spec; n_cages; n_atoms

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_site_table

Purpose: Per-atom site kind from a site-table spec.

Parameters: sys; spec; kinds; len

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless site codes

Example: `docs/examples/minimal_fortran.f90`

### seams_steinhardt

Purpose: Local q_l and the Lechner-Dellago average.

Parameters: sys; cutoff; typ; order_l; ql; ql_len; qlbar; qlbar_len

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: cutoff in angstrom; q_l is dimensionless

Example: `docs/examples/minimal_fortran.f90`

### seams_system_count

Purpose: Number of points in the frame.

Parameters: sys

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless count

Example: `docs/examples/minimal_fortran.f90`

### seams_system_create

Purpose: Allocate an empty frame handle.

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_system_destroy

Purpose: Free a frame handle.

Parameters: sys

Returns: Nothing.

Ownership: Frees the c_ptr the caller owns.

Errors: None.

Units: handle

Example: `docs/examples/minimal_fortran.f90`

### seams_system_read

Purpose: Read one frame from a trajectory path.

Parameters: sys; path; frame; type_filter

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: frame index is 1-based; coordinates in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_system_set_frame

Purpose: Store points and a nine-number cell.

Parameters: sys; xyz; n; box; types; ids

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_system_set_ortho

Purpose: Store points and an orthorhombic cell.

Parameters: sys; xyz; n; ox; oy; oz; lx; ly; lz; types; ids

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: origin and edge lengths in angstrom

Example: `docs/examples/minimal_fortran.f90`

### seams_version_string

Purpose: Return the version as text.

Parameters: none

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: text

Example: `docs/examples/minimal_fortran.f90`

### signature_count

Purpose: Closed-cage count for one signature on the cutoff graph.

Parameters: sys; spec; n_cages; n_atoms

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless counts

Example: `docs/examples/minimal_fortran.f90`

### site_kinds

Purpose: Per-atom site kind for a Fortran spec string.

Parameters: sys; spec; kinds; n

Returns: Integer(c_int) status. SEAMS_OK is zero.

Ownership: A c_ptr from seams_system_create is freed with seams_system_destroy. Allocatable character buffers from c_chars are owned by the caller.

Errors: A non-zero seams_status. seams_last_error and seams_handle_error hold the detail.

Units: dimensionless site codes

Example: `docs/examples/minimal_fortran.f90`

<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/fortran/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug fortran` writes the 2400 by 1800 PNG.
```

