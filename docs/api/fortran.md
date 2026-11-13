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
use seams
type(seams_frame) :: frame
real(seams_real) :: xyz(3, n)
frame = new_frame()
call frame%set_ortho(xyz, lx, ly, lz)
call frame%chill_plus(labels, counts)
call frame%close()
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
call frame%set_ortho(xyz, lx, ly, lz, stat=stat, errmsg=err)
if (stat /= SEAMS_OK) print *, err
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

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CAGE_CUBIC

Purpose: Double-diamond cage only.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CAGE_HEX

Purpose: Hexagonal cage only.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CAGE_WATER

Purpose: Atom in neither cage.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CHILL_OVITO

Purpose: OVITO CHILL+ rule. Open c_ij edges.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_CHILL_PUBLISHED

Purpose: Published CHILL+ rule. Inclusive c_ij edges.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_ARG

Purpose: A caller argument is the wrong shape.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_BOX

Purpose: The cell is missing or has no volume.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_INTERNAL

Purpose: The library hit an internal failure.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_NOMEM

Purpose: An allocation failed.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_NULL

Purpose: A required pointer is null.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_RANGE

Purpose: A number is outside the accepted range.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_ERR_STATE

Purpose: The handle has no frame yet.

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

### SEAMS_GRAPH_SEEDED

Purpose: Cage graph of mutual seeds and union completions.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_CLATHRATE

Purpose: CHILL clathrate. The population is counts(5).

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_CUBIC

Purpose: CHILL cubic ice. counts(1) is this population.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_HEXAGONAL

Purpose: CHILL hexagonal ice.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_INTERCLATHRATE

Purpose: CHILL interfacial clathrate.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_INTERFACIAL

Purpose: CHILL interfacial ice.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_RECUBIC

Purpose: CHILL reclassified cubic ice.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_REHEX

Purpose: CHILL reclassified hexagonal ice.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_UNCLASSIFIED

Purpose: CHILL atom that matched no class.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_LABEL_WATER

Purpose: CHILL liquid water.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### SEAMS_OK

Purpose: The call succeeded.

Parameters: none

Units: dimensionless code

Example: `docs/examples/minimal_fortran.f90`

### brine_assign

Purpose: Share a brine front.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### brine_close

Purpose: Release one owner's reference to a brine front.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### brine_final

Purpose: Release one owner's reference to a brine front.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### brine_fronts

Purpose: The two planar crossings of the ice slab.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### brine_pockets

Purpose: Brine pockets and the liquid reservoir.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### brine_profile

Purpose: Concentration profile along the front axis.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### brine_rejection

Purpose: Molal rejection coefficient.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### cage_assign

Purpose: Share a cage set.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### cage_close

Purpose: Release one owner's reference to a cage set.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### cage_counts

Purpose: Disjoint cage population counts.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### cage_final

Purpose: Release one owner's reference to a cage set.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### cage_labels

Purpose: One cage code per atom.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_assign

Purpose: Share a frame. The last finalizer destroys it.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_brine

Purpose: Planar ice front from an ice assignment.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_cages

Purpose: Hexagonal and double-diamond cage membership.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_chill

Purpose: Original CHILL labels and nine population counts.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_chill_plus

Purpose: CHILL+ labels and nine population counts.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_close

Purpose: Release one owner's reference to a frame.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_contact_pairs

Purpose: Mutual nearest unlike contacts of a site table.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_coordination

Purpose: Coordination number integrated to rmax.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_count

Purpose: Number of points in the frame.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_density

Purpose: Number density along one cell axis.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_domains

Purpose: Largest polar or apolar domain.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_final

Purpose: Release one owner's reference to a frame.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_fingerprint

Purpose: Topological fingerprint of the bonded graph.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_hbonds

Purpose: Hydrogen-bond count.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_knn

Purpose: Build a k-nearest neighbour list.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_monolayer

Purpose: Primitive-ring histogram of a monolayer.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_neighbours

Purpose: Build a cutoff or mutual neighbour list.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_prisms

Purpose: Perfect and imperfect prism counts by ring size.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_q4q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_q6bar

Purpose: Lechner-Dellago q6-bar and the ice count above a cut.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_ql

Purpose: Local q_l and, on request, its neighbour average.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_rdf

Purpose: Partial radial distribution.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_read

Purpose: Read one 1-based frame from a trajectory path.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_rings

Purpose: Find primitive rings up to a maximum size.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_running_cn

Purpose: Running coordination number along the partial RDF.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_set_frame

Purpose: Store points and a nine-number LAMMPS cell.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_set_ortho

Purpose: Store points in an orthorhombic cell. xyz has shape (3, n).

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_signature

Purpose: Count closed cages of one ring census.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_site_kinds

Purpose: Per-atom site kind from a site-table spec.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_steinhardt

Purpose: Local q_l on the cutoff shell, with the average on request.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### frame_symmetry

Purpose: Per-atom ice shell label.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### label_name

Purpose: Name of a CHILL label code.

Parameters: see the type-bound procedure

Units: text

Example: `docs/examples/minimal_fortran.f90`

### neighbour_assign

Purpose: Share a neighbour list.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### neighbour_close

Purpose: Release one owner's reference to a neighbour list.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### neighbour_columns

Purpose: CSR columns copied into an allocatable array.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### neighbour_count

Purpose: Number of neighbour rows.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### neighbour_final

Purpose: Release one owner's reference to a neighbour list.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### neighbour_offsets

Purpose: CSR offsets copied into an allocatable array.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### new_frame

Purpose: Allocate an empty frame. The caller owns the result.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### num_threads

Purpose: Width the next analysis will request.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### path_nframes

Purpose: Number of frames in a trajectory path.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_assign

Purpose: Share a ring set.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_close

Purpose: Release one owner's reference to a ring set.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_count

Purpose: Number of primitive rings.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_final

Purpose: Release one owner's reference to a ring set.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_fraction

Purpose: Fraction of atoms on a ring. Size 0 counts every ring.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_histogram

Purpose: Ring counts by size. Index k is rings of size k.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_sizes

Purpose: Size of each ring.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_symmetry

Purpose: Symmetry code of each ring.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### ring_vertices

Purpose: Atom indices of one ring. The index is 1-based.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### set_num_threads

Purpose: Set the OpenMP width. Zero restores the default.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### status_text

Purpose: Name of a status code.

Parameters: see the type-bound procedure

Units: text

Example: `docs/examples/minimal_fortran.f90`

### topo_assign

Purpose: Share a fingerprint.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### topo_atom_key

Purpose: Local key of one atom. The index is 1-based.

Parameters: see the type-bound procedure

Units: text

Example: `docs/examples/minimal_fortran.f90`

### topo_close

Purpose: Release one owner's reference to a fingerprint.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### topo_final

Purpose: Release one owner's reference to a fingerprint.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### topo_key

Purpose: Frame key of a fingerprint.

Parameters: see the type-bound procedure

Units: text

Example: `docs/examples/minimal_fortran.f90`

### topo_rings

Purpose: Primitive-ring census stored on a fingerprint.

Parameters: optional stat and errmsg; see the procedure

Units: angstrom where a length is taken; counts are dimensionless

Example: `docs/examples/minimal_fortran.f90`

### version_text

Purpose: Library version text.

Parameters: see the type-bound procedure

Units: text

Example: `docs/examples/minimal_fortran.f90`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/fortran/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug fortran` writes the 2400 by 1800 PNG.
```

