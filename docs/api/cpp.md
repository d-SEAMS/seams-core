# C++

d-SEAMS 2.14.0. C++ callers include `seams.h` and link `libseams`. The
header is C linkage. Coordinates are Angstrom. A cell is nine LAMMPS
numbers: `ox, oy, oz, lx, ly, lz, xy, xz, yz`.

`libseams_lang` is the wider C ABI used by the language fronts
(`seams_frame_*`). It is not required for a program that only needs
`seams.h`.

`docs/examples/minimal_cpp.cpp` reads the
GenIce sI cell with `seams_frame_read`, then calls
`seams_frame_neighbours`, `seams_frame_steinhardt`,
`seams_frame_chill_plus`, and `seams_frame_rings`.
`meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_cpp.cpp
:language: cpp
```

Build:

```bash
meson setup build -Dwith_tests=true
meson compile -C build
```

CI runs the C API test, which is the same binary a C++ caller links.
The executable is `build/src/test_seams_api`.

```bash
meson test -C build seams_c_api --print-errorlogs
meson test -C build seams_lang_fixtures --print-errorlogs
```

`bindings/cxx/test_lang.cpp` is the `seams_lang_fixtures` test. It calls
every function in `seams_lang.h` on the in-tree trajectories.
`meson test docs_api_cpp` compiles and runs the eight-atom program below.
It calls `seams_set_num_threads`, `seams_chill_plus`, and `seams_q4q6`.

```{literalinclude} ../../../docs/examples/ice.cpp
:language: cpp
```

```{literalinclude} ../../../bindings/cxx/test_lang.cpp
:language: cpp
:lines: 33-40,53-75,108-120,141-159,184-208,301-357,449-453,495-502,1058-1075,1119-1192
```

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### seams_frame

Purpose: Opaque frame used by the language bindings.

Parameters: none

Units: handle

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_box

Purpose: Copy the nine LAMMPS cell numbers.

Parameters: frame; box9

Units: lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_cages

Purpose: Score hexagonal and double-diamond cage membership.

Parameters: frame; cutoff; type; k; graph; complete; per_atom; n_hex; n_cubic; n_water; n_both

Units: cutoff in angstrom; coordinates in angstrom; k is a dimensionless neighbour count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill

Purpose: Write one original CHILL label per atom.

Parameters: frame; cutoff; type; graph; labels; counts9

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill_plus

Purpose: Write one CHILL+ label per atom.

Parameters: frame; cutoff; type; graph; labels; counts9

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill_plus_rule

Purpose: Write CHILL+ labels for a named rule.

Parameters: frame; cutoff; type; graph; rule; labels; counts9

Units: cutoff in angstrom; coordinates in angstrom; rule is text

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_chill_rings

Purpose: Run CHILL+ and primitive rings on one neighbour list.

Parameters: frame; cutoff; type; graph; max_ring; skin; labels; counts9; ring_counts; ring_counts_len; times

Units: cutoff in angstrom; coordinates in angstrom; skin in angstrom; phase times in seconds

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_free

Purpose: Free a language-binding frame.

Parameters: frame

Units: handle

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_from_arrays

Purpose: Build a frame from coordinates and a cell.

Parameters: xyz; n; box9; types; ids

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_from_species

Purpose: Build a frame and keep molecule ids.

Parameters: xyz; n; box9; types; ids; mols

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_ids

Purpose: Pointer to the atom ids stored on the frame.

Parameters: frame

Units: dimensionless ids

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_index

Purpose: Return the frame index that was read.

Parameters: frame

Units: dimensionless index, 1-based in the file

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_monolayer

Purpose: Count rings in sheets of a given area.

Parameters: frame; cutoff; type; max_depth; sheet_area; output_dir; counts; counts_len

Units: cutoff in angstrom; coordinates in angstrom; sheet_area in square angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_neighbour_pairs

Purpose: Return the number of directed neighbour pairs.

Parameters: frame

Units: dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_neighbours

Purpose: Build the neighbour list and discard the rows.

Parameters: frame; cutoff; type; graph

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_nop

Purpose: Return the number of atoms in the frame.

Parameters: frame

Units: dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_positions

Purpose: Pointer to the coordinates stored on the frame.

Parameters: frame

Units: angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_prisms

Purpose: Count perfect and imperfect prism stacks.

Parameters: frame; cutoff; type; max_depth; output_dir; rows; n_rows

Units: cutoff in angstrom; coordinates in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_q4_q6

Purpose: Write local and averaged q4 and q6.

Parameters: frame; cutoff; type; q4; q6; q4bar; q6bar

Units: cutoff in angstrom; coordinates in angstrom; order parameters are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_q4_q6_chill

Purpose: Write q4, q6, and CHILL+ from one cutoff list.

Parameters: frame; cutoff; type; q4; q6; q4bar; q6bar; labels; counts9

Units: cutoff in angstrom; coordinates in angstrom; order parameters are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf

Purpose: Radial distribution of one type pair.

Parameters: frame; rmax; nbins; type_i; type_j; r; g; count; volume

Units: rmax and r in angstrom; g is dimensionless; volume in cubic angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf2d

Purpose: Two-dimensional radial distribution in a monolayer.

Parameters: frame; cutoff; binwidth; output_dir; r; g

Units: cutoff, binwidth, and r in angstrom; g is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf2d_nbin

Purpose: Number of bins in the monolayer RDF.

Parameters: cutoff; binwidth

Units: cutoff and binwidth in angstrom; result is a dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rdf_nbins

Purpose: Number of bins in the bulk RDF.

Parameters: rmax; nbins

Units: rmax in angstrom; nbins is a dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_read

Purpose: Read one frame from a trajectory path.

Parameters: path; frame; type_filter

Units: frame index is dimensionless; coordinates in the file are angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_rings

Purpose: Histogram of primitive rings.

Parameters: frame; cutoff; type; graph; max_depth; counts; counts_len

Units: cutoff in angstrom; coordinates in angstrom; counts[k] is a dimensionless ring count

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_select

Purpose: Keep selected oxygens for later classifiers.

Parameters: frame; type; mols; nmol

Units: type and molecule ids are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_set_positions

Purpose: Replace coordinates on an array-built frame.

Parameters: frame; xyz; n; box9

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_signature

Purpose: Count closed cages whose ring census matches a spec.

Parameters: frame; cutoff; type; graph; spec; n_cages; n_atoms

Units: cutoff in angstrom; coordinates in angstrom; n_cages and n_atoms are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_signature_edge

Purpose: Same cage census with an explicit face-edge window.

Parameters: frame; cutoff; type; graph; spec; edge; n_cages; n_atoms

Units: cutoff in angstrom; coordinates in angstrom; edge lengths in angstrom

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_steinhardt

Purpose: Local Steinhardt q_l on the cutoff shell.

Parameters: frame; cutoff; type; order_l; ql; qlbar

Units: cutoff in angstrom; coordinates in angstrom; ql is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_steinhardt_graph

Purpose: Local Steinhardt q_l on a named graph.

Parameters: frame; cutoff; type; order_l; graph; ql; qlbar

Units: cutoff in angstrom; coordinates in angstrom; ql is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_symmetry

Purpose: Per-atom shell label for Ih, Ic, or stacking disorder.

Parameters: frame; shell_rmsd; labels

Units: shell_rmsd is a dimensionless residual; labels are dimensionless codes

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_types

Purpose: Pointer to the atom types stored on the frame.

Parameters: frame

Units: dimensionless types

Example: `docs/examples/minimal_cpp.cpp`

### seams_frame_wl

Purpose: Normalised Steinhardt w_l for degree 4 or 6.

Parameters: frame; cutoff; type; order_l; wl

Units: cutoff in angstrom; coordinates in angstrom; wl is dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_frames_jump_tau90

Purpose: Time for an H-H vector to rotate 90 degrees.

Parameters: first; second; dt; oxygen_type; hydrogen_type; tau

Units: dt and tau in the caller's time unit; types are dimensionless

Example: `docs/examples/minimal_cpp.cpp`

### seams_lang_version

Purpose: Version string of the language-binding library.

Parameters: none

Units: text

Example: `docs/examples/minimal_cpp.cpp`

### seams_last_error

Purpose: Return the thread-local error detail.

Parameters: none

Units: text

Example: `docs/examples/minimal_cpp.cpp`

### seams_path_nframes

Purpose: Count frames in a trajectory path.

Parameters: path

Units: dimensionless count

Example: `docs/examples/minimal_cpp.cpp`

### seams_phase_times

Purpose: Seconds spent in neighbours, CHILL+, and rings.

Parameters: neighbours; chill; rings; reused

Units: times in seconds; reused is 0 or 1

Example: `docs/examples/minimal_cpp.cpp`

### seams_prism_row

Purpose: One row of a prism stack: size, perfect, imperfect.

Parameters: size; perfect; imperfect

Units: dimensionless counts

Example: `docs/examples/minimal_cpp.cpp`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/cpp/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug cpp` writes the 2400 by 1800 PNG.
```

