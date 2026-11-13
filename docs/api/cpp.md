# C++

`include/seams.hpp` is the header-only C++ API. It lives in namespace
`seams`, includes `seams.h`, and links `libseams`. Meson installs the
header next to `seams.h`. `pkg-config --cflags seams` and
`find_package(seams)` both put that directory on the include path.
`SEAMS_CXX_HEADER` in `seams-config.cmake` is the installed `seams.hpp`.

## Quickstart

```bash
c++ -std=c++20 -Iinclude docs/examples/cpp_api.cpp -Lbuild/src -lseams -o cpp_api
./cpp_api "$PWD"
```

`meson test docs_api_hpp` compiles that program and runs it on this tree.
The eight-atom diamond is cubic. Ice Ih is hexagonal, 128 atoms, with mean
q6-bar 0.562448843702.

```{literalinclude} ../../../docs/examples/cpp_api.cpp
:language: cpp
```

The wider language ABI remains `libseams_lang` (`seams_frame_*`). A program
that only needs `seams.hpp` does not link it. `docs/examples/minimal_cpp.cpp`
is that ABI on the GenIce sI cell. `meson test docs_examples` runs it.

```{literalinclude} ../../../docs/examples/minimal_cpp.cpp
:language: cpp
```

`docs/examples/ice.cpp` is the eight-atom C call. `meson test docs_api_cpp`
runs it.

```{literalinclude} ../../../docs/examples/ice.cpp
:language: cpp
```

## Concepts

`System`, `Neighbours`, `Rings`, `Cages`, `Fingerprint`, and `Brine` are
move-only. The destructor frees the C handle. Copying is deleted.

A status other than success throws `seams::Error`. `what()` is
`seams_handle_error`. `code()` is the status integer. `status()` is the
`Status` enumerator.

Coordinates enter as `std::span<const double>`. Paths and site-table text
enter as `std::string_view`. Results come back as values, `std::vector`,
`std::array`, or `std::string`. `Neighbours::pairs` returns the directed
pairs. `Rings` is a range: `begin` and `end` walk the rings.

`Graph`, `ChillRule`, `Label`, `CageLabel`, `SiteKind`, `RingSymmetry`,
`Axis`, and `DomainSubset` are `enum class`. Analyses with several knobs
take an options aggregate. `ChillOptions`, `SteinhardtOptions`,
`CageOptions`, `DistributionOptions`, `BondOptions`, and the rest default
to the published cutoff of 3.5 angstroms.

Methods that return a value are `[[nodiscard]]`. A method that reads a
frame takes a non-const `System` because the C call clears the handle
error. `count` is const.

A handle requires exclusive use. The Steinhardt scratch inside the library
is process-global, so two frames must not run `steinhardt`, `ql`, `q4q6`,
or `q6bar` at the same time. `set_num_threads` is process-global.

Doxygen comments on `seams.hpp` and `seams.h` feed the same doxyrest tree
as the engine headers (`Doxyfile_seams.cfg`, the rgpot layout). The
generated pages are `api/index`.

## Ice Ih

`input/traj/ice1h.lammpstrj`, frame 1, every atom, cutoff 3.5 angstrom, one
thread. `tests/test_seams_hpp.cpp` (`meson test seams_hpp_match`) requires
these numbers to match the C++ engine.

```
ih n 128
ih chill0 0 128 0 0 0 0 0 0 0
ih q3 0.744029547233
ih q4 0.509334928671
ih q6 0.625855471508
ih q6bar 0.562448843702
ih q8 0.229702115021
ih q12 0.421988011385
ih vert 6 0
ih cages 128 0 0 0
ih knn 512
ih rdf 8 4158.664606579 0
ih rcn 26.204724409449
ih cn 4.031496062992
ih hbond 0
ih rho 4 0.030779111111
ih fp 16 256
ih site 8
ih pairs 0 128 0
ih domain 0 0 0
ih prism 0
ih sheet 256
```

`chill0` is the nine original CHILL populations in label order. q3 through
q12 are means. q6bar is the mean of the Lechner-Dellago average. `vert` is
the length of ring 0 and its first atom. `cages` is hexagonal, cubic,
water, both. `knn` is the mutual four-neighbour pair count. `rdf` is the
bin count, the volume in cubic angstroms, and g of the first bin. `rcn` is
the last running coordination. `cn` is the coordination to 3.5 angstroms.
`fp` is the fingerprint key length and the ring-census sum. `site` is atom
0 under `1=waterO`. `pairs` is contacts, cations, and anions under
`1=cationHead`. `domain` is the apolar subset. `prism` is perfect
six-stacks. `sheet` is six-rings in the monolayer histogram.

sI has 46 atoms and 2 cages of signature `512` covering 40 atoms. The XYZ
file has 256 atoms and coordination 43.325490196078 out to 6 angstroms.
The CON file has 4 atoms and 0 hydrogen bonds. Triclinic sH has 136 atoms
and mean q6 0.625601998379.

`symmetry` writes one shell label per atom. `brine` measures a planar front
from an ice flag per atom. `docs/examples/cpp_api.cpp` prints:

```
ice Ih atoms 128 hexagonal 128 q6bar 0.562448843702 shell_ih 128 pockets 0
```

The test rejects any other hexagonal count or q6-bar mean.

## Errors and ownership

```cpp
try {
  seams::System system;
  system.read("missing.lammpstrj", 1, -1);
} catch (const seams::Error &err) {
  // err.what() is seams_handle_error. err.code() is the status.
}
```

`System` and the result handles are move-only. A move leaves the source
empty. The destructor of an empty handle does nothing. The caller does not
call `seams_*_destroy`.

`bindings/cxx/test_lang.cpp` is the `seams_lang` fixture. It is a separate
ABI from `seams.hpp`.

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

