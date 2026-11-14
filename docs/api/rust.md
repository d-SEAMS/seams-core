# Rust

`seams-sys` is the raw FFI crate. `seams` is the safe wrapper. Build both with
`SEAMS_LIB_DIR` pointing at the directory that contains `libseams.so`.

On crates.io the packages are `dseams-sys` and `dseams`, because `seams`
there is an unrelated project. The library names stay `seams_sys` and
`seams`, so the code below is unchanged:

```toml
[dependencies]
dseams = "2.14.0"
```

Without `SEAMS_LIB_DIR`, `dseams-sys` compiles a bundled copy of the C API
and Eigen with the `cc` crate (x86_64, a C++20 compiler with OpenMP). That
copy reads LAMMPS dumps and XYZ. The `system` feature links an installed
libseams through pkg-config instead, with every reader it was built with.
`python3 bindings/rust/seams-sys/vendor.py` fills `seams-sys/vendor/`
before `cargo package`.

```rust
fn main() -> Result<(), seams::Error> {
    let a = 2.75 * 4.0 / 3.0_f64.sqrt();
    let sites = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
        (0.25, 0.25, 0.25),
        (0.75, 0.75, 0.25),
        (0.75, 0.25, 0.75),
        (0.25, 0.75, 0.75),
    ];
    let mut xyz = Vec::new();
    for (x, y, z) in sites {
        xyz.extend([x * a, y * a, z * a]);
    }
    let edge = 3.0 * a;
    let mut sys = seams::System::new()?;
    sys.set_ortho(&xyz, edge, edge, edge)?;
    let (labels, counts) = sys.chill_plus(3.5, seams_sys::SEAMS_GRAPH_CUTOFF, seams_sys::SEAMS_CHILL_PUBLISHED)?;
    let (q4, q6) = sys.q4q6(3.5)?;
    let degree = sys.neighbour_degree(3.5)?;
    let (six, _fraction) = sys.six_ring_fraction(3.5)?;
    assert_eq!(labels.len(), 8);
    assert_eq!(counts.iter().sum::<i32>(), 8);
    assert_eq!(q4.len(), 8);
    assert_eq!(q6.len(), 8);
    assert_eq!(degree.len(), 8);
    assert!(six >= 0);
    Ok(())
}
```

`set_frame` takes the nine-number cell, including a non-zero tilt.
`seams::set_num_threads`, `seams::num_threads`, and `seams::version` match the C API.
`seams::Error` carries the status `code` and the `message` from `seams_handle_error`.
`System::read` loads XYZ, CON, chemfiles, and triclinic dumps.
`chill`, `steinhardt` (including the Lechner-Dellago average), `cages`,
`signature_cages`, `knn`, `rdf`, `running_cn`, `coordination`, `hbonds`,
`density_axis`, `fingerprint`, `contact_pairs`, `domains`, `prisms`,
`monolayer`, `rings` with `vertices`, and `site_table` are the remaining
safe wrappers. `Neighbours`, `Rings`, `Cages`, and `Fingerprint` free their
handles on drop. `bindings/rust/seams/tests/cpp_match.rs` checks those
values against the C++ path on ice Ih, sI, the XYZ file, the CON file,
and the triclinic sH dump. On ice Ih it prints `rustcall hexagonal 128`,
`rustcall q6bar 0.562448843702`, `rustcall vert 6 0`,
`rustcall cages 128 0 0 0`, `rustcall knn 512`, and `rustcall sheet 256`.
The sI signature prints `rustcall sig 2 40`.
`System::count` is the number of points in the frame.
`bindings/rust/seams/tests/ice.rs` is the `cargo test` that calls every
safe wrapper. `meson test docs_api_rust` runs that test.
`docs/examples/minimal_rust.rs` reads the
GenIce sI cell and calls `System::set_ortho`, `count`, `neighbour_degree`,
`q4q6`, `chill_plus`, and `six_ring_fraction`. `meson test docs_examples`
runs it.

```{literalinclude} ../../../docs/examples/minimal_rust.rs
:language: rust
```

`examples/ice.rs` prints the cubic count.

```{literalinclude} ../../../bindings/rust/seams/tests/ice.rs
:language: rust
:lines: 25-57
```

`bindings/rust/seams-sys` publishes each name from `include/seams.h`. The one-line summary, parameters, ownership, and error codes for those names are the C symbol reference.

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### Cages

Purpose: Owned hexagonal and double-diamond membership.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### Error

Purpose: Status code and message from a failed call.

Parameters: code; message

Returns: The status code and the detail string.

Ownership: The caller owns the String.

Errors: This is the error type. code is the seams_status value and message is seams_handle_error or seams_last_error.

Units: code is dimensionless; message is text

Example: `docs/examples/minimal_rust.rs`

### Fingerprint

Purpose: Owned topological fingerprint.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### Neighbours

Purpose: Owned neighbour-list handle.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### Rings

Purpose: Owned primitive-ring handle.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### System

Purpose: Owned frame handle.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### cages

Purpose: Hexagonal and double-diamond cage membership.

Parameters: cutoff; graph; k; complete

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; codes are dimensionless

Example: `docs/examples/minimal_rust.rs`

### chill

Purpose: Original CHILL labels and nine population counts.

Parameters: cutoff; type_filter; graph; rule

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; codes are dimensionless

Example: `docs/examples/minimal_rust.rs`

### chill_plus

Purpose: CHILL+ labels and nine population counts.

Parameters: cutoff; graph; rule

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; graph and rule are dimensionless codes

Example: `docs/examples/minimal_rust.rs`

### contact_pairs

Purpose: Mutual nearest unlike contacts of a site table.

Parameters: spec

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: dimensionless counts

Example: `docs/examples/minimal_rust.rs`

### coordination

Purpose: Coordination number integrated to rmax.

Parameters: rmax; bins

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: rmax in angstrom; coordination is dimensionless

Example: `docs/examples/minimal_rust.rs`

### count

Purpose: Number of points in the system.

Parameters: self

Returns: The value named in the purpose.

Ownership: A returned String is owned by the caller. No C handle is allocated.

Errors: None for the integer queries. version copies library text.

Units: dimensionless count

Example: `docs/examples/minimal_rust.rs`

### counts

Purpose: Disjoint hexagonal, cubic, water, and both counts.

Parameters: self

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: dimensionless counts

Example: `docs/examples/minimal_rust.rs`

### density_axis

Purpose: Number density along one cell axis.

Parameters: type_filter; bins; axis

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: density in inverse cubic angstrom; centres in angstrom

Example: `docs/examples/minimal_rust.rs`

### domains

Purpose: Largest polar or apolar domain.

Parameters: cutoff; spec; subset

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; percolation is dimensionless

Example: `docs/examples/minimal_rust.rs`

### fingerprint

Purpose: Topological fingerprint of the bonded graph.

Parameters: cutoff; hops

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; hops is a dimensionless depth

Example: `docs/examples/minimal_rust.rs`

### hbonds

Purpose: Hydrogen-bond count.

Parameters: cutoff; oxygen_type; hydrogen_type; dist; angle

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: distances in angstrom; angle in degrees; count is dimensionless

Example: `docs/examples/minimal_rust.rs`

### key

Purpose: Frame fingerprint text.

Parameters: self

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: text

Example: `docs/examples/minimal_rust.rs`

### knn

Purpose: K-nearest neighbour list.

Parameters: k; cutoff; mutual

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; k is a dimensionless count

Example: `docs/examples/minimal_rust.rs`

### monolayer

Purpose: Primitive-ring histogram of a monolayer.

Parameters: cutoff; max_depth; sheet_area

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; area in square angstrom

Example: `docs/examples/minimal_rust.rs`

### neighbour_degree

Purpose: Degree of each atom on the cutoff graph.

Parameters: cutoff

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; degrees are dimensionless

Example: `docs/examples/minimal_rust.rs`

### new

Purpose: Allocate an empty System.

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### nnz

Purpose: Number of neighbour pairs, omitting self.

Parameters: self

Returns: The value named in the purpose.

Ownership: A returned String is owned by the caller. No C handle is allocated.

Errors: None for the integer queries. version copies library text.

Units: dimensionless count

Example: `docs/examples/minimal_rust.rs`

### num_threads

Purpose: OpenMP width the next call will use.

Parameters: none

Returns: The value named in the purpose.

Ownership: A returned String is owned by the caller. No C handle is allocated.

Errors: None for the integer queries. version copies library text.

Units: dimensionless thread count

Example: `docs/examples/minimal_rust.rs`

### prisms

Purpose: Perfect and imperfect prism counts by ring size.

Parameters: cutoff

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_rust.rs`

### q4q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: cutoff

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_rust.rs`

### rdf

Purpose: Partial radial distribution and the cell volume.

Parameters: rmax; bins

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: rmax in angstrom; g is dimensionless; volume in cubic angstrom

Example: `docs/examples/minimal_rust.rs`

### read

Purpose: Load one frame from a trajectory path.

Parameters: path; frame; type_filter

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: frame index is 1-based; coordinates in angstrom

Example: `docs/examples/minimal_rust.rs`

### rings

Purpose: Primitive-ring census stored on a fingerprint, or a new ring handle.

Parameters: cutoff; max_depth

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; sizes are dimensionless

Example: `docs/examples/minimal_rust.rs`

### running_cn

Purpose: Running coordination number along the partial RDF.

Parameters: rmax; bins

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: rmax in angstrom; coordination is dimensionless

Example: `docs/examples/minimal_rust.rs`

### set_frame

Purpose: Store points and a nine-number cell.

Parameters: xyz; box9

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_rust.rs`

### set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: dimensionless thread count

Example: `docs/examples/minimal_rust.rs`

### set_ortho

Purpose: Store points in an orthorhombic cell at the origin.

Parameters: xyz; lx; ly; lz

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: coordinates and edge lengths in angstrom

Example: `docs/examples/minimal_rust.rs`

### signature_cages

Purpose: Closed cages whose ring census matches a signature.

Parameters: cutoff; graph; spec

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_rust.rs`

### site_table

Purpose: Per-atom site kind from a site-table spec.

Parameters: spec

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: dimensionless site codes

Example: `docs/examples/minimal_rust.rs`

### six_ring_fraction

Purpose: Fraction of atoms on a primitive six-ring.

Parameters: cutoff

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; fraction is dimensionless

Example: `docs/examples/minimal_rust.rs`

### steinhardt

Purpose: Local q_l and, on request, the Lechner-Dellago average.

Parameters: cutoff; order; average

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; q_l is dimensionless

Example: `docs/examples/minimal_rust.rs`

### version

Purpose: Library version text.

Parameters: none

Returns: The value named in the purpose.

Ownership: A returned String is owned by the caller. No C handle is allocated.

Errors: None for the integer queries. version copies library text.

Units: text

Example: `docs/examples/minimal_rust.rs`

### vertices

Purpose: Atom indices of one ring.

Parameters: index

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: dimensionless indices

Example: `docs/examples/minimal_rust.rs`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/rust/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug rust` writes the 2400 by 1800 PNG.
```

