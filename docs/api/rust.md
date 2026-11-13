# Rust

`seams-sys` is the raw FFI crate. `seams` is the safe wrapper. Build both with
`SEAMS_LIB_DIR` pointing at the directory that contains `libseams.so`.

```rust
let mut sys = seams::System::new()?;
sys.set_ortho(&xyz, lx, ly, lz)?;
let (labels, counts) = sys.chill_plus(3.5, seams_sys::SEAMS_GRAPH_CUTOFF, seams_sys::SEAMS_CHILL_PUBLISHED)?;
let (q4, q6) = sys.q4q6(3.5)?;
let degree = sys.neighbour_degree(3.5)?;
let (six, fraction) = sys.six_ring_fraction(3.5)?;
```

`set_frame` takes the nine-number cell, including a non-zero tilt.
`seams::set_num_threads`, `seams::num_threads`, and `seams::version` match the C API.
`seams::Error` carries the status `code` and the `message` from `seams_last_error`.
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

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### Error

Purpose: Status code and message from a failed call.

Parameters: code; message

Units: code is dimensionless; message is text

Example: `docs/examples/minimal_rust.rs`

### System

Purpose: Owned frame handle.

Parameters: none

Units: handle

Example: `docs/examples/minimal_rust.rs`

### chill_plus

Purpose: CHILL+ labels and nine population counts.

Parameters: cutoff; graph; rule

Units: cutoff in angstrom; graph and rule are dimensionless codes

Example: `docs/examples/minimal_rust.rs`

### count

Purpose: Number of points in the system.

Parameters: self

Units: dimensionless count

Example: `docs/examples/minimal_rust.rs`

### neighbour_degree

Purpose: Degree of each atom on the cutoff graph.

Parameters: cutoff

Units: cutoff in angstrom; degrees are dimensionless

Example: `docs/examples/minimal_rust.rs`

### new

Purpose: Allocate an empty System.

Parameters: none

Units: handle

Example: `docs/examples/minimal_rust.rs`

### num_threads

Purpose: OpenMP width the next call will use.

Parameters: none

Units: dimensionless thread count

Example: `docs/examples/minimal_rust.rs`

### q4q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: cutoff

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_rust.rs`

### set_frame

Purpose: Store points and a nine-number cell.

Parameters: xyz; box9

Units: coordinates and box lengths in angstrom

Example: `docs/examples/minimal_rust.rs`

### set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: n

Units: dimensionless thread count

Example: `docs/examples/minimal_rust.rs`

### set_ortho

Purpose: Store points in an orthorhombic cell at the origin.

Parameters: xyz; lx; ly; lz

Units: coordinates and edge lengths in angstrom

Example: `docs/examples/minimal_rust.rs`

### six_ring_fraction

Purpose: Fraction of atoms on a primitive six-ring.

Parameters: cutoff

Units: cutoff in angstrom; fraction is dimensionless

Example: `docs/examples/minimal_rust.rs`

### version

Purpose: Library version text.

Parameters: none

Units: text

Example: `docs/examples/minimal_rust.rs`
<!-- API_SYMBOLS_END -->
