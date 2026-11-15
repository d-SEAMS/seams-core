# Rust

`seams-sys` is the raw FFI crate. `seams` is the safe crate. A handle is an
owned value with `Drop`. A fallible call returns `Result<T, seams::Error>`.
The error text is `seams_handle_error`. Coordinates and cutoffs are angstroms.
A cell is nine LAMMPS numbers: origin, edge lengths, then `xy`, `xz`, `yz`.

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

To link a meson build instead, point `SEAMS_LIB_DIR` at the directory that
contains `libseams.so` before `cargo build`. The same directory belongs on
`LD_LIBRARY_PATH` at run time.

## Quickstart

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
    let ice = sys.chill_plus(&seams::ChillOptions::default())?;
    let order = sys.q4q6(3.5)?;
    let degree = sys.neighbour_degree(3.5)?;
    let (six, _fraction) = sys.six_ring_fraction(3.5)?;
    assert_eq!(ice.labels.len(), 8);
    assert_eq!(ice.counts.iter().sum::<i32>(), 8);
    assert_eq!(order.q4.len(), 8);
    assert_eq!(order.q6.len(), 8);
    assert_eq!(degree.len(), 8);
    assert!(six >= 0);
    Ok(())
}
```

```bash
export SEAMS_LIB_DIR=$PWD/build/src
export LD_LIBRARY_PATH=$PWD/build/src
cargo test --manifest-path bindings/rust/Cargo.toml -- --test-threads=1
cargo test --manifest-path bindings/rust/Cargo.toml --examples -- --test-threads=1
cargo test --manifest-path bindings/rust/Cargo.toml --doc
cargo clippy --manifest-path bindings/rust/Cargo.toml --all-targets -- -D warnings
cargo fmt --manifest-path bindings/rust/Cargo.toml --all -- --check
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

The eight-atom cubic cell is `seams::diamond8`. CHILL+ labels every atom cubic.

```rust
let (xyz, edge) = seams::diamond8();
seams::set_num_threads(1)?;
let mut system = seams::System::new()?;
system.set_ortho(&xyz, edge, edge, edge)?;
let ice = system.chill_plus(&seams::ChillOptions::default())?;
assert_eq!(ice.counts[seams::Label::Cubic as usize], 8);
```

`docs/examples/minimal_rust.rs` reads the GenIce sI cell (46 atoms, all
clathrate, 54 six-rings). `cargo test` runs it through
`docs/examples/run_binding.py`.

```{literalinclude} ../../../docs/examples/minimal_rust.rs
:language: rust
```

The 64-atom diamond in `bindings/rust/seams/tests/ice.rs` checks both CHILL+
rules, a q6 mean above 0.3, four neighbours, and a six-ring fraction above
one half.

```{literalinclude} ../../../bindings/rust/seams/tests/ice.rs
:language: rust
:lines: 29-63
```

`bindings/rust/seams-sys` publishes each name from `include/seams.h`. The one-line summary, parameters, ownership, and error codes for those names are the C symbol reference.

## Concepts

A `System` owns one frame. `set_ortho` stores an orthorhombic cell at the
origin. `set_frame` stores the nine-number cell, including a tilt. `read`
loads one 1-based frame from XYZ, CON, chemfiles, or a LAMMPS dump. `count`
is the number of points.

Analyses with several knobs take an options struct and `Default`. `Graph`
is `Cutoff`, `Mutual`, or `Seeded`. `ChillRule` is `Published` or `Ovito`.
`Label`, `CageLabel`, `SiteKind`, `RingSymmetry`, `Axis`, and `DomainSubset`
are the C enumerators as Rust enums.

`Neighbours::pairs` yields `NeighbourPair`. `Rings::iter` yields each ring's
vertex indices. `Cages`, `Fingerprint`, and `Brine` are owned handles. Their
destructors free the C object.

A handle is `Send`. It is not `Sync`: the C library requires exclusive use of
a handle. Two frames on two threads may run `steinhardt`, `ql`, `q4q6`, or
`q6bar` at the same time; the Steinhardt scratch belongs to the calling thread.
`set_num_threads` is process-global.

`cargo doc --manifest-path bindings/rust/Cargo.toml --no-deps` is the rustdoc.
With `SEAMS_LIB_DIR` set and `sphinxcontrib-rust` installed, Sphinx renders
that rustdoc under `crates/seams/lib`, the same extension rgpot uses for
`rgpot-core`.

## Ice Ih

`input/traj/ice1h.lammpstrj`, frame 1, every atom, cutoff 3.5 angstrom, one
thread. The safe crate and the C++ engine print the same census. Hexagonal
is 128. The mean Lechner-Dellago q6-bar is 0.562448843702.

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

`chill0` is the nine CHILL populations in label order. Cubic, hexagonal,
water, interfacial, clathrate, interfacial clathrate, unclassified,
reclassified cubic, reclassified hexagonal. q3 through q12 are means of the
local parameter. q6bar is the mean of the averaged q6. `vert` is the length
of ring 0 and its first atom. `cages` is hexagonal, cubic, water, both.
`knn` is the mutual four-neighbour pair count. `rdf` is the bin count, the
cell volume in cubic angstroms, and g of the first bin. `rcn` is the last
running coordination. `cn` is the coordination to 3.5 angstroms. `fp` is the
fingerprint key length and the ring-census sum. `site` is the kind of atom 0
under `1=waterO` (water oxygen, code 8). `pairs` is contacts, cations, anions
under `1=cationHead`. `domain` is the apolar subset size, the largest
component, and the percolation measure. `prism` is perfect six-stacks.
`sheet` is six-rings in the monolayer histogram.

sI (`input/traj/genice_sI.lammpstrj`, oxygen only) has 46 atoms and 2 cages
of signature `512` covering 40 atoms. The XYZ file has 256 atoms and a
coordination of 43.325490196078 out to 6 angstroms. The CON file has 4 atoms
and 0 hydrogen bonds at the default geometry. Triclinic sH has 136 atoms and
a mean q6 of 0.625601998379.

`symmetry(0.0)` on this ice Ih frame labels all 128 shells as ice Ih.
`brine` with every atom marked ice reports 0 pockets. The C++ example
prints `shell_ih 128 pockets 0` for the same frame. Both calls are on
`System`. The doctests cover them on the eight-atom cell.

## Errors and ownership

`Error` is an enum. `Library` carries a `Status` and the handle message.
`Argument` is a value the C boundary cannot accept, such as a nul inside a
path. `code()` is the status integer. `Display` prints the message and the
status name.

```rust
let mut system = seams::System::new()?;
let err = system.chill_plus(&seams::ChillOptions::default()).unwrap_err();
assert_eq!(err.status(), Some(seams::Status::State));
```

`System`, `Neighbours`, `Rings`, `Cages`, `Fingerprint`, and `Brine` own
their C pointers. Moving a handle transfers that ownership. Dropping it
calls the matching `seams_*_destroy`. The raw pointer is not part of the
public API.

`bindings/rust/seams/tests/cpp_match.rs` reruns every analysis above against
`tests/binding_oracle.cpp` and requires the same numbers.

<!-- API_SYMBOLS_BEGIN -->

## Symbol reference

### Axis

Purpose: Axis

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### BondOptions

Purpose: BondOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Brine

Purpose: Brine

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### BrineFronts

Purpose: BrineFronts

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### BrineOptions

Purpose: BrineOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### BrineProfile

Purpose: Concentration profile along the brine axis.

Parameters: named in the rustdoc on this item.

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### CageCounts

Purpose: CageCounts

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### CageLabel

Purpose: CageLabel

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### CageOptions

Purpose: CageOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Cages

Purpose: Owned hexagonal and double-diamond membership.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### ChillOptions

Purpose: ChillOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### ChillRule

Purpose: ChillRule

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Classification

Purpose: Classification

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Contacts

Purpose: Contacts

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### DensityOptions

Purpose: DensityOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### DensityProfile

Purpose: DensityProfile

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Distribution

Purpose: Distribution

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### DistributionOptions

Purpose: DistributionOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### DomainOptions

Purpose: DomainOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### DomainReport

Purpose: DomainReport

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### DomainSubset

Purpose: DomainSubset

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Error

Purpose: Cell axis.

Parameters: named in the rustdoc on this item.

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

### FingerprintOptions

Purpose: FingerprintOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Graph

Purpose: Graph

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### KnnOptions

Purpose: KnnOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Label

Purpose: Label

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### NeighbourOptions

Purpose: NeighbourOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### NeighbourPair

Purpose: NeighbourPair

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### NeighbourPairs

Purpose: NeighbourPairs

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Neighbours

Purpose: Owned neighbour-list handle.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### OrderOptions

Purpose: OrderOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Pockets

Purpose: Pockets

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### PrismCensus

Purpose: PrismCensus

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### PrismOptions

Purpose: PrismOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Q4Q6

Purpose: Q4Q6

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Q6Bar

Purpose: Q6Bar

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Q6BarOptions

Purpose: Q6BarOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Rejection

Purpose: Rejection

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### RingIter

Purpose: RingIter

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### RingOptions

Purpose: RingOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### RingSymmetry

Purpose: RingSymmetry

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Rings

Purpose: Owned primitive-ring handle.

Parameters: none

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### RunningCoordination

Purpose: RunningCoordination

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### SheetOptions

Purpose: SheetOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### SignatureMatch

Purpose: SignatureMatch

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### SignatureOptions

Purpose: SignatureOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### SiteKind

Purpose: SiteKind

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Status

Purpose: Status

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### Steinhardt

Purpose: Steinhardt

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### SteinhardtOptions

Purpose: SteinhardtOptions

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### System

Purpose: Knobs for a hydrogen-bond count.

Parameters: named in the rustdoc on this item.

Returns: The owned handle.

Ownership: Drop frees the C object (seams_system_destroy or the matching destroy function).

Errors: Constructors return Err(Error) when the C call fails. Error.code is the seams_status.

Units: handle

Example: `docs/examples/minimal_rust.rs`

### abi

Purpose: Abi

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### atom

Purpose: Atom

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### brine

Purpose: Brine

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

Purpose: Owned planar ice-front measurement.

Parameters: named in the rustdoc on this item.

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; graph and rule are dimensionless codes

Example: `docs/examples/minimal_rust.rs`

### code

Purpose: Code

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

Purpose: The two planar crossings of an ice slab.

Parameters: named in the rustdoc on this item.

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

### degrees

Purpose: Degrees

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### density

Purpose: Density

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### diamond8

Purpose: Diamond8

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

Parameters: named in the rustdoc on this item.

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; counts are dimensionless

Example: `docs/examples/minimal_rust.rs`

### fraction

Purpose: Fraction

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### fronts

Purpose: Fronts

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### hbonds

Purpose: Hydrogen-bond count.

Parameters: cutoff; oxygen_type; hydrogen_type; dist; angle

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: distances in angstrom; angle in degrees; count is dimensionless

Example: `docs/examples/minimal_rust.rs`

### histogram

Purpose: Histogram

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### iter

Purpose: Iter

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

### label_name

Purpose: Label_name

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### labels

Purpose: Labels

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### message

Purpose: Message

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

Purpose: Knobs for a planar ice front.

Parameters: named in the rustdoc on this item.

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; degrees are dimensionless

Example: `docs/examples/minimal_rust.rs`

### neighbours

Purpose: Neighbours

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### new

Purpose: New

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### nframes

Purpose: Nframes

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

### pairs

Purpose: Pairs

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### pockets

Purpose: Pockets

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### prisms

Purpose: Prisms

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### profile

Purpose: Profile

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### q4q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: cutoff

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; q4 and q6 are dimensionless

Example: `docs/examples/minimal_rust.rs`

### q6bar

Purpose: Q6bar

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### ql

Purpose: Ql

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

### rejection

Purpose: Rejection

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### ring_symmetry_name

Purpose: Ring_symmetry_name

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

### set_atoms

Purpose: Set_atoms

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

Purpose: Per-atom site kind.

Parameters: named in the rustdoc on this item.

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

### size

Purpose: Size

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### status

Purpose: Status

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

Example: `docs/examples/minimal_rust.rs`

### steinhardt

Purpose: Local q_l and, on request, the Lechner-Dellago average.

Parameters: cutoff; order; average

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: cutoff in angstrom; q_l is dimensionless

Example: `docs/examples/minimal_rust.rs`

### symmetry

Purpose: Symmetry

Parameters: none

Returns: Result of the value, or Err(Error).

Ownership: Returned Vec and String values are owned by the caller. Handles returned as structs are freed by Drop.

Errors: Error.code is the seams_status. Error.message is seams_handle_error when the call had a handle, otherwise seams_last_error.

Units: as named by the parameters

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

