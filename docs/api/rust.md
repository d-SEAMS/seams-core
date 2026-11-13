# Rust

`seams-sys` is the raw FFI crate. `seams` is the safe crate. A handle is an
owned value with `Drop`. A fallible call returns `Result<T, seams::Error>`.
The error text is `seams_handle_error`. Coordinates and cutoffs are angstroms.
A cell is nine LAMMPS numbers: origin, edge lengths, then `xy`, `xz`, `yz`.

Point `SEAMS_LIB_DIR` at the directory that contains `libseams.so` before
`cargo build`. The same directory belongs on `LD_LIBRARY_PATH` at run time.

## Quickstart

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

A handle is `Send`. It is not `Sync`. The C library requires exclusive use of
a handle, and the Steinhardt scratch is process-global, so two frames must
not run `steinhardt`, `ql`, `q4q6`, or `q6bar` at the same time.
`set_num_threads` is process-global as well.

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

Purpose: Cell axis.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### BondOptions

Purpose: Knobs for a hydrogen-bond count.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Brine

Purpose: Owned planar ice-front measurement.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### BrineFronts

Purpose: The two planar crossings of an ice slab.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### BrineOptions

Purpose: Knobs for a planar ice front.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### BrineProfile

Purpose: Concentration profile along the brine axis.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### CageCounts

Purpose: Disjoint cage populations.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### CageLabel

Purpose: Per-atom cage membership.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### CageOptions

Purpose: Knobs for hexagonal and double-diamond cages.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Cages

Purpose: Owned hexagonal and double-diamond membership.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### ChillOptions

Purpose: Knobs for CHILL and CHILL+.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### ChillRule

Purpose: CHILL+ threshold table.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Classification

Purpose: CHILL or CHILL+ labels and the nine population counts.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Contacts

Purpose: Mutual nearest unlike contacts of a site table.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### DensityOptions

Purpose: Knobs for the axial number density.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### DensityProfile

Purpose: Number density along one cell axis.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Distribution

Purpose: Partial radial distribution.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### DistributionOptions

Purpose: Knobs for a partial RDF, its running coordination, or the integral.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### DomainOptions

Purpose: Knobs for a domain search.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### DomainReport

Purpose: Largest polar or apolar domain.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### DomainSubset

Purpose: Domain subset.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Error

Purpose: A failed call.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Fingerprint

Purpose: Owned topological fingerprint.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### FingerprintOptions

Purpose: Knobs for a topological fingerprint.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Graph

Purpose: Neighbour graph.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### KnnOptions

Purpose: Knobs for a k-nearest neighbour list.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Label

Purpose: CHILL and CHILL+ atom label.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### NeighbourOptions

Purpose: Knobs for a neighbour list.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### NeighbourPair

Purpose: One directed neighbour pair.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### NeighbourPairs

Purpose: Iterator of directed neighbour pairs.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Neighbours

Purpose: Owned neighbour list.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### OrderOptions

Purpose: Knobs for q_l on a named graph.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Pockets

Purpose: Brine pockets trapped in the ice.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### PrismCensus

Purpose: Perfect and imperfect prism counts by ring size.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### PrismOptions

Purpose: Knobs for prism stacks.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Q4Q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Q6Bar

Purpose: Lechner-Dellago q6-bar and the ice count above a cut.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Q6BarOptions

Purpose: Knobs for the q6-bar ice summary.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Rejection

Purpose: Molal rejection coefficient.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### RingIter

Purpose: Iterator of primitive rings.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### RingOptions

Purpose: Knobs for primitive rings.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### RingSymmetry

Purpose: Primitive-ring symmetry class.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Rings

Purpose: Owned primitive-ring set.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### RunningCoordination

Purpose: Running coordination number along a partial RDF.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### SheetOptions

Purpose: Knobs for a monolayer ring histogram.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### SignatureMatch

Purpose: Closed cages whose ring census matches a signature.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### SignatureOptions

Purpose: Knobs for a signature cage census.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### SiteKind

Purpose: Site-table kind.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Status

Purpose: Status codes.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### Steinhardt

Purpose: Local Steinhardt parameter, and the average when it was requested.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### SteinhardtOptions

Purpose: Knobs for a local Steinhardt parameter on the cutoff shell.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### System

Purpose: One frame: coordinates, a cell, and the analyses that read them.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### abi

Purpose: ABI number compiled into this library.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### atom

Purpose: Local key of one atom.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### brine

Purpose: Planar ice front from an ice assignment.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### cages

Purpose: Hexagonal and double-diamond cage membership.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### chill

Purpose: Original CHILL labels and nine population counts.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### chill_plus

Purpose: CHILL+ labels and nine population counts.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### code

Purpose: Status integer.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### contact_pairs

Purpose: Mutual nearest unlike contacts.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### coordination

Purpose: Coordination number integrated to rmax.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### count

Purpose: Number of points in the frame.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### counts

Purpose: Disjoint population counts.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### degrees

Purpose: Degree of each atom.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### density

Purpose: Number density along one cell axis.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### diamond8

Purpose: Eight-atom cubic diamond used by the crate examples.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### domains

Purpose: Largest polar or apolar domain.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### fingerprint

Purpose: Topological fingerprint of the bonded graph.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### fraction

Purpose: Fraction of atoms that sit on a ring of one size.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### fronts

Purpose: The two planar crossings.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### hbonds

Purpose: Hydrogen-bond count.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### histogram

Purpose: Histogram of ring size.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### iter

Purpose: Walk the rings.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### key

Purpose: Frame fingerprint text.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### knn

Purpose: k-nearest neighbour list.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### label_name

Purpose: Short name of a CHILL+ label.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### labels

Purpose: One cage code per atom.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### message

Purpose: Detail text.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### monolayer

Purpose: Primitive-ring histogram of a monolayer.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### neighbour_degree

Purpose: Degree of each atom on the cutoff graph.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### neighbours

Purpose: Cutoff or mutual neighbour list.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### new

Purpose: Allocate an empty frame.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### nframes

Purpose: Count frames in a trajectory path.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### nnz

Purpose: Number of neighbour pairs, omitting self.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### num_threads

Purpose: OpenMP width the next call will use.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### pairs

Purpose: Directed pairs.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### pockets

Purpose: Pockets trapped in the ice, and the liquid reservoir.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### prisms

Purpose: Perfect and imperfect prism counts by ring size.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### profile

Purpose: Concentration profile along the front axis.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### q4q6

Purpose: Local q4 and q6 on the cutoff shell.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### q6bar

Purpose: q6-bar, its mean and median, and the ice count above the cut.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### ql

Purpose: Local q_l on a named graph.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### rdf

Purpose: Partial radial distribution and the cell volume.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### read

Purpose: Read one frame from a trajectory path.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### rejection

Purpose: Molal rejection coefficient.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### ring_symmetry_name

Purpose: Short name of a ring-symmetry code.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### rings

Purpose: Primitive-ring census stored on the fingerprint.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### running_cn

Purpose: Running coordination number along the partial RDF.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### set_atoms

Purpose: Store points, a cell, and optional species and ids.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### set_frame

Purpose: Store points and a general nine-number cell.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### set_num_threads

Purpose: Set the OpenMP width for the next call.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### set_ortho

Purpose: Store points in an orthorhombic cell whose origin is zero.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### signature_cages

Purpose: Closed cages of one ring census.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### site_table

Purpose: Per-atom site kind.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### six_ring_fraction

Purpose: Count of six-rings and the fraction of atoms that sit on one.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### size

Purpose: Number of primitive rings.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### status

Purpose: Status when the failure came from the library and the code is known.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### steinhardt

Purpose: Local q_l on the cutoff shell.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### symmetry

Purpose: Per-atom ice-shell label.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### version

Purpose: Library version text, for example 2.14.0.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`

### vertices

Purpose: Atom indices of one ring.

Parameters: named in the rustdoc on this item.

Units: named in the rustdoc on this item. Lengths are angstroms.

Example: `docs/examples/minimal_rust.rs`
<!-- API_SYMBOLS_END -->


## Structure still

```{figure} ../tutorials/molecules/rust/overlay-web.png
:alt: Molecule still for this front
:width: 480px

Ice Ih coloured by published CHILL+, drawn with chemparseplot. `python docs/scripts/render_molecules.py --slug rust` writes the 2400 by 1800 PNG.
```

