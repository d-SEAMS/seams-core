<!-- https://docs.dseams.info/reference/rust-api.md -->

# Rust {#rust}

Item docs from the `dseams` crate, taken from [dseams/src/lib.rs](https://github.com/d-SEAMS/dseams-rs/blob/main/dseams/src/lib.rs). `dseams` is the safe API. `dseams-sys` is the raw FFI, re-exported as `seams::ffi`.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### seams::Axis {#seams-axis}

```rust
pub enum Axis
```

### seams::BondOptions {#seams-bondoptions}

```rust
pub struct BondOptions
```

### seams::Brine {#seams-brine}

```rust
pub struct Brine
```

Owned planar ice-front measurement.  ``` let (xyz, edge) = seams::diamond8(); let mut system = seams::System::new()?; system.set_ortho(&xyz, edge, edge, edge)?; let ice = vec!\[1; 8\]; let front = system.brine(&ice, &seams::BrineOptions::default())?; let pockets = front.pockets()?; assert!(pockets.pockets >= 0); # Ok::<(), seams::Error>(()) ```

### seams::BrineFronts {#seams-brinefronts}

```rust
pub struct BrineFronts
```

### seams::BrineOptions {#seams-brineoptions}

```rust
pub struct BrineOptions
```

### seams::BrineProfile {#seams-brineprofile}

```rust
pub struct BrineProfile
```

### seams::CageCounts {#seams-cagecounts}

```rust
pub struct CageCounts
```

### seams::CageLabel {#seams-cagelabel}

```rust
pub enum CageLabel
```

### seams::CageOptions {#seams-cageoptions}

```rust
pub struct CageOptions
```

### seams::Cages {#seams-cages}

```rust
pub struct Cages
```

Owned hexagonal and double-diamond membership.  ``` let (xyz, edge) = seams::diamond8(); let mut system = seams::System::new()?; system.set_ortho(&xyz, edge, edge, edge)?; let cages = system.cages(&seams::CageOptions::default())?; let counts = cages.counts()?; assert_eq!(counts.hexagonal + counts.cubic + counts.water + counts.both, 8); # Ok::<(), seams::Error>(()) ```

### seams::ChillOptions {#seams-chilloptions}

```rust
pub struct ChillOptions
```

### seams::ChillRule {#seams-chillrule}

```rust
pub enum ChillRule
```

### seams::Classification {#seams-classification}

```rust
pub struct Classification
```

### seams::Contacts {#seams-contacts}

```rust
pub struct Contacts
```

### seams::DensityOptions {#seams-densityoptions}

```rust
pub struct DensityOptions
```

### seams::DensityProfile {#seams-densityprofile}

```rust
pub struct DensityProfile
```

### seams::Distribution {#seams-distribution}

```rust
pub struct Distribution
```

### seams::DistributionOptions {#seams-distributionoptions}

```rust
pub struct DistributionOptions
```

### seams::DomainOptions {#seams-domainoptions}

```rust
pub struct DomainOptions
```

### seams::DomainReport {#seams-domainreport}

```rust
pub struct DomainReport
```

### seams::DomainSubset {#seams-domainsubset}

```rust
pub enum DomainSubset
```

### seams::Error {#seams-error}

```rust
pub enum Error
```

### seams::Fingerprint {#seams-fingerprint}

```rust
pub struct Fingerprint
```

Owned topological fingerprint.  ``` let (xyz, edge) = seams::diamond8(); let mut system = seams::System::new()?; system.set_ortho(&xyz, edge, edge, edge)?; let fp = system.fingerprint(&seams::FingerprintOptions::default())?; assert!(!fp.key()?.is_empty()); assert_eq!(fp.rings()?.len(), 13); # Ok::<(), seams::Error>(()) ```

### seams::FingerprintOptions {#seams-fingerprintoptions}

```rust
pub struct FingerprintOptions
```

### seams::Graph {#seams-graph}

```rust
pub enum Graph
```

### seams::KnnOptions {#seams-knnoptions}

```rust
pub struct KnnOptions
```

### seams::Label {#seams-label}

```rust
pub enum Label
```

### seams::NeighbourOptions {#seams-neighbouroptions}

```rust
pub struct NeighbourOptions
```

### seams::NeighbourPair {#seams-neighbourpair}

```rust
pub struct NeighbourPair
```

### seams::NeighbourPairs {#seams-neighbourpairs}

```rust
pub struct NeighbourPairs
```

Iterator of directed neighbour pairs.  Yields \[`NeighbourPair`\]. Created by \[`Neighbours::pairs`\].

### seams::Neighbours {#seams-neighbours}

```rust
pub struct Neighbours
```

Owned neighbour list.  The handle is freed on drop. It is \[`Send`\] and not \[`Sync`\].  ``` let (xyz, edge) = seams::diamond8(); let mut system = seams::System::new()?; system.set_ortho(&xyz, edge, edge, edge)?; let list = system.neighbours(&seams::NeighbourOptions::default())?; assert!(list.pairs().all(|pair| pair.atom != pair.neighbour)); assert!(list.degrees()?.iter().all(|degree| *degree == 4)); # Ok::<(), seams::Error>(()) ```

### seams::OrderOptions {#seams-orderoptions}

```rust
pub struct OrderOptions
```

### seams::Pockets {#seams-pockets}

```rust
pub struct Pockets
```

### seams::PrismCensus {#seams-prismcensus}

```rust
pub struct PrismCensus
```

### seams::PrismOptions {#seams-prismoptions}

```rust
pub struct PrismOptions
```

### seams::Q4Q6 {#seams-q4q6}

```rust
pub struct Q4Q6
```

### seams::Q6Bar {#seams-q6bar}

```rust
pub struct Q6Bar
```

### seams::Q6BarOptions {#seams-q6baroptions}

```rust
pub struct Q6BarOptions
```

### seams::Rejection {#seams-rejection}

```rust
pub struct Rejection
```

### seams::RingIter {#seams-ringiter}

```rust
pub struct RingIter
```

Iterator of primitive rings. Each item is a \[`Result`\] of vertex indices.

### seams::RingOptions {#seams-ringoptions}

```rust
pub struct RingOptions
```

### seams::RingSymmetry {#seams-ringsymmetry}

```rust
pub enum RingSymmetry
```

### seams::Rings {#seams-rings}

```rust
pub struct Rings
```

Owned primitive-ring set.  ``` let (xyz, edge) = seams::diamond8(); let mut system = seams::System::new()?; system.set_ortho(&xyz, edge, edge, edge)?; let rings = system.rings(&seams::RingOptions::default())?; assert!(rings.size() > 0); assert!(!rings.vertices(0)?.is_empty()); # Ok::<(), seams::Error>(()) ```

### seams::RunningCoordination {#seams-runningcoordination}

```rust
pub struct RunningCoordination
```

### seams::SheetOptions {#seams-sheetoptions}

```rust
pub struct SheetOptions
```

### seams::SignatureMatch {#seams-signaturematch}

```rust
pub struct SignatureMatch
```

### seams::SignatureOptions {#seams-signatureoptions}

```rust
pub struct SignatureOptions
```

### seams::SiteKind {#seams-sitekind}

```rust
pub enum SiteKind
```

### seams::Status {#seams-status}

```rust
pub enum Status
```

### seams::Steinhardt {#seams-steinhardt}

```rust
pub struct Steinhardt
```

### seams::SteinhardtOptions {#seams-steinhardtoptions}

```rust
pub struct SteinhardtOptions
```

### seams::System {#seams-system}

```rust
pub struct System
```

One frame: coordinates, a cell, and the analyses that read them.  The handle is freed on drop. It is \[`Send`\] and not \[`Sync`\].  ``` let mut system = seams::System::new()?; assert_eq!(system.count(), 0); let err = system.read(“missing-trajectory.lammpstrj”, 1, -1).expect_err(“no such file”); assert!(!err.message().is_empty()); # Ok::<(), seams::Error>(()) ```

### seams::abi {#seams-abi}

```rust
pub fn abi () -> i32
```

### seams::diamond8 {#seams-diamond8}

```rust
pub fn diamond8 () -> (Vec<f64>, f64)
```

### seams::has_chemfiles {#seams-has-chemfiles}

```rust
pub fn has_chemfiles () -> bool
```

### seams::has_readcon {#seams-has-readcon}

```rust
pub fn has_readcon () -> bool
```

### seams::label_name {#seams-label-name}

```rust
pub fn label_name (label: Label) -> String
```

### seams::nframes {#seams-nframes}

```rust
pub fn nframes (path: &str) -> Result<usize, Error>
```

Count frames in a trajectory path.  ``` let err = seams::nframes(”\\0”).expect_err(“a nul cannot be a path”); assert_eq!(err.code(), seams_sys::SEAMS_ERR_ARG); ```

### seams::num_threads {#seams-num-threads}

```rust
pub fn num_threads () -> i32
```

### seams::ring_symmetry_name {#seams-ring-symmetry-name}

```rust
pub fn ring_symmetry_name (symmetry: RingSymmetry) -> String
```

### seams::set_num_threads {#seams-set-num-threads}

```rust
pub fn set_num_threads (n: i32) -> Result<(), Error>
```

Set the OpenMP width for the next call. Zero restores the library default.  The width is process-global. Do not call this while another thread is inside the library.  ``` seams::set_num_threads(1)?; assert_eq!(seams::num_threads(), 1); # Ok::<(), seams::Error>(()) ```

### seams::version {#seams-version}

```rust
pub fn version () -> String
```

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| seams::Axis | no direct C call |
| seams::BondOptions | no direct C call |
| seams::Brine | seams_brine, seams_brine_destroy |
| seams::BrineFronts | no direct C call |
| seams::BrineOptions | no direct C call |
| seams::BrineProfile | no direct C call |
| seams::CageCounts | no direct C call |
| seams::CageLabel | no direct C call |
| seams::CageOptions | no direct C call |
| seams::Cages | seams_cage_set, seams_cages_destroy |
| seams::ChillOptions | no direct C call |
| seams::ChillRule | no direct C call |
| seams::Classification | no direct C call |
| seams::Contacts | no direct C call |
| seams::DensityOptions | no direct C call |
| seams::DensityProfile | no direct C call |
| seams::Distribution | no direct C call |
| seams::DistributionOptions | no direct C call |
| seams::DomainOptions | no direct C call |
| seams::DomainReport | no direct C call |
| seams::DomainSubset | seams_handle_error |
| seams::Error | no direct C call |
| seams::Fingerprint | seams_fingerprint_destroy, seams_topo |
| seams::FingerprintOptions | no direct C call |
| seams::Graph | no direct C call |
| seams::KnnOptions | no direct C call |
| seams::Label | no direct C call |
| seams::NeighbourOptions | no direct C call |
| seams::NeighbourPair | no direct C call |
| seams::NeighbourPairs | no direct C call |
| seams::Neighbours | seams_neighbours |
| seams::OrderOptions | no direct C call |
| seams::Pockets | no direct C call |
| seams::PrismCensus | no direct C call |
| seams::PrismOptions | no direct C call |
| seams::Q4Q6 | no direct C call |
| seams::Q6Bar | no direct C call |
| seams::Q6BarOptions | no direct C call |
| seams::Rejection | no direct C call |
| seams::RingIter | no direct C call |
| seams::RingOptions | no direct C call |
| seams::RingSymmetry | no direct C call |
| seams::Rings | seams_rings, seams_rings_destroy |
| seams::RunningCoordination | no direct C call |
| seams::SheetOptions | no direct C call |
| seams::SignatureMatch | no direct C call |
| seams::SignatureOptions | no direct C call |
| seams::SiteKind | no direct C call |
| seams::Status | no direct C call |
| seams::Steinhardt | no direct C call |
| seams::SteinhardtOptions | no direct C call |
| seams::System | seams_system, seams_system_destroy |
| seams::abi | seams_abi |
| seams::diamond8 | no direct C call |
| seams::has_chemfiles | no direct C call |
| seams::has_readcon | no direct C call |
| seams::label_name | seams_label_name |
| seams::nframes | seams_path_nframes |
| seams::num_threads | seams_get_num_threads |
| seams::ring_symmetry_name | seams_ring_symmetry_name |
| seams::set_num_threads | seams_set_num_threads |
| seams::version | seams_version_string |
