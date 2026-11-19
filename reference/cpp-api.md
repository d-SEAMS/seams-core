<!-- https://docs.dseams.info/reference/cpp-api.md -->

# C++ {#c}

`include/seams.hpp` is the native C++ API: move-only handles in namespace `seams`, with results returned by value. [C++](../stable/cpp.md) is the narrative for that header. This page lists the public types and methods, and the `seams.h` functions each method calls.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### seams::BondOptions {#seams-bondoptions}

```cpp
struct BondOptions
```

Knobs for a hydrogen-bond count.

### seams::Brine {#seams-brine}

```cpp
class Brine : public detail::Handle<seams_brine, seams_brine_destroy>
```

Owned planar ice-front measurement.

### seams::Brine::fronts {#seams-brine-fronts}

```cpp
[[nodiscard]] BrineFronts fronts() const
```

### seams::Brine::pockets {#seams-brine-pockets}

```cpp
[[nodiscard]] Pockets pockets() const
```

### seams::Brine::profile {#seams-brine-profile}

```cpp
[[nodiscard]] BrineProfile profile() const
```

### seams::Brine::rejection {#seams-brine-rejection}

```cpp
[[nodiscard]] Rejection rejection() const
```

### seams::BrineFronts {#seams-brinefronts}

```cpp
struct BrineFronts
```

The two planar crossings of an ice slab.

### seams::BrineOptions {#seams-brineoptions}

```cpp
struct BrineOptions
```

Knobs for a planar ice front.

### seams::BrineProfile {#seams-brineprofile}

```cpp
struct BrineProfile
```

Concentration profile along the brine axis.

### seams::CageCounts {#seams-cagecounts}

```cpp
struct CageCounts
```

Disjoint cage populations.

### seams::CageOptions {#seams-cageoptions}

```cpp
struct CageOptions
```

Knobs for hexagonal and double-diamond cages.

### seams::Cages {#seams-cages}

```cpp
class Cages : public detail::Handle<seams_cage_set, seams_cages_destroy>
```

Owned hexagonal and double-diamond membership.

### seams::Cages::counts {#seams-cages-counts}

```cpp
[[nodiscard]] CageCounts counts() const
```

Disjoint population counts.

### seams::Cages::labels {#seams-cages-labels}

```cpp
[[nodiscard]] std::vector<int> labels(std::size_t n_atoms) const
```

One cage code per atom.

### seams::ChillOptions {#seams-chilloptions}

```cpp
struct ChillOptions
```

Knobs for CHILL and CHILL+.

### seams::ChordlessCycles {#seams-chordlesscycles}

```cpp
struct ChordlessCycles
```

Chordless cycles of one length.

### seams::ChordlessCycles::count {#seams-chordlesscycles-count}

```cpp
[[nodiscard]] std::size_t count() const
```

### seams::ChordlessCycles::cycle {#seams-chordlesscycles-cycle}

```cpp
[[nodiscard]] std::vector<int> cycle(std::size_t index) const
```

### seams::Classification {#seams-classification}

```cpp
struct Classification
```

CHILL or CHILL+ labels and the nine population counts.

### seams::Contacts {#seams-contacts}

```cpp
struct Contacts
```

Mutual nearest unlike contacts of a site table.

### seams::DensityOptions {#seams-densityoptions}

```cpp
struct DensityOptions
```

Knobs for the axial number density.

### seams::DensityProfile {#seams-densityprofile}

```cpp
struct DensityProfile
```

Number density along one cell axis.

### seams::Distribution {#seams-distribution}

```cpp
struct Distribution
```

Partial radial distribution.

### seams::DistributionOptions {#seams-distributionoptions}

```cpp
struct DistributionOptions
```

Knobs for a partial RDF, its running coordination, or the integral.

### seams::DomainOptions {#seams-domainoptions}

```cpp
struct DomainOptions
```

Knobs for a domain search.

### seams::DomainReport {#seams-domainreport}

```cpp
struct DomainReport
```

Largest polar or apolar domain.

### seams::Error {#seams-error}

```cpp
class Error : public std::runtime_error
```

Failure from the C API, or from a rejected argument.

### seams::Error::Error {#seams-error-error}

```cpp
Error(int status, std::string message) : std::runtime_error(std::move(message)), code_(status)
```

Build an error from a status and a message.

### seams::Error::code {#seams-error-code}

```cpp
[[nodiscard]] int code() const noexcept
```

Status integer. Same values as `seams_status`.

### seams::Error::status {#seams-error-status}

```cpp
[[nodiscard]] Status status() const noexcept
```

Status as the `Status` enumeration when the code is known.

### seams::Fingerprint {#seams-fingerprint}

```cpp
class Fingerprint : public detail::Handle<seams_topo, seams_fingerprint_destroy>
```

Owned topological fingerprint.

### seams::Fingerprint::atom {#seams-fingerprint-atom}

```cpp
[[nodiscard]] std::string atom(int index) const
```

Local key of one atom.

### seams::Fingerprint::key {#seams-fingerprint-key}

```cpp
[[nodiscard]] std::string key() const
```

Frame fingerprint text.

### seams::Fingerprint::rings {#seams-fingerprint-rings}

```cpp
[[nodiscard]] std::vector<int> rings() const
```

Primitive-ring census stored on the fingerprint.

### seams::FingerprintOptions {#seams-fingerprintoptions}

```cpp
struct FingerprintOptions
```

Knobs for a topological fingerprint.

### seams::KnnOptions {#seams-knnoptions}

```cpp
struct KnnOptions
```

Knobs for a k-nearest neighbour list.

### seams::NeighbourOptions {#seams-neighbouroptions}

```cpp
struct NeighbourOptions
```

Knobs for a neighbour list.

### seams::Neighbours {#seams-neighbours}

```cpp
class Neighbours : public detail::Handle<seams_neighbours, seams_neighbours_destroy>
```

Owned neighbour list.

### seams::Neighbours::chordless {#seams-neighbours-chordless}

```cpp
[[nodiscard]] ChordlessCycles chordless(int size) const
```

Chordless cycles of one length on this graph.

### seams::Neighbours::degrees {#seams-neighbours-degrees}

```cpp
[[nodiscard]] std::vector<int> degrees() const
```

Degree of each atom.

### seams::Neighbours::nnz {#seams-neighbours-nnz}

```cpp
[[nodiscard]] std::size_t nnz() const
```

Number of neighbour pairs, omitting self.

### seams::Neighbours::pairs {#seams-neighbours-pairs}

```cpp
[[nodiscard]] std::vector<std::pair<int, int>> pairs() const
```

Directed pairs `(atom, neighbour)`.

### seams::OrderOptions {#seams-orderoptions}

```cpp
struct OrderOptions
```

Knobs for q_l on a named graph.

### seams::Pockets {#seams-pockets}

```cpp
struct Pockets
```

Brine pockets trapped in the ice.

### seams::PrismCensus {#seams-prismcensus}

```cpp
struct PrismCensus
```

Perfect and imperfect prism counts by ring size.

### seams::PrismOptions {#seams-prismoptions}

```cpp
struct PrismOptions
```

Knobs for prism stacks.

### seams::Q4Q6 {#seams-q4q6}

```cpp
struct Q4Q6
```

Local q4 and q6 on the cutoff shell.

### seams::Q6Bar {#seams-q6bar}

```cpp
struct Q6Bar
```

Lechner-Dellago q6-bar and the ice count above a cut.

### seams::Q6BarOptions {#seams-q6baroptions}

```cpp
struct Q6BarOptions
```

Knobs for the q6-bar ice summary.

### seams::Rejection {#seams-rejection}

```cpp
struct Rejection
```

Molal rejection coefficient.

### seams::RingOptions {#seams-ringoptions}

```cpp
struct RingOptions
```

Knobs for primitive rings.

### seams::Rings {#seams-rings}

```cpp
class Rings : public detail::Handle<seams_rings, seams_rings_destroy>
```

Owned primitive-ring set.

### seams::Rings::fraction {#seams-rings-fraction}

```cpp
[[nodiscard]] double fraction(int ring_size) const
```

Fraction of atoms that sit on a ring of one size.

### seams::Rings::histogram {#seams-rings-histogram}

```cpp
[[nodiscard]] std::vector<int> histogram(int up_to) const
```

Histogram of ring size.

### seams::Rings::iterator {#seams-rings-iterator}

```cpp
iterator(const Rings *owner, int index) : owner_(owner), index_(index)
```

### seams::Rings::size {#seams-rings-size}

```cpp
[[nodiscard]] int size() const
```

Number of primitive rings.

### seams::Rings::vertices {#seams-rings-vertices}

```cpp
[[nodiscard]] std::vector<int> vertices(int index) const
```

Atom indices of one ring.

### seams::RunningCoordination {#seams-runningcoordination}

```cpp
struct RunningCoordination
```

Running coordination number along a partial RDF.

### seams::SheetOptions {#seams-sheetoptions}

```cpp
struct SheetOptions
```

Knobs for a monolayer ring histogram.

### seams::SignatureMatch {#seams-signaturematch}

```cpp
struct SignatureMatch
```

Closed cages whose ring census matches a signature.

### seams::SignatureOptions {#seams-signatureoptions}

```cpp
struct SignatureOptions
```

Knobs for a signature cage census.

### seams::Steinhardt {#seams-steinhardt}

```cpp
struct Steinhardt
```

Local Steinhardt parameter, and the average when it was asked for.

### seams::SteinhardtOptions {#seams-steinhardtoptions}

```cpp
struct SteinhardtOptions
```

Knobs for a local Steinhardt parameter on the cutoff shell.

### seams::System {#seams-system}

```cpp
class System : public detail::Handle<seams_system, seams_system_destroy>
```

One frame: coordinates, a cell, and the analyses that read them.

### seams::System::System {#seams-system-system}

```cpp
System() : detail::Handle<seams_system, seams_system_destroy>(seams_system_create())
```

Allocate an empty frame.

### seams::System::brine {#seams-system-brine}

```cpp
[[nodiscard]] Brine brine(std::span<const int> ice, const BrineOptions &options)
```

Planar ice front from an ice assignment.

### seams::System::cages {#seams-system-cages}

```cpp
[[nodiscard]] Cages cages(const CageOptions &options)
```

### seams::System::chill {#seams-system-chill}

```cpp
std::pair<std::vector<int>, std::vector<int>> chill(double cutoff, int type, int graph, int rule)
```

### seams::System::chill_plus {#seams-system-chill-plus}

```cpp
[[nodiscard]] Classification chill_plus(const ChillOptions &options)
```

### seams::System::contact_pairs {#seams-system-contact-pairs}

```cpp
[[nodiscard]] Contacts contact_pairs(std::string_view spec)
```

### seams::System::coordination {#seams-system-coordination}

```cpp
[[nodiscard]] double coordination(const DistributionOptions &options)
```

### seams::System::count {#seams-system-count}

```cpp
[[nodiscard]] int count() const
```

Number of points in the frame.

### seams::System::density {#seams-system-density}

```cpp
[[nodiscard]] DensityProfile density(const DensityOptions &options)
```

### seams::System::domains {#seams-system-domains}

```cpp
[[nodiscard]] DomainReport domains(const DomainOptions &options)
```

### seams::System::fingerprint {#seams-system-fingerprint}

```cpp
[[nodiscard]] Fingerprint fingerprint(const FingerprintOptions &options)
```

### seams::System::hbonds {#seams-system-hbonds}

```cpp
[[nodiscard]] int hbonds(const BondOptions &options)
```

### seams::System::knn {#seams-system-knn}

```cpp
[[nodiscard]] Neighbours knn(const KnnOptions &options)
```

### seams::System::monolayer {#seams-system-monolayer}

```cpp
[[nodiscard]] std::vector<int> monolayer(const SheetOptions &options)
```

Primitive-ring histogram of a monolayer.

### seams::System::neighbour_degree {#seams-system-neighbour-degree}

```cpp
[[nodiscard]] std::vector<int> neighbour_degree(double cutoff, int species = 0)
```

Degree of each atom on the cutoff graph.

### seams::System::neighbours {#seams-system-neighbours}

```cpp
[[nodiscard]] Neighbours neighbours(const NeighbourOptions &options)
```

### seams::System::prisms {#seams-system-prisms}

```cpp
[[nodiscard]] PrismCensus prisms(const PrismOptions &options)
```

### seams::System::q4q6 {#seams-system-q4q6}

```cpp
[[nodiscard]] Q4Q6 q4q6(double cutoff, int species = 0)
```

### seams::System::q6bar {#seams-system-q6bar}

```cpp
[[nodiscard]] Q6Bar q6bar(const Q6BarOptions &options)
```

### seams::System::ql {#seams-system-ql}

```cpp
[[nodiscard]] Steinhardt ql(const OrderOptions &options)
```

### seams::System::rdf {#seams-system-rdf}

```cpp
[[nodiscard]] Distribution rdf(const DistributionOptions &options)
```

### seams::System::read {#seams-system-read}

```cpp
void read(std::string_view path, int frame, int type_filter)
```

Read one frame from a trajectory path.

### seams::System::rings {#seams-system-rings}

```cpp
[[nodiscard]] Rings rings(const RingOptions &options)
```

### seams::System::rings_chordless {#seams-system-rings-chordless}

```cpp
[[nodiscard]] ChordlessCycles rings_chordless(const NeighbourOptions &options, int size)
```

Chordless cycles of one length on a neighbour graph.

### seams::System::running_cn {#seams-system-running-cn}

```cpp
[[nodiscard]] RunningCoordination running_cn(const DistributionOptions &options)
```

### seams::System::set_ortho {#seams-system-set-ortho}

```cpp
void set_ortho(std::span<const double> xyz, double lx, double ly, double lz)
```

Store points in an orthorhombic cell whose origin is zero.

### seams::System::signature_cages {#seams-system-signature-cages}

```cpp
[[nodiscard]] SignatureMatch signature_cages(const SignatureOptions &options)
```

### seams::System::site_table {#seams-system-site-table}

```cpp
[[nodiscard]] std::vector<int> site_table(std::string_view spec)
```

Per-atom site kind.

### seams::System::steinhardt {#seams-system-steinhardt}

```cpp
std::pair<std::vector<double>, std::vector<double>> steinhardt(double cutoff, int order, bool average)
```

### seams::System::symmetry {#seams-system-symmetry}

```cpp
[[nodiscard]] std::vector<int> symmetry(double shell_rmsd)
```

Per-atom ice-shell label.

### seams::abi {#seams-abi}

```cpp
[[nodiscard]] inline int abi()
```

ABI number compiled into this library.

### seams::check {#seams-check}

```cpp
inline void check(int status, const void *handle) { if (status == SEAMS_OK)
```

### seams::label_name {#seams-label-name}

```cpp
[[nodiscard]] inline std::string_view label_name(Label label)
```

Short name of a CHILL+ label.

### seams::nframes {#seams-nframes}

```cpp
[[nodiscard]] inline std::size_t nframes(std::string_view path)
```

Count frames in a trajectory path.

### seams::num_threads {#seams-num-threads}

```cpp
[[nodiscard]] inline int num_threads()
```

OpenMP width the next call will use.

### seams::ring_symmetry_name {#seams-ring-symmetry-name}

```cpp
[[nodiscard]] inline std::string_view ring_symmetry_name(RingSymmetry symmetry)
```

Short name of a ring-symmetry code.

### seams::set_num_threads {#seams-set-num-threads}

```cpp
inline void set_num_threads(int n)
```

Set the OpenMP width for the next call.

### seams::version {#seams-version}

```cpp
[[nodiscard]] inline std::string version()
```

Library version text, for example `2.14.0`.

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| seams::BondOptions | no direct C call |
| seams::Brine | seams_brine_fronts, seams_brine_pockets, seams_brine_profile, seams_brine_rejection |
| seams::Brine::fronts | seams_brine_fronts |
| seams::Brine::pockets | seams_brine_pockets |
| seams::Brine::profile | seams_brine_profile |
| seams::Brine::rejection | seams_brine_rejection |
| seams::BrineFronts | no direct C call |
| seams::BrineOptions | no direct C call |
| seams::BrineProfile | no direct C call |
| seams::CageCounts | no direct C call |
| seams::CageOptions | no direct C call |
| seams::Cages | seams_cages_counts, seams_cages_labels |
| seams::Cages::counts | seams_cages_counts |
| seams::Cages::labels | seams_cages_labels |
| seams::ChillOptions | no direct C call |
| seams::ChordlessCycles | no direct C call |
| seams::ChordlessCycles::count | no direct C call |
| seams::ChordlessCycles::cycle | no direct C call |
| seams::Classification | no direct C call |
| seams::Contacts | no direct C call |
| seams::DensityOptions | no direct C call |
| seams::DensityProfile | no direct C call |
| seams::Distribution | no direct C call |
| seams::DistributionOptions | no direct C call |
| seams::DomainOptions | no direct C call |
| seams::DomainReport | no direct C call |
| seams::Error | no direct C call |
| seams::Error::Error | no direct C call |
| seams::Error::code | no direct C call |
| seams::Error::status | no direct C call |
| seams::Fingerprint | seams_fingerprint_atom, seams_fingerprint_key, seams_fingerprint_rings |
| seams::Fingerprint::atom | seams_fingerprint_atom |
| seams::Fingerprint::key | seams_fingerprint_key |
| seams::Fingerprint::rings | seams_fingerprint_rings |
| seams::FingerprintOptions | no direct C call |
| seams::KnnOptions | no direct C call |
| seams::NeighbourOptions | no direct C call |
| seams::Neighbours | seams_neighbours_copy, seams_neighbours_count, seams_rings_chordless |
| seams::Neighbours::chordless | seams_neighbours_count, seams_rings_chordless |
| seams::Neighbours::degrees | seams_neighbours_copy |
| seams::Neighbours::nnz | seams_neighbours_copy |
| seams::Neighbours::pairs | seams_neighbours_copy |
| seams::OrderOptions | no direct C call |
| seams::Pockets | no direct C call |
| seams::PrismCensus | no direct C call |
| seams::PrismOptions | no direct C call |
| seams::Q4Q6 | no direct C call |
| seams::Q6Bar | no direct C call |
| seams::Q6BarOptions | no direct C call |
| seams::Rejection | no direct C call |
| seams::RingOptions | no direct C call |
| seams::Rings | seams_rings_count, seams_rings_fraction, seams_rings_histogram, seams_rings_vertices |
| seams::Rings::fraction | seams_rings_fraction |
| seams::Rings::histogram | seams_rings_histogram |
| seams::Rings::iterator | no direct C call |
| seams::Rings::size | seams_rings_count |
| seams::Rings::vertices | seams_rings_vertices |
| seams::RunningCoordination | no direct C call |
| seams::SheetOptions | no direct C call |
| seams::SignatureMatch | no direct C call |
| seams::SignatureOptions | no direct C call |
| seams::Steinhardt | no direct C call |
| seams::SteinhardtOptions | no direct C call |
| seams::System | seams_brine_analyse, seams_cages, seams_chill, seams_chill_plus, seams_contact_pairs, seams_coordination, seams_density_axis, seams_domains, seams_fingerprint, seams_hbonds, seams_monolayer, seams_neighbours_build, seams_neighbours_knn, seams_prisms, seams_q4q6, seams_q6bar, seams_ql, seams_rdf, seams_rings_find, seams_running_cn, seams_signature_cages, seams_site_table, seams_steinhardt, seams_system_count, seams_system_create, seams_system_read, seams_system_set_frame, seams_system_set_ortho, seams_system_symmetry |
| seams::System::System | no direct C call |
| seams::System::brine | seams_brine_analyse |
| seams::System::cages | seams_cages |
| seams::System::chill | seams_system_count |
| seams::System::chill_plus | seams_chill, seams_chill_plus, seams_system_count |
| seams::System::contact_pairs | seams_contact_pairs |
| seams::System::coordination | seams_coordination |
| seams::System::count | seams_system_count |
| seams::System::density | seams_density_axis |
| seams::System::domains | seams_domains |
| seams::System::fingerprint | seams_fingerprint |
| seams::System::hbonds | seams_hbonds |
| seams::System::knn | seams_neighbours_knn |
| seams::System::monolayer | seams_monolayer |
| seams::System::neighbour_degree | seams_neighbours_build |
| seams::System::neighbours | seams_neighbours_build |
| seams::System::prisms | seams_prisms |
| seams::System::q4q6 | seams_q4q6, seams_system_count |
| seams::System::q6bar | seams_q6bar, seams_system_count |
| seams::System::ql | seams_ql, seams_system_count |
| seams::System::rdf | seams_rdf |
| seams::System::read | seams_system_read |
| seams::System::rings | seams_rings_find |
| seams::System::rings_chordless | seams_neighbours_build |
| seams::System::running_cn | seams_running_cn |
| seams::System::set_ortho | no direct C call |
| seams::System::signature_cages | seams_signature_cages |
| seams::System::site_table | seams_site_table, seams_system_count |
| seams::System::steinhardt | seams_ql, seams_steinhardt, seams_system_count |
| seams::System::symmetry | seams_system_count, seams_system_symmetry |
| seams::abi | seams_abi |
| seams::check | no direct C call |
| seams::label_name | seams_label_name |
| seams::nframes | seams_path_nframes |
| seams::num_threads | seams_get_num_threads |
| seams::ring_symmetry_name | seams_ring_symmetry_name |
| seams::set_num_threads | seams_set_num_threads |
| seams::version | seams_version_string |
