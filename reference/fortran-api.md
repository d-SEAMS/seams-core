<!-- https://docs.dseams.info/reference/fortran-api.md -->

# Fortran {#fortran}

FORD-style comments (`!>`) and the public type-bound procedures in [dseams-fortran](https://github.com/d-SEAMS/dseams-fortran). `seams_c` is the `bind(C)` layer. `seams` is the module a program uses.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### seams::atom_key {#seams-atom-key}

```fortran
atom_key => topo_atom_key
```

### seams::brine {#seams-brine}

```fortran
brine => frame_brine
```

### seams::cages {#seams-cages}

```fortran
cages => frame_cages
```

### seams::chill {#seams-chill}

```fortran
chill => frame_chill
```

### seams::chill_plus {#seams-chill-plus}

```fortran
chill_plus => frame_chill_plus
```

### seams::close {#seams-close}

```fortran
close => frame_close
```

### seams::columns {#seams-columns}

```fortran
columns => neighbour_columns
```

### seams::contact_pairs {#seams-contact-pairs}

```fortran
contact_pairs => frame_contact_pairs
```

### seams::coordination {#seams-coordination}

```fortran
coordination => frame_coordination
```

### seams::count {#seams-count}

```fortran
count => frame_count
```

### seams::counts {#seams-counts}

```fortran
counts => cage_counts
```

### seams::density {#seams-density}

```fortran
density => frame_density
```

### seams::domains {#seams-domains}

```fortran
domains => frame_domains
```

### seams::fingerprint {#seams-fingerprint}

```fortran
fingerprint => frame_fingerprint
```

### seams::fraction {#seams-fraction}

```fortran
fraction => ring_fraction
```

### seams::fronts {#seams-fronts}

```fortran
fronts => brine_fronts
```

### seams::hbonds {#seams-hbonds}

```fortran
hbonds => frame_hbonds
```

### seams::histogram {#seams-histogram}

```fortran
histogram => ring_histogram
```

### seams::integer {#seams-integer}

```fortran
integer
```

The call succeeded.

### seams::key {#seams-key}

```fortran
key => topo_key
```

### seams::knn {#seams-knn}

```fortran
knn => frame_knn
```

### seams::labels {#seams-labels}

```fortran
labels => cage_labels
```

### seams::monolayer {#seams-monolayer}

```fortran
monolayer => frame_monolayer
```

### seams::neighbours {#seams-neighbours}

```fortran
neighbours => frame_neighbours
```

### seams::offsets {#seams-offsets}

```fortran
offsets => neighbour_offsets
```

### seams::pockets {#seams-pockets}

```fortran
pockets => brine_pockets
```

### seams::prisms {#seams-prisms}

```fortran
prisms => frame_prisms
```

### seams::profile {#seams-profile}

```fortran
profile => brine_profile
```

### seams::q4q6 {#seams-q4q6}

```fortran
q4q6 => frame_q4q6
```

### seams::q6bar {#seams-q6bar}

```fortran
q6bar => frame_q6bar
```

### seams::ql {#seams-ql}

```fortran
ql => frame_ql
```

### seams::rdf {#seams-rdf}

```fortran
rdf => frame_rdf
```

### seams::read {#seams-read}

```fortran
read => frame_read
```

### seams::rejection {#seams-rejection}

```fortran
rejection => brine_rejection
```

### seams::rings {#seams-rings}

```fortran
rings => frame_rings
```

### seams::running_cn {#seams-running-cn}

```fortran
running_cn => frame_running_cn
```

### seams::seams_brine {#seams-seams-brine}

```fortran
seams_brine
```

Planar ice/brine front.

### seams::seams_cages {#seams-seams-cages}

```fortran
seams_cages
```

Hexagonal and double-diamond membership.

### seams::seams_frame {#seams-seams-frame}

```fortran
seams_frame
```

One frame. The finalizer frees the C frame.

### seams::seams_neighbours {#seams-seams-neighbours}

```fortran
seams_neighbours
```

Cutoff or neighbour graph. `offsets` and `columns` allocate copies.

### seams::seams_rings {#seams-seams-rings}

```fortran
seams_rings
```

Primitive rings. `histogram` is 0-based: index k is rings of size k.

### seams::seams_topo {#seams-seams-topo}

```fortran
seams_topo
```

Topological fingerprint.

### seams::set_frame {#seams-set-frame}

```fortran
set_frame => frame_set_frame
```

### seams::set_ortho {#seams-set-ortho}

```fortran
set_ortho => frame_set_ortho
```

### seams::signature {#seams-signature}

```fortran
signature => frame_signature
```

### seams::site_kinds {#seams-site-kinds}

```fortran
site_kinds => frame_site_kinds
```

### seams::sizes {#seams-sizes}

```fortran
sizes => ring_sizes
```

### seams::steinhardt {#seams-steinhardt}

```fortran
steinhardt => frame_steinhardt
```

### seams::symmetry {#seams-symmetry}

```fortran
symmetry => frame_symmetry
```

### seams::vertices {#seams-vertices}

```fortran
vertices => ring_vertices
```

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| seams::atom_key | seams_fingerprint_atom, seams_fingerprint_key |
| seams::brine | seams_brine_analyse, seams_system_count |
| seams::cages | seams_system_count |
| seams::chill | seams_chill, seams_chill_plus, seams_system_count |
| seams::chill_plus | seams_chill, seams_chill_plus, seams_system_count |
| seams::close | no direct C call |
| seams::columns | seams_neighbours_copy |
| seams::contact_pairs | seams_contact_pairs |
| seams::coordination | seams_coordination |
| seams::count | seams_system_count |
| seams::counts | seams_cages_counts |
| seams::density | seams_density_axis |
| seams::domains | seams_domains |
| seams::fingerprint | seams_fingerprint |
| seams::fraction | seams_rings_fraction |
| seams::fronts | seams_brine_fronts |
| seams::hbonds | seams_hbonds |
| seams::histogram | seams_rings_histogram |
| seams::integer | no direct C call |
| seams::key | seams_fingerprint_atom, seams_fingerprint_key |
| seams::knn | seams_neighbours_knn |
| seams::labels | seams_cages_labels |
| seams::monolayer | seams_monolayer |
| seams::neighbours | seams_neighbours_build |
| seams::offsets | seams_neighbours_copy |
| seams::pockets | seams_brine_pockets |
| seams::prisms | seams_prisms |
| seams::profile | seams_brine_profile |
| seams::q4q6 | seams_q4q6, seams_system_count |
| seams::q6bar | seams_q6bar, seams_system_count |
| seams::ql | seams_ql, seams_system_count |
| seams::rdf | seams_rdf |
| seams::read | seams_system_read |
| seams::rejection | seams_brine_rejection |
| seams::rings | seams_rings_find |
| seams::running_cn | seams_running_cn |
| seams::seams_brine | seams_brine_analyse, seams_brine_destroy, seams_brine_fronts, seams_brine_pockets, seams_brine_profile, seams_brine_rejection |
| seams::seams_cages | seams_cages, seams_cages_counts, seams_cages_destroy, seams_cages_labels |
| seams::seams_frame | no direct C call |
| seams::seams_neighbours | seams_neighbours_build, seams_neighbours_copy, seams_neighbours_count, seams_neighbours_csr, seams_neighbours_destroy, seams_neighbours_knn |
| seams::seams_rings | seams_rings_count, seams_rings_destroy, seams_rings_find, seams_rings_fraction, seams_rings_histogram, seams_rings_labels, seams_rings_vertices |
| seams::seams_topo | no direct C call |
| seams::set_frame | seams_system_set_frame |
| seams::set_ortho | seams_system_set_ortho |
| seams::signature | seams_signature_cages |
| seams::site_kinds | seams_site_table, seams_system_count |
| seams::sizes | seams_rings_count, seams_rings_labels |
| seams::steinhardt | seams_steinhardt, seams_system_count |
| seams::symmetry | seams_system_count, seams_system_symmetry |
| seams::vertices | seams_rings_vertices |
