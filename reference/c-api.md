<!-- https://docs.dseams.info/reference/c-api.md -->

# C {#c}

`include/seams.h` is the native C API and the ABI the other languages call. The pages under [C API](../capi/index.md) group the same declarations. This page lists every `SEAMS_API` function the header declares, with the signature and the parameter notes from that header.

The C column is the function itself: this is the C API.

## Native API {#native-api}

Generated from the sources named below. A name on this page is a public call in that language.

### seams_abi {#seams-abi}

```c
int seams_abi(void)
```

ABI number compiled into this library.

Returns

SEAMS_ABI, which is 1.

### seams_brine_analyse {#seams-brine-analyse}

```c
seams_status seams_brine_analyse(seams_system *sys, const int *ice, size_t n_ice, int water_type, int ion_type, int n_bins, int axis, seams_brine **out)
```

Planar ice/brine front from an ice assignment.

Parameters

- **sys**: Frame handle.
- **ice**: Per-atom ice flags, length n_ice.
- **n_ice**: Must equal the frame count.
- **water_type**: Species kept as water, or 0.
- **ion_type**: Species kept as ions, or 0 for none.
- **n_bins**: Bin count along axis, or 0 for 20.
- **axis**: 0, 1, or 2.
- **out**: Receives a new handle. On failure *out is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_brine_destroy {#seams-brine-destroy}

```c
void seams_brine_destroy(seams_brine *front)
```

Free a brine-front handle.

Parameters

- **front**: Handle from seams_brine_analyse(), or NULL. NULL is ignored.

### seams_brine_fronts {#seams-brine-fronts}

```c
seams_status seams_brine_fronts(const seams_brine *front, double *z_enter, double *z_exit, double *thickness, int *found)
```

The two planar crossings of the ice slab.

Parameters

- **front**: Brine handle.
- **z_enter**: Where the ice fraction crosses the cut on entering the ice along the positive axis, or NULL.
- **z_exit**: The crossing on leaving the ice, or NULL.
- **thickness**: Ice run along the axis, or NULL.
- **found**: Receives 1 when both crossings exist, otherwise 0, or NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or every output pointer is NULL.

### seams_brine_pockets {#seams-brine-pockets}

```c
seams_status seams_brine_pockets(const seams_brine *front, int *n_pockets, int *reservoir_waters)
```

Brine pockets trapped in the ice, excluding the liquid reservoir.

Parameters

- **front**: Brine handle.
- **n_pockets**: Receives the pocket count, or NULL.
- **reservoir_waters**: Receives the water count of the liquid reservoir, or NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or both output pointers are NULL.

### seams_brine_profile {#seams-brine-profile}

```c
seams_status seams_brine_profile(const seams_brine *front, double *z_mid, double *molality, double *ice_fraction, size_t cap, size_t *n_out)
```

Concentration profile along the front axis.

Parameters

- **front**: Brine handle.
- **z_mid**: Bin centers, or NULL to skip the column.
- **molality**: Bin molalities, or NULL to skip the column.
- **ice_fraction**: Bin ice fractions, or NULL to skip the column.
- **cap**: Capacity of each non-NULL column.
- **n_out**: Receives the bin count.

Returns

::SEAMS_OK when every non-NULL column was filled, or when every column pointer was NULL and n_out was written. ::SEAMS_ERR_RANGE when a non-NULL column is short. ::SEAMS_ERR_NULL when front or n_out is NULL.

### seams_brine_rejection {#seams-brine-rejection}

```c
seams_status seams_brine_rejection(const seams_brine *front, double *c_ice, double *c_liquid, double *k, int *defined)
```

Molal rejection coefficient k = c_ice / c_liquid.

Parameters

- **front**: Brine handle.
- **c_ice**: Bulk ice molality, or NULL.
- **c_liquid**: Bulk liquid molality, or NULL.
- **k**: The ratio, or NULL. Unchanged when the ratio is undefined.
- **defined**: Receives 1 when the liquid bulk contains salt or the ice bulk does not, otherwise 0, or NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or every output pointer is NULL.

### seams_cages {#seams-cages}

```c
seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)
```

Hexagonal and double-diamond cages. Same search as seams_cages_find().

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **graph**: ::SEAMS_GRAPH_CUTOFF, ::SEAMS_GRAPH_MUTUAL, or ::SEAMS_GRAPH_SEEDED.
- **k**: Neighbour count for the seeded graph. Zero selects 4.
- **complete**: Non-zero runs ring-adjacent completion on seeded labels.
- **out**: Receives a ::seams_cage_set. On failure *out is NULL.

Returns

Same codes as seams_cages_find().

### seams_cages_counts {#seams-cages-counts}

```c
seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex, int *n_cubic, int *n_water, int *n_both)
```

Disjoint cage population counts.

Parameters

- **cages**: Cage handle.
- **n_hex**: Hexagonal only, or NULL.
- **n_cubic**: Double-diamond only, or NULL.
- **n_water**: Neither cage, or NULL.
- **n_both**: Both cages, or NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL when cages is NULL or every count pointer is NULL.

### seams_cages_destroy {#seams-cages-destroy}

```c
void seams_cages_destroy(seams_cage_set *cages)
```

Free a cage-membership handle.

Parameters

- **cages**: Handle from seams_cages_find() or seams_cages(), or NULL. NULL is ignored.

### seams_cages_find {#seams-cages-find}

```c
seams_status seams_cages_find(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)
```

Hexagonal and double-diamond cage membership.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom. The seeded graph uses cutoff + 1.5 as the candidate shell.
- **type**: Species filter.
- **graph**: ::SEAMS_GRAPH_CUTOFF, ::SEAMS_GRAPH_MUTUAL, or ::SEAMS_GRAPH_SEEDED.
- **k**: Neighbour count for the seeded graph. Zero selects 4.
- **complete**: Non-zero runs ring-adjacent completion on the seeded labels. Ignored for the cutoff and mutual graphs.
- **out**: Receives a new handle. On failure *out is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_cages_labels {#seams-cages-labels}

```c
seams_status seams_cages_labels(const seams_cage_set *cages, int *per_atom, size_t len)
```

Per-atom cage code.

Parameters

- **cages**: Cage handle.
- **per_atom**: Caller buffer of ::seams_cage_label values, one per atom in the frame. Atoms the species filter dropped are ::SEAMS_CAGE_WATER.
- **len**: Length of per_atom. Must be at least the frame count stored on the handle.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

### seams_cages_signature {#seams-cages-signature}

```c
seams_status seams_cages_signature(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)
```

Count closed cages whose ring census matches a signature.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **spec**: Signature text. Must be non-NULL.
- **n_cages**: Receives the cage count.
- **n_atoms**: Receives the vertex count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_chill {#seams-chill}

```c
seams_status seams_chill(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, size_t labels_len, int counts9[9])
```

Original CHILL labels, published or OVITO edges.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter. Zero and negative keep every atom.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **rule**: ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.
- **labels**: Caller buffer of labels_len ints.
- **labels_len**: Must be at least the frame count.
- **counts9**: Nine ints that receive the population of each label, or NULL to skip the histogram.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_chill_plus {#seams-chill-plus}

```c
seams_status seams_chill_plus(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, int counts9[9])
```

CHILL+ labels, without a buffer length.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter. Zero and negative keep every atom.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **rule**: ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.
- **labels**: Caller buffer of one label per atom.
- **counts9**: Nine ints, or NULL.

Returns

Same codes as seams_chill_plus_n().

### seams_chill_plus_n {#seams-chill-plus-n}

```c
seams_status seams_chill_plus_n(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, size_t labels_len, int counts9[9])
```

CHILL+ labels into a buffer of known length.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter. Zero and negative keep every atom.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **rule**: ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.
- **labels**: Caller buffer of labels_len ints.
- **labels_len**: Must be at least the frame count.
- **counts9**: Nine ints that receive the population of each label, or NULL.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys or labels; ::SEAMS_ERR_STATE when sys holds no frame; ::SEAMS_ERR_RANGE for a short labels buffer or a non-positive or non-finite cutoff; ::SEAMS_ERR_ARG for an unknown rule or graph; ::SEAMS_ERR_NOMEM; ::SEAMS_ERR_INTERNAL.

### seams_contact_pairs {#seams-contact-pairs}

```c
seams_status seams_contact_pairs(seams_system *sys, const char *site_spec, size_t *n_pairs, int *n_cation, int *n_anion)
```

Mutual nearest unlike contacts of a site table.

Parameters

- **sys**: Frame handle.
- **site_spec**: Text such as “1=cationHead,2=anion”. Must be non-NULL.
- **n_pairs**: Receives the contact count, or NULL.
- **n_cation**: Receives the cation molecule count, or NULL.
- **n_anion**: Receives the anion molecule count, or NULL. At least one of the three outputs must be non-NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_coordination {#seams-coordination}

```c
seams_status seams_coordination(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *cn)
```

Coordination number integrated to rmax.

Parameters

- **sys**: Frame handle.
- **rmax**: Positive finite outer radius, Angstrom.
- **bins**: Bin count. Zero selects about rmax / 0.1.
- **type_i**: Central species, as in seams_rdf().
- **type_j**: Counted species, as in seams_rdf().
- **cn**: Receives the coordination number.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_density_axis {#seams-density-axis}

```c
seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis, double *z, double *rho, size_t cap, size_t *n_out)
```

Number density along one cell axis.

Parameters

- **sys**: Frame handle.
- **type**: Species filter. Zero keeps every atom.
- **bins**: Bin count. Zero selects about one bin per 0.1 Angstrom.
- **axis**: 0, 1, or 2.
- **z**: Bin centers, or NULL to query.
- **rho**: Density per bin, or NULL to query.
- **cap**: Capacity of each non-NULL column.
- **n_out**: Receives the bin count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_density_fractional {#seams-density-fractional}

```c
seams_status seams_density_fractional(seams_system *sys, int type, int bins, int axis, const double *direction, size_t direction_len, double origin, const double *masses, size_t masses_len, double *coordinate, double *number_density, double *mass_density, double *bin_volume, size_t cap, size_t *n_out)
```

Number and mass density along a cell axis or a Cartesian direction.

Parameters

- **sys**: Frame handle.
- **type**: Species filter. Zero keeps every atom.
- **bins**: Positive bin count.
- **axis**: 0, 1, or 2. Ignored when direction is not NULL.
- **direction**: Three Cartesian components, or NULL for a cell axis.
- **direction_len**: 3 when direction is not NULL. Ignored otherwise.
- **origin**: Fractional origin. Any finite value is wrapped into \[0, 1).
- **masses**: One mass per stored atom, or NULL for mass 1.
- **masses_len**: Length of masses. Must equal the atom count when masses is not NULL.
- **coordinate**: Bin centres, or NULL to skip.
- **number_density**: Number density, or NULL to skip.
- **mass_density**: Mass density in the caller’s mass units per cubic Angstrom, or NULL to skip. An amu mass gives amu per cubic Angstrom; times 1.66054 that is grams per cubic centimetre.
- **bin_volume**: Bin volume in cubic Angstrom, or NULL to skip.
- **cap**: Capacity of each non-NULL column.
- **n_out**: Receives the bin count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_directed_hbonds {#seams-directed-hbonds}

```c
seams_status seams_directed_hbonds( seams_system *sys, double pair_cutoff, double dist_cutoff, double angle_cutoff, const int *donor_types, size_t n_donor_types, const int *acceptor_types, size_t n_acceptor_types, int hydrogen_type, int exclude_same_mol, seams_directed_hbond *edges, size_t cap, size_t *n_out)
```

Directed hydrogen-bond edges for selected donor and acceptor types.

Parameters

- **sys**: Frame handle.
- **pair_cutoff**: Positive finite donor-acceptor cutoff, Angstrom.
- **dist_cutoff**: Positive finite acceptor-hydrogen cutoff, Angstrom.
- **angle_cutoff**: Positive finite acceptor-first angle, degrees.
- **donor_types**: Species that may donate. NULL only when the count is 0.
- **n_donor_types**: Length of donor_types. Zero writes no edges.
- **acceptor_types**: Species that may accept. NULL only when the count is 0.
- **n_acceptor_types**: Length of acceptor_types. Zero writes no edges.
- **hydrogen_type**: Positive hydrogen species, distinct from both type sets.
- **exclude_same_mol**: 1 drops a shared molID, 0 keeps it.
- **edges**: Caller buffer, or NULL to query the count.
- **cap**: Capacity of edges. Ignored when edges is NULL.
- **n_out**: Receives the edge count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_domains {#seams-domains}

```c
seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec, int subset, int *n_subset, int *largest, double *percolation)
```

Largest domain of a polar or apolar site subset.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cluster cutoff, Angstrom.
- **site_spec**: Site table text. Must be non-NULL.
- **subset**: 0 is apolar, 1 is polar.
- **n_subset**: Receives the number of sites in the subset, or NULL.
- **largest**: Receives the size of the largest domain, or NULL.
- **percolation**: Receives the largest domain size over the subset size, or NULL. At least one output must be non-NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_fingerprint {#seams-fingerprint}

```c
seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k, int graph, int hops, int color_types, seams_topo **out)
```

Topological fingerprint of the bonded graph.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **k**: Neighbour count for the mutual graph. Zero selects 4.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **hops**: Neighbourhood depth. Zero selects 2.
- **color_types**: Non-zero colors vertices by LAMMPS type.
- **out**: Receives a ::seams_topo. On failure *out is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_fingerprint_atom {#seams-fingerprint-atom}

```c
seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf, size_t cap, size_t *n_out)
```

Local key of one atom. Same buffer rule as seams_fingerprint_key().

Parameters

- **fp**: Fingerprint handle.
- **index**: Atom index in the frame.
- **buf**: Caller buffer of cap bytes, or NULL to query.
- **cap**: Capacity of buf, including the terminating NUL.
- **n_out**: Receives the key length without the terminating NUL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE for a bad index or a short buffer.

### seams_fingerprint_destroy {#seams-fingerprint-destroy}

```c
void seams_fingerprint_destroy(seams_topo *fp)
```

Free a fingerprint.

Parameters

- **fp**: Handle from seams_fingerprint(), or NULL. NULL is ignored.

### seams_fingerprint_key {#seams-fingerprint-key}

```c
seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap, size_t *n_out)
```

Frame key text.

Parameters

- **fp**: Fingerprint handle.
- **buf**: Caller buffer of cap bytes, or NULL to query.
- **cap**: Capacity of buf, including the terminating NUL.
- **n_out**: Receives the key length without the terminating NUL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when buf is short.

### seams_fingerprint_rings {#seams-fingerprint-rings}

```c
seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts, size_t counts_len)
```

Primitive-ring census stored on the fingerprint.

Parameters

- **fp**: Fingerprint handle.
- **counts**: Caller buffer. counts\[k\] receives the number of rings of size k. Sizes past the buffer are dropped.
- **counts_len**: Length of counts. Must be positive.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when counts_len is 0.

### seams_get_num_threads {#seams-get-num-threads}

```c
int seams_get_num_threads(void)
```

Team width the next analysis call will request.

Returns

A positive count. After seams_set_num_threads(0) this is the saved default, which is omp_get_max_threads() when OpenMP is linked and 1 otherwise.

### seams_handle_error {#seams-handle-error}

```c
const char *seams_handle_error(const void *handle)
```

Detail string stored on a handle.

Parameters

- **handle**: A ::seams_system, ::seams_neighbours, ::seams_rings, ::seams_cage_set, ::seams_topo, or ::seams_brine, or NULL.

Returns

The handle’s string, or the thread-local string when handle is NULL. Never NULL.

### seams_hbond_edges {#seams-hbond-edges}

```c
seams_status seams_hbond_edges(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int fields, seams_hbond_edge *edges, size_t cap, size_t *n_out)
```

Hydrogen-bond edges for one frame.

Parameters

- **sys**: Frame handle.
- **cutoff**: Donor-acceptor cutoff, Angstrom.
- **oxygen_type**: Heavy-atom species filter as in seams_hbonds().
- **hydrogen_type**: Hydrogen species (positive).
- **dist_cutoff**: Acceptor-H distance, Angstrom.
- **angle_cutoff**: O-O-H angle, degrees.
- **fields**: Bitmask: SEAMS_HBOND_EDGE_WITH_H and/or SEAMS_HBOND_EDGE_WITH_GEOMETRY.
- **edges**: Caller buffer with caller-set size on each element, or NULL to query.
- **cap**: Capacity of @p edges. Ignored when @p edges is NULL.
- **n_out**: Receives the edge count.

Returns

::SEAMS_OK, or the same error classes as seams_hbonds(), plus ::SEAMS_ERR_ARG when fields is unknown or size is invalid, and ::SEAMS_ERR_RANGE when cap is smaller than *n_out.

### seams_hbond_tcf_create {#seams-hbond-tcf-create}

```c
seams_status seams_hbond_tcf_create(size_t max_lag, size_t origin_stride, double dt, int directed, seams_hbond_tcf **out)
```

Create a streaming H-bond lifetime accumulator from per-frame edges.

Parameters

- **max_lag**: Maximum lag in frames, or 0 to retain all lags. Above PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated, which is ::SEAMS_ERR_RANGE.
- **origin_stride**: Frame stride between time origins, at least 1.
- **dt**: Positive finite frame spacing.
- **directed**: Non-zero treats (donor, acceptor) as directed; zero canonicalizes pairs.
- **out**: Receives the new handle; set to NULL on failure.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_hbond_tcf_curve {#seams-hbond-tcf-curve}

```c
seams_status seams_hbond_tcf_curve(const seams_hbond_tcf *tcf, double *lag_time, double *c_intermittent, double *s_continuous, size_t *n_origins, size_t cap, size_t *n_lags)
```

Read intermittent C_HB(t) and continuous S_HB(t) lifetime curves.

Parameters

- **tcf**: H-bond TCF handle.
- **lag_time**: Receives lag times k * dt, or NULL.
- **c_intermittent**: Receives intermittent autocorrelation C_HB, or NULL.
- **s_continuous**: Receives unbroken survival probability S_HB, or NULL.
- **n_origins**: Receives origin counts per lag, or NULL.
- **cap**: Capacity of non-NULL output arrays.
- **n_lags**: Receives the available lag count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_NOMEM.

### seams_hbond_tcf_destroy {#seams-hbond-tcf-destroy}

```c
void seams_hbond_tcf_destroy(seams_hbond_tcf *tcf)
```

Free an H-bond lifetime accumulator.

Parameters

- **tcf**: Handle from seams_hbond_tcf_create(), or NULL.

### seams_hbond_tcf_fit {#seams-hbond-tcf-fit}

```c
seams_status seams_hbond_tcf_fit(const seams_hbond_tcf *tcf, double t_start, double t_end, double *tau_c_int, double *tau_c_exp, double *tau_s_int, double *tau_s_exp)
```

Integrated and single-exponential lifetimes for C_HB and S_HB.

Parameters

- **tcf**: H-bond TCF handle with at least 2 lag points.
- **t_start**: Fit window start time, >= 0.
- **t_end**: Fit window end time, > t_start (or both 0 for full range).
- **tau_c_int**: Receives trapezoidal integral of C_HB, or NULL.
- **tau_c_exp**: Receives exponential decay time of C_HB, or NULL.
- **tau_s_int**: Receives trapezoidal integral of S_HB, or NULL.
- **tau_s_exp**: Receives exponential decay time of S_HB, or NULL.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_hbond_tcf_push {#seams-hbond-tcf-push}

```c
seams_status seams_hbond_tcf_push(seams_hbond_tcf *tcf, const int *donors, const int *acceptors, size_t n_edges)
```

Push one frame’s hydrogen-bond edge set into the lifetime accumulator.

Parameters

- **tcf**: H-bond TCF handle.
- **donors**: Donor ids of length n_edges (may be NULL when n_edges is 0).
- **acceptors**: Acceptor ids of length n_edges (may be NULL when n_edges is 0).
- **n_edges**: Number of edges in this frame.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_NOMEM.

### seams_hbonds {#seams-hbonds}

```c
seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int *n_bonds)
```

Hydrogen-bond count on one frame.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite donor-acceptor cutoff, Angstrom.
- **oxygen_type**: Heavy-atom species. Zero keeps every atom that is not hydrogen_type.
- **hydrogen_type**: Positive hydrogen species.
- **dist_cutoff**: Acceptor-H distance, Angstrom.
- **angle_cutoff**: O-O-H angle, degrees.
- **n_bonds**: Receives the number of hydrogen bonds.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_host_guest {#seams-host-guest}

```c
seams_status seams_host_guest(seams_system *sys, int guest_type, int framework_type, double cutoff, int *guest_index, int *nearest_index, int *nearest_type, double *nearest_distance, size_t guest_cap, size_t *n_guest, int *framework_types, size_t type_cap, size_t *n_types, int *contact_counts, size_t contact_cap)
```

Nearest framework site and contact counts for each guest.

Parameters

- **sys**: Frame handle.
- **guest_type**: Positive guest species.
- **framework_type**: Framework species. Zero keeps every other type. A value equal to guest_type is rejected.
- **cutoff**: Contact distance in Angstrom. Zero counts only a site at the same point. Negative is rejected.
- **guest_index**: Guest indices, or NULL to skip.
- **nearest_index**: Nearest framework index, or -1. NULL skips the column.
- **nearest_type**: Type of that site, or NULL to skip.
- **nearest_distance**: Distance in Angstrom, or NULL to skip.
- **guest_cap**: Capacity of each non-NULL guest column.
- **n_guest**: Receives the guest count.
- **framework_types**: Sorted framework types, or NULL to skip.
- **type_cap**: Capacity of framework_types when that pointer is not NULL.
- **n_types**: Receives the framework type count.
- **contact_counts**: Row-major counts, or NULL to skip.
- **contact_cap**: Capacity of contact_counts. Must hold n_guest * n_types when contact_counts is not NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_label_name {#seams-label-name}

```c
const char *seams_label_name(int label)
```

Name of a CHILL label.

Parameters

- **label**: A ::seams_label value.

Returns

Static string. Unknown codes return “unknown”.

### seams_last_error {#seams-last-error}

```c
const char *seams_last_error(void)
```

Thread-local detail of the last failing call on this thread.

Returns

Pointer to a thread-local string. Never NULL.

### seams_last_error_message {#seams-last-error-message}

```c
const char *seams_last_error_message(void)
```

Thread-local detail of the last failing call on this thread.

Returns

Pointer to a thread-local string. Never NULL.

### seams_lsi {#seams-lsi}

```c
seams_status seams_lsi(seams_system *sys, const int *centres, size_t n_centres, const int *partners, size_t n_partners, double cutoff, double *lsi, size_t lsi_len)
```

Local Structure Index for arbitrary selections.

Parameters

- **sys**: Frame handle.
- **centres**: Centre indices, or NULL when n_centres is zero.
- **n_centres**: Centre count. Zero means every atom.
- **partners**: Partner indices, or NULL when n_partners is zero.
- **n_partners**: Partner count. Zero means every atom.
- **cutoff**: Positive finite shell cutoff.
- **lsi**: Caller buffer, one value per selected centre.
- **lsi_len**: Length of lsi.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_molecule_centres {#seams-molecule-centres}

```c
seams_status seams_molecule_centres(const seams_system *sys, const int *mols, size_t mols_len, const double *masses, size_t masses_len, int *mol_ids_out, double *centres_out, size_t cap, size_t *n_mols)
```

Intact periodic centres of mass or geometry for molecules in a frame.

Parameters

- **sys**: Frame handle.
- **mols**: Optional molecule id per atom, or NULL to use stored molecule ids.
- **mols_len**: Length of mols when mols is non-NULL; must equal n, or 0 when NULL.
- **masses**: Optional positive weight per atom, or NULL for equal weights.
- **masses_len**: Length of masses when masses is non-NULL; must equal n, or 0 when NULL.
- **mol_ids_out**: Receives one molecule id per molecule, or NULL.
- **centres_out**: Receives 3 * n_mols coordinates in Angstrom, or NULL.
- **cap**: Molecule capacity of non-NULL output buffers.
- **n_mols**: Receives the molecule count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.

### seams_molecule_dipoles {#seams-molecule-dipoles}

```c
seams_status seams_molecule_dipoles( const seams_system *sys, const int *mols, size_t mols_len, const double *charges, size_t charges_len, int *mol_ids_out, double *dipoles, double *unit_vectors, size_t cap, size_t *n_mols)
```

Generic per-molecule dipole or geometric orientation vectors.

Parameters

- **sys**: Frame handle.
- **mols**: Optional molecule id per atom, or NULL to use stored molecule ids.
- **mols_len**: Length of mols when non-NULL, or 0 when NULL.
- **charges**: Optional partial charge per atom, or NULL for geometric mode.
- **charges_len**: Length of charges when non-NULL, or 0 when NULL.
- **mol_ids_out**: Receives molecule ids, length cap, or NULL.
- **dipoles**: Receives 3 * cap raw dipole vectors, or NULL.
- **unit_vectors**: Receives 3 * cap normalized orientation vectors, or NULL.
- **cap**: Molecule capacity of non-NULL output arrays.
- **n_mols**: Receives the molecule count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.

### seams_mols {#seams-mols}

```c
seams_status seams_mols(const seams_system *sys, int *mols, size_t len)
```

Copy the stored molecule ids out.

Parameters

- **sys**: Frame handle.
- **mols**: Caller buffer for n ints, in stored order.
- **len**: Length of mols. Must be at least n.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

### seams_monolayer {#seams-monolayer}

```c
seams_status seams_monolayer(seams_system *sys, double cutoff, int type, int max_depth, double sheet_area, int *counts, size_t counts_len)
```

Primitive-ring histogram of a monolayer.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **max_depth**: Largest ring size, at least 3.
- **sheet_area**: Sheet area, square Angstrom. Must be positive.
- **counts**: counts\[k\] is the number of rings of size k.
- **counts_len**: Must be at least max_depth + 1.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_msd_create {#seams-msd-create}

```c
seams_status seams_msd_create(size_t natoms, size_t max_lag, size_t origin_stride, double dt, const double axis[3], seams_msd **out)
```

Create a streaming multi-origin MSD accumulator with axial/lateral split.

Parameters

- **natoms**: Particle or molecule count per frame, at least 1. Above PTRDIFF_MAX / (3 * sizeof(double)) the coordinates cannot be allocated, which is ::SEAMS_ERR_RANGE.
- **max_lag**: Maximum lag in frames, or 0 to retain all lags. Above PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated, which is ::SEAMS_ERR_RANGE.
- **origin_stride**: Frame stride between time origins, at least 1.
- **dt**: Positive finite frame spacing in physical time units.
- **axis**: Non-zero 3-vector defining the axial direction.
- **out**: Receives the new handle; set to NULL on failure.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_msd_curve {#seams-msd-curve}

```c
seams_status seams_msd_curve(const seams_msd *msd, double *lag_time, double *total, double *axial, double *lateral, size_t *n_origins, size_t cap, size_t *n_lags)
```

Read accumulated lag times and total, axial, and lateral MSD curves.

Parameters

- **msd**: MSD handle.
- **lag_time**: Receives lag times k * dt, or NULL.
- **total**: Receives 3D total MSD, or NULL.
- **axial**: Receives 1D axial MSD along axis, or NULL.
- **lateral**: Receives 2D lateral MSD perpendicular to axis, or NULL.
- **n_origins**: Receives time-origin counts per lag, or NULL.
- **cap**: Capacity of non-NULL output arrays.
- **n_lags**: Receives the available lag count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_NOMEM.

### seams_msd_destroy {#seams-msd-destroy}

```c
void seams_msd_destroy(seams_msd *msd)
```

Free an MSD accumulator.

Parameters

- **msd**: Handle from seams_msd_create(), or NULL.

### seams_msd_fit {#seams-msd-fit}

```c
seams_status seams_msd_fit(const seams_msd *msd, double t_start, double t_end, double *d_total, double *d_axial, double *d_lateral, double *alpha)
```

Linear Einstein diffusion fit on an MSD accumulator.

Parameters

- **msd**: MSD handle with at least 2 lag points.
- **t_start**: Window start time, >= 0.
- **t_end**: Window end time, > t_start (or both 0 for default window).
- **d_total**: Receives 3D diffusion coefficient, or NULL.
- **d_axial**: Receives 1D axial diffusion coefficient, or NULL.
- **d_lateral**: Receives 2D lateral diffusion coefficient, or NULL.
- **alpha**: Receives the log-log anomalous exponent, or NULL. Receives NaN when the window holds fewer than two lags with t > 0 and MSD_total > 1e-20, so the log-log fit was skipped; test it with isnan(). The diffusion coefficients are still the linear fit.

Returns

::SEAMS_OK when the linear fit succeeded, including a skipped log-log fit (alpha NaN); ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_msd_push {#seams-msd-push}

```c
seams_status seams_msd_push(seams_msd *msd, const double *xyz, size_t natoms, const double *box9, size_t box_len)
```

Push one frame into a streaming MSD accumulator.

Parameters

- **msd**: MSD handle.
- **xyz**: Coordinates of 3 * natoms doubles in Angstrom.
- **natoms**: Particle count; must match the count passed to seams_msd_create().
- **box9**: Optional 9-number periodic cell, or NULL for unwrapped input.
- **box_len**: Must be 9 when box9 is non-NULL, or 0 when box9 is NULL.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX; ::SEAMS_ERR_ARG; ::SEAMS_ERR_STATE; ::SEAMS_ERR_NOMEM.

### seams_neighbours_build {#seams-neighbours-build}

```c
seams_status seams_neighbours_build(seams_system *sys, double cutoff, int type, int graph, seams_neighbours **out)
```

Build a cutoff or mutual neighbour list.

Parameters

- **sys**: Frame handle that already holds a frame.
- **cutoff**: Positive finite cutoff, Angstrom. The mutual graph keeps the four nearest neighbours inside cutoff + 1.5.
- **type**: Species to keep. Zero and negative keep every atom.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL. ::SEAMS_GRAPH_SEEDED is rejected here.
- **out**: Receives a new handle. On failure *out is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_neighbours_copy {#seams-neighbours-copy}

```c
seams_status seams_neighbours_copy(const seams_neighbours *nb, int *offsets, size_t offsets_cap, size_t *n_offsets, int *cols, size_t cols_cap, size_t *nnz)
```

Copy the neighbour CSR into caller buffers.

Parameters

- **nb**: Neighbour handle.
- **offsets**: Caller buffer of offsets_cap ints, or NULL to query.
- **offsets_cap**: Capacity of offsets. Ignored when offsets is NULL.
- **n_offsets**: Receives the required offset count (n + 1).
- **cols**: Caller buffer of cols_cap ints, or NULL to query.
- **cols_cap**: Capacity of cols. Ignored when cols is NULL.
- **nnz**: Receives the required column count.

Returns

::SEAMS_OK when every non-NULL buffer was filled, or when both buffers were NULL and the lengths were written. ::SEAMS_ERR_RANGE when a non-NULL buffer is short. ::SEAMS_ERR_NULL when nb, n_offsets, or nnz is NULL.

### seams_neighbours_count {#seams-neighbours-count}

```c
int seams_neighbours_count(const seams_neighbours *nb)
```

Number of neighbour rows, as an int.

Parameters

- **nb**: Neighbour handle, or NULL.

Returns

The row count, or 0 when nb is NULL.

### seams_neighbours_csr {#seams-neighbours-csr}

```c
seams_status seams_neighbours_csr(const seams_neighbours *nb, const int **offsets, int *n_offsets, const int **cols, int *nnz)
```

Borrow the neighbour list as CSR.

Parameters

- **nb**: Neighbour handle.
- **offsets**: Receives the offset pointer.
- **n_offsets**: Receives the offset count, n + 1.
- **cols**: Receives the column pointer. May be NULL when the list has no pairs; *nnz is then 0.
- **nnz**: Receives the column count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL.

### seams_neighbours_destroy {#seams-neighbours-destroy}

```c
void seams_neighbours_destroy(seams_neighbours *nb)
```

Free a neighbour list.

Parameters

- **nb**: Handle from seams_neighbours_build() or seams_neighbours_knn(), or NULL. NULL is ignored.

### seams_neighbours_knn {#seams-neighbours-knn}

```c
seams_status seams_neighbours_knn(seams_system *sys, int k, double cutoff, int type, int mutual, seams_neighbours **out)
```

k-nearest neighbour list.

Parameters

- **sys**: Frame handle.
- **k**: Positive neighbour count.
- **cutoff**: Nominal positive cutoff, Angstrom. Candidates are taken inside cutoff + 1.5.
- **type**: Species filter. Zero and negative keep every atom.
- **mutual**: Non-zero keeps the mutual graph. Zero keeps the one-way k-nearest graph.
- **out**: Receives a neighbour handle. On failure *out is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_neighbours_nrows {#seams-neighbours-nrows}

```c
seams_status seams_neighbours_nrows(const seams_neighbours *nb, size_t *n)
```

Number of neighbour rows.

Parameters

- **nb**: Neighbour handle.
- **n**: Receives the row count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL when nb or n is NULL.

### seams_orientation_order {#seams-orientation-order}

```c
seams_status seams_orientation_order(const double *vectors, size_t n_vectors, const double axis[3], double *cos_theta, double *p2, size_t cap, double *mean_p1, double *mean_p2)
```

Project molecular vectors onto an axis for P1 and P2 order parameters.

Parameters

- **vectors**: Input 3-vectors, 3 * n_vectors doubles (or NULL when n_vectors is 0).
- **n_vectors**: Number of 3-vectors.
- **axis**: Non-zero reference axis of 3 finite doubles.
- **cos_theta**: Optional output buffer of cap doubles for cos(theta), or NULL.
- **p2**: Optional output buffer of cap doubles for P2(cos(theta)), or NULL.
- **cap**: Capacity of non-NULL per-vector buffers; must be at least n_vectors.
- **mean_p1**: Receives ensemble mean of cos(theta), or NULL.
- **mean_p2**: Receives ensemble mean of P2(cos(theta)), or NULL.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_path_nframes {#seams-path-nframes}

```c
seams_status seams_path_nframes(const char *path, size_t *nframes)
```

Number of frames in a trajectory path.

Parameters

- **path**: Filesystem path.
- **nframes**: Receives the count. Must be non-NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_INTERNAL. *nframes is unchanged on failure.

### seams_pore_frame {#seams-pore-frame}

```c
seams_status seams_pore_frame(seams_system *sys, int framework_type, int axis, const double *direction, size_t direction_len, const double *centre_frac, size_t centre_len, double *axis_out, double *centre_out, double *centre_frac_out, int *axis_index)
```

Pore axis and centre, supplied or taken from the framework.

Parameters

- **sys**: Frame handle.
- **framework_type**: Species used to place an automatic centre. Zero keeps every atom.
- **axis**: 0, 1, 2, or -1 for the longest edge when direction is NULL. Any other negative value is rejected.
- **direction**: Three Cartesian components parallel to a cell edge, or NULL.
- **direction_len**: 3 when direction is not NULL. Ignored otherwise.
- **centre_frac**: Three fractional coordinates, or NULL to detect the centre.
- **centre_len**: 3 when centre_frac is not NULL. Ignored otherwise.
- **axis_out**: Receives the unit Cartesian axis. Three doubles.
- **centre_out**: Receives the Cartesian centre. Three doubles.
- **centre_frac_out**: Receives the fractional centre. Three doubles.
- **axis_index**: Receives 0, 1, or 2.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_pore_profiles {#seams-pore-profiles}

```c
seams_status seams_pore_profiles( seams_system *sys, int type, int radial_bins, int azimuthal_bins, int axial_bins, double r_max, int axis, const double *direction, size_t direction_len, const double *centre_frac, size_t centre_len, int framework_type, double origin, const double *masses, size_t masses_len, double *radius, double *radial_number, double *radial_mass, double *radial_volume, size_t radial_cap, size_t *n_radial, double *theta, double *azimuthal_number, double *azimuthal_mass, double *azimuthal_volume, size_t azimuthal_cap, size_t *n_azimuthal, double *axial, double *axial_number, double *axial_mass, double *axial_volume, size_t axial_cap, size_t *n_axial)
```

Radial, azimuthal, and axial densities about a pore axis.

Parameters

- **sys**: Frame handle.
- **type**: Species filter. Zero keeps every atom.
- **radial_bins**: Positive radial bin count.
- **azimuthal_bins**: Positive azimuthal bin count.
- **axial_bins**: Positive axial bin count.
- **r_max**: Outer radius in Angstrom. Zero covers the face.
- **axis**: 0, 1, 2, or -1 when direction is NULL. Any other negative value is rejected.
- **direction**: Three Cartesian components parallel to a cell edge, or NULL.
- **direction_len**: 3 when direction is not NULL. Ignored otherwise.
- **centre_frac**: Three fractional coordinates, or NULL to detect the centre.
- **centre_len**: 3 when centre_frac is not NULL. Ignored otherwise.
- **framework_type**: Species for an automatic centre. Zero keeps every atom.
- **origin**: Fractional axial origin, wrapped into \[0, 1).
- **masses**: One mass per stored atom, or NULL for mass 1.
- **masses_len**: Length of masses. Must equal the atom count when masses is not NULL.
- **radius**: Radial bin centres, or NULL to skip.
- **radial_number**: Radial number density, or NULL to skip.
- **radial_mass**: Radial mass density in the caller’s mass units per cubic Angstrom, or NULL to skip. An amu mass gives amu per cubic Angstrom; times 1.66054 that is grams per cubic centimetre.
- **radial_volume**: Radial shell volume, or NULL to skip.
- **radial_cap**: Capacity of each non-NULL radial column.
- **n_radial**: Receives the radial bin count.
- **theta**: Azimuthal bin centres in radians, or NULL to skip.
- **azimuthal_number**: Azimuthal number density, or NULL to skip.
- **azimuthal_mass**: Azimuthal mass density in the caller’s mass units per cubic Angstrom, or NULL to skip. An amu mass gives amu per cubic Angstrom; times 1.66054 that is grams per cubic centimetre.
- **azimuthal_volume**: Azimuthal sector volume, or NULL to skip.
- **azimuthal_cap**: Capacity of each non-NULL azimuthal column.
- **n_azimuthal**: Receives the azimuthal bin count.
- **axial**: Axial bin centres, or NULL to skip.
- **axial_number**: Axial number density, or NULL to skip.
- **axial_mass**: Axial mass density in the caller’s mass units per cubic Angstrom, or NULL to skip. An amu mass gives amu per cubic Angstrom; times 1.66054 that is grams per cubic centimetre.
- **axial_volume**: Axial slab volume, or NULL to skip.
- **axial_cap**: Capacity of each non-NULL axial column.
- **n_axial**: Receives the axial bin count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_prisms {#seams-prisms}

```c
seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect, int *imperfect, size_t len)
```

Perfect and imperfect prism counts by ring size.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **perfect**: Caller buffer. perfect\[k\] receives the perfect k-prisms.
- **imperfect**: Caller buffer. imperfect\[k\] receives the deformed k-prisms.
- **len**: Length of each buffer. At least 7, so that sizes 3 through 6 are filled.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_q4q6 {#seams-q4q6}

```c
seams_status seams_q4q6(seams_system *sys, double cutoff, int type, double *q4, double *q6)
```

Local q4 and q6 on the cutoff shell, without buffer lengths.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **q4**: Local q4, or NULL.
- **q6**: Local q6, or NULL. q4 and q6 may not both be NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_q4q6_n {#seams-q4q6-n}

```c
seams_status seams_q4q6_n(seams_system *sys, double cutoff, int type, double *q4, size_t q4_len, double *q6, size_t q6_len)
```

Local q4 and q6 on the cutoff shell into buffers of known length.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **q4**: Local q4, or NULL.
- **q4_len**: Length of q4. Ignored when q4 is NULL.
- **q6**: Local q6, or NULL. q4 and q6 may not both be NULL.
- **q6_len**: Length of q6. Ignored when q6 is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE for a short buffer or a bad cutoff, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_q6bar {#seams-q6bar}

```c
seams_status seams_q6bar(seams_system *sys, double cutoff, int type, double ice_cut, double *q6bar, size_t q6bar_len, double *mean, double *median, int *n_ice)
```

Lechner-Dellago q6-bar and the ice count above a cut.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **ice_cut**: Finite threshold on q6-bar.
- **q6bar**: Per-atom average, or NULL to keep only the summary.
- **q6bar_len**: Length of q6bar. Ignored when q6bar is NULL.
- **mean**: Receives the mean, or NULL.
- **median**: Receives the median, or NULL.
- **n_ice**: Receives the count above ice_cut, or NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL. At least one of q6bar, mean, median, and n_ice must be non-NULL.

### seams_ql {#seams-ql}

```c
seams_status seams_ql(seams_system *sys, double cutoff, int type, int order_l, int graph, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len)
```

Local Steinhardt q_l, and the Lechner-Dellago average on request.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter. Zero and negative keep every atom.
- **order_l**: Spherical-harmonic degree.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **ql**: Local q_l, or NULL.
- **ql_len**: Length of ql. Ignored when ql is NULL.
- **qlbar**: Averaged q_l, or NULL to skip the average.
- **qlbar_len**: Length of qlbar. Ignored when qlbar is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_quenched_analysis {#seams-quenched-analysis}

```c
seams_status seams_quenched_analysis( seams_system *sys, seams_quench_callback callback, void *user_data, double cutoff, int type, int graph, int max_depth, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int *raw_labels, size_t raw_labels_len, int *quenched_labels, size_t quenched_labels_len, int counts_raw[9], int counts_quenched[9], int *rings_raw, int *rings_quenched, size_t ring_counts_len, int *hbond_raw, int *hbond_quenched, int *hbond_edge_delta, int *ring_count_delta, int *class_change_count)
```

Run CHILL+, ring census, and H-bonds on raw and quenched frames.

Parameters

- **sys**: Frame handle.
- **callback**: Quench callback that writes output coordinates.
- **user_data**: Opaque callback payload.
- **cutoff**: Positive finite neighbour cutoff, Angstrom.
- **type**: Species filter for neighbour, CHILL+, and rings.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **max_depth**: Primitive-ring depth, from 3 through 12.
- **oxygen_type**: Heavy-atom species filter for H-bonds.
- **hydrogen_type**: Hydrogen species filter for H-bonds.
- **dist_cutoff**: Acceptor-H distance, Angstrom.
- **angle_cutoff**: O-O-H angle, degrees.
- **raw_labels**: Raw-frame CHILL+ labels, length at least atom count.
- **raw_labels_len**: Length of @p raw_labels.
- **quenched_labels**: Quenched-frame CHILL+ labels, length at least atom count.
- **quenched_labels_len**: Length of @p quenched_labels.
- **counts_raw**: Raw-frame CHILL+ population counts.
- **counts_quenched**: Quenched-frame CHILL+ population counts.
- **rings_raw**: Raw-frame primitive-ring histogram.
- **rings_quenched**: Quenched-frame primitive-ring histogram.
- **ring_counts_len**: Length of @p rings_raw and @p rings_quenched.
- **hbond_raw**: Receives raw-frame H-bond edge count.
- **hbond_quenched**: Receives quenched-frame H-bond edge count.
- **hbond_edge_delta**: Receives symmetric-difference edge count between raw and quenched H-bond edge sets.
- **ring_count_delta**: Receives sum over ring sizes of absolute raw-minus- quenched count differences.
- **class_change_count**: Receives the atom count whose CHILL+ class changed.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_rdf {#seams-rdf}

```c
seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *g, int *count, size_t cap, size_t *n_out, double *volume)
```

Partial radial distribution.

Parameters

- **sys**: Frame handle.
- **rmax**: Positive finite outer radius, Angstrom.
- **bins**: Bin count. Zero selects about rmax / 0.1.
- **type_i**: First species. Zero and negative keep every atom of a single-species frame, or every atom as one species otherwise.
- **type_j**: Second species. Zero and negative reuse the type_i species.
- **r**: Bin centers, or NULL to query.
- **g**: g(r) per bin, or NULL to query.
- **count**: Pair count per bin, or NULL to query.
- **cap**: Capacity of each non-NULL column.
- **n_out**: Receives the bin count.
- **volume**: Receives the cell volume, or NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_reorient_create {#seams-reorient-create}

```c
seams_status seams_reorient_create(size_t n_vectors, size_t max_lag, size_t origin_stride, double dt, seams_reorient **out)
```

Create a streaming orientational correlation accumulator for C1 and C2.

Parameters

- **n_vectors**: Number of molecular vectors per frame, at least 1. Above PTRDIFF_MAX / (3 * sizeof(double)) the vectors cannot be allocated, which is ::SEAMS_ERR_RANGE.
- **max_lag**: Maximum lag in frames, or 0 to retain all lags. Above PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated, which is ::SEAMS_ERR_RANGE.
- **origin_stride**: Frame stride between time origins, at least 1.
- **dt**: Positive finite frame spacing.
- **out**: Receives the new handle; set to NULL on failure.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_reorient_curve {#seams-reorient-curve}

```c
seams_status seams_reorient_curve(const seams_reorient *reorient, double *lag_time, double *c1, double *c2, size_t *n_origins, size_t cap, size_t *n_lags)
```

Read accumulated C1(t) = <P1(u(0).u(t))> and C2(t) = <P2(u(0).u(t))>.

Parameters

- **reorient**: Reorientation handle.
- **lag_time**: Receives lag times k * dt, or NULL.
- **c1**: Receives first Legendre correlation C1, or NULL.
- **c2**: Receives second Legendre correlation C2, or NULL.
- **n_origins**: Receives origin counts per lag, or NULL.
- **cap**: Capacity of non-NULL output arrays.
- **n_lags**: Receives the available lag count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_NOMEM.

### seams_reorient_destroy {#seams-reorient-destroy}

```c
void seams_reorient_destroy(seams_reorient *reorient)
```

Free a reorientation accumulator.

Parameters

- **reorient**: Handle from seams_reorient_create(), or NULL.

### seams_reorient_fit {#seams-reorient-fit}

```c
seams_status seams_reorient_fit(const seams_reorient *reorient, double t_start, double t_end, double *tau1_int, double *tau1_exp, double *tau2_int, double *tau2_exp)
```

Integrated and single-exponential reorientation times tau_1 and tau_2.

Parameters

- **reorient**: Reorientation handle with at least 2 lag points.
- **t_start**: Fit window start time, >= 0.
- **t_end**: Fit window end time, > t_start (or both 0 for full range).
- **tau1_int**: Receives trapezoidal integral of C1, or NULL.
- **tau1_exp**: Receives exponential decay time of C1, or NULL.
- **tau2_int**: Receives trapezoidal integral of C2, or NULL.
- **tau2_exp**: Receives exponential decay time of C2, or NULL.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_reorient_push {#seams-reorient-push}

```c
seams_status seams_reorient_push(seams_reorient *reorient, const double *vectors, size_t n_vectors)
```

Push one frame of molecular 3-vectors into a reorientation accumulator.

Parameters

- **reorient**: Reorientation handle.
- **vectors**: Array of 3 * n_vectors finite doubles.
- **n_vectors**: Vector count; must match seams_reorient_create().

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_ring_symmetry_name {#seams-ring-symmetry-name}

```c
const char *seams_ring_symmetry_name(int symmetry)
```

Name of a ring-symmetry code.

Parameters

- **symmetry**: A ::seams_ring_symmetry value.

Returns

Static string. Unknown codes return “unknown”.

### seams_rings_chordless {#seams-rings-chordless}

```c
seams_status seams_rings_chordless(const seams_neighbours *nb, int size, int *rings, size_t rings_len, size_t *n_rings, int *per_atom, size_t per_atom_len)
```

Chordless cycles of one size.

Parameters

- **nb**: Graph from seams_neighbours_build().
- **size**: Cycle length, from 3 through 12.
- **rings**: Caller buffer of size times the cycle count ints, or NULL to skip the list.
- **rings_len**: Length of rings in ints. Ignored when rings is NULL.
- **n_rings**: Receives the cycle count. Required when rings is not NULL.
- **per_atom**: Caller buffer of one int per atom, or NULL. Receives the number of cycles through each atom.
- **per_atom_len**: Length of per_atom. Ignored when per_atom is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL. rings, n_rings and per_atom may not all be NULL. A short rings buffer returns ::SEAMS_ERR_RANGE after n_rings and per_atom are written.

### seams_rings_count {#seams-rings-count}

```c
int seams_rings_count(const seams_rings *rings)
```

Number of primitive rings, as an int.

Parameters

- **rings**: Ring handle, or NULL.

Returns

The ring count, or 0 when rings is NULL.

### seams_rings_destroy {#seams-rings-destroy}

```c
void seams_rings_destroy(seams_rings *rings)
```

Free a ring set.

Parameters

- **rings**: Handle from seams_rings_find(), or NULL. NULL is ignored.

### seams_rings_find {#seams-rings-find}

```c
seams_status seams_rings_find(seams_system *sys, double cutoff, int type, int graph, int max_depth, seams_rings **out)
```

Primitive rings (Franzblau) up to a maximum size.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **max_depth**: Largest ring size, from 3 through 12.
- **out**: Receives a new handle. On failure *out is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_rings_fraction {#seams-rings-fraction}

```c
seams_status seams_rings_fraction(const seams_rings *rings, int size, double *fraction)
```

Fraction of atoms that sit on a ring.

Parameters

- **rings**: Ring handle.
- **size**: Ring size to count. Zero counts atoms that sit on any ring.
- **fraction**: Receives a value in \[0, 1\].

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_NOMEM.

### seams_rings_from_csr {#seams-rings-from-csr}

```c
seams_status seams_rings_from_csr( const int *offsets, size_t offsets_len, const int *columns, size_t columns_len, int max_depth, const seams_system *geometry, seams_rings **out)
```

Primitive rings from an explicit caller-supplied CSR graph.

Parameters

- **offsets**: CSR offsets.
- **offsets_len**: Number of offsets; at least one.
- **columns**: CSR columns, or NULL when columns_len is zero.
- **columns_len**: Number of columns.
- **max_depth**: Largest ring size, from 3 through 12.
- **geometry**: Optional frame supplying coordinates and periodic cell.
- **out**: Receives a new ring handle.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_rings_histogram {#seams-rings-histogram}

```c
seams_status seams_rings_histogram(const seams_rings *rings, int *counts, int counts_len)
```

Count primitive rings by size, with an int length.

Parameters

- **rings**: Ring handle.
- **counts**: Caller buffer. counts\[k\] receives the number of rings of size k. Sizes outside the buffer are ignored.
- **counts_len**: Length of counts. Must be positive.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

### seams_rings_histogram_n {#seams-rings-histogram-n}

```c
seams_status seams_rings_histogram_n(const seams_rings *rings, int *counts, size_t counts_len)
```

Count primitive rings by size into a buffer of known length.

Parameters

- **rings**: Ring handle.
- **counts**: Caller buffer. counts\[k\] receives the number of rings of size k. Sizes outside the buffer are ignored.
- **counts_len**: Length of counts. Must be positive.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when counts_len is 0.

### seams_rings_labels {#seams-rings-labels}

```c
seams_status seams_rings_labels(const seams_rings *rings, int *sizes, int *symmetry)
```

Size and symmetry of each ring, without buffer lengths.

Parameters

- **rings**: Ring handle.
- **sizes**: One int per ring, or NULL.
- **symmetry**: One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_rings_labels_n {#seams-rings-labels-n}

```c
seams_status seams_rings_labels_n(const seams_rings *rings, int *sizes, size_t sizes_len, int *symmetry, size_t symmetry_len)
```

Size and symmetry of each ring into buffers of known length.

Parameters

- **rings**: Ring handle.
- **sizes**: One int per ring, or NULL.
- **sizes_len**: Length of sizes. Ignored when sizes is NULL.
- **symmetry**: One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.
- **symmetry_len**: Length of symmetry. Ignored when symmetry is NULL.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when a buffer is shorter than the ring count, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_rings_nrings {#seams-rings-nrings}

```c
seams_status seams_rings_nrings(const seams_rings *rings, size_t *n)
```

Number of primitive rings.

Parameters

- **rings**: Ring handle.
- **n**: Receives the ring count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL when rings or n is NULL.

### seams_rings_vertices {#seams-rings-vertices}

```c
seams_status seams_rings_vertices(const seams_rings *rings, int index, int *atoms, size_t cap, size_t *n_out)
```

Atom indices of one ring.

Parameters

- **rings**: Ring handle.
- **index**: Ring index, from 0 to seams_rings_count(rings) - 1.
- **atoms**: Caller buffer, or NULL to query.
- **cap**: Capacity of atoms. Ignored when atoms is NULL.
- **n_out**: Receives the vertex count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

### seams_rings_winding {#seams-rings-winding}

```c
seams_status seams_rings_winding( const seams_rings *rings, int axis, int *axis_winding, size_t axis_winding_len, int *contractible, size_t contractible_len)
```

Periodic-axis winding and contractibility of every ring.

Parameters

- **rings**: Ring handle.
- **axis**: Periodic cell axis, 0, 1, or 2.
- **axis_winding**: Signed winding output, or NULL.
- **axis_winding_len**: Length of axis_winding.
- **contractible**: One for a zero three-component winding, or NULL.
- **contractible_len**: Length of contractible.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE.

### seams_running_cn {#seams-running-cn}

```c
seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *cn, size_t cap, size_t *n_out)
```

Running coordination number along the partial RDF.

Parameters

- **sys**: Frame handle.
- **rmax**: Positive finite outer radius, Angstrom.
- **bins**: Bin count. Zero selects about rmax / 0.1.
- **type_i**: Central species, as in seams_rdf().
- **type_j**: Counted species, as in seams_rdf().
- **r**: Bin centers, Angstrom, or NULL to query.
- **cn**: Running coordination number, or NULL to query.
- **cap**: Capacity of each non-NULL column.
- **n_out**: Receives the bin count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_set_mols {#seams-set-mols}

```c
seams_status seams_set_mols(seams_system *sys, const int *mols, size_t len)
```

Replace the stored molecule ids on a loaded frame.

Parameters

- **sys**: Frame handle.
- **mols**: Caller buffer of n molecule ids.
- **len**: Length of mols. Must equal the stored atom count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE when len differs from the atom count.

### seams_set_num_threads {#seams-set-num-threads}

```c
seams_status seams_set_num_threads(int n)
```

Set the OpenMP team used by later analysis calls.

Parameters

- **n**: Positive width, or 0 to restore the width captured at first use. A negative width is rejected.

Returns

::SEAMS_OK, or ::SEAMS_ERR_RANGE when n is negative, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_signature_cages {#seams-signature-cages}

```c
seams_status seams_signature_cages(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)
```

Closed cages of one ring census. Same search as seams_cages_signature().

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **spec**: Signature text such as “512” or “hc”. Must be non-NULL.
- **n_cages**: Receives the cage count.
- **n_atoms**: Receives the vertex count.

Returns

Same codes as seams_cages_signature().

### seams_signature_cages_edge {#seams-signature-cages-edge}

```c
seams_status seams_signature_cages_edge(seams_system *sys, double cutoff, int type, int graph, const char *spec, const char *edge, size_t *n_cages, size_t *n_atoms)
```

Closed cages of one ring census with an explicit face-edge window.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **graph**: ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.
- **spec**: Signature text such as “512” or “hc”. Must be non-NULL.
- **edge**: Face-edge window, or NULL. “none” turns the window off.
- **n_cages**: Receives the cage count.
- **n_atoms**: Receives the vertex count.

Returns

Same codes as seams_cages_signature().

### seams_site_table {#seams-site-table}

```c
seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds, size_t len)
```

Per-atom site kind from a site-table spec.

Parameters

- **sys**: Frame handle.
- **spec**: Site table text such as “1=cationHead,2=anion”. Must be non-NULL.
- **kinds**: One ::seams_site_kind per atom.
- **len**: Length of kinds. Must be at least the frame count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_status_string {#seams-status-string}

```c
const char *seams_status_string(seams_status status)
```

Short name of a status code.

Parameters

- **status**: A ::seams_status value.

Returns

Pointer to a static string. Unknown codes return “unknown status”.

### seams_steinhardt {#seams-steinhardt}

```c
seams_status seams_steinhardt(seams_system *sys, double cutoff, int type, int order_l, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len)
```

Local Steinhardt q_l on the cutoff shell, and the average on request.

Parameters

- **sys**: Frame handle.
- **cutoff**: Positive finite cutoff, Angstrom.
- **type**: Species filter.
- **order_l**: 3, 4, 6, 8, or 12.
- **ql**: Local q_l, or NULL.
- **ql_len**: Length of ql. Ignored when ql is NULL.
- **qlbar**: Averaged q_l, or NULL to skip the Lechner-Dellago average.
- **qlbar_len**: Length of qlbar. Ignored when qlbar is NULL.

Returns

Same codes as seams_ql().

### seams_system_box {#seams-system-box}

```c
seams_status seams_system_box(const seams_system *sys, double *box9, size_t box_len)
```

Copy the stored cell out as nine LAMMPS numbers.

Parameters

- **sys**: Frame handle.
- **box9**: Caller buffer: ox, oy, oz, lx, ly, lz, xy, xz, yz. The values match the cell that was stored, to rounding.
- **box_len**: Must be 9.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE when box_len != 9; ::SEAMS_ERR_STATE when no frame is stored.

### seams_system_count {#seams-system-count}

```c
int seams_system_count(const seams_system *sys)
```

Number of points in the frame, as an int.

Parameters

- **sys**: Frame handle, or NULL.

Returns

The atom count, or 0 when sys is NULL or holds no frame.

### seams_system_create {#seams-system-create}

```c
seams_system *seams_system_create(void)
```

Allocate an empty frame.

Returns

A new handle, or NULL on allocation failure. seams_last_error_message() then explains the failure.

### seams_system_destroy {#seams-system-destroy}

```c
void seams_system_destroy(seams_system *sys)
```

Free a frame.

Parameters

- **sys**: Handle from seams_system_create(), or NULL. NULL is ignored.

### seams_system_ids {#seams-system-ids}

```c
seams_status seams_system_ids(const seams_system *sys, int *ids, size_t len)
```

Copy the stored atom ids out.

Parameters

- **sys**: Frame handle.
- **ids**: Caller buffer for n ints, in stored order.
- **len**: Length of ids. Must be at least n.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

### seams_system_natoms {#seams-system-natoms}

```c
seams_status seams_system_natoms(const seams_system *sys, size_t *n)
```

Number of points in the frame.

Parameters

- **sys**: Frame handle.
- **n**: Receives the atom count.

Returns

::SEAMS_OK, ::SEAMS_ERR_NULL when sys or n is NULL, or ::SEAMS_ERR_STATE when sys holds no frame. *n is unchanged on failure.

### seams_system_positions {#seams-system-positions}

```c
seams_status seams_system_positions(const seams_system *sys, double *xyz, size_t len)
```

Copy the stored coordinates out.

Parameters

- **sys**: Frame handle.
- **xyz**: Caller buffer for 3 n doubles, Angstrom, in stored order.
- **len**: Length of xyz in doubles. Must be at least 3 n.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

### seams_system_read {#seams-system-read}

```c
seams_status seams_system_read(seams_system *sys, const char *path, int frame, int type_filter)
```

Read one frame from a trajectory file.

Parameters

- **sys**: Frame handle. On failure the previously stored frame is kept.
- **path**: Filesystem path. CON, LAMMPS dump, and the formats the linked readers accept.
- **frame**: Frame index, 1-based, matching the file order.
- **type_filter**: Negative keeps every atom. Zero guesses the oxygen species. A positive value keeps that species number.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE (frame index or empty frame), ::SEAMS_ERR_ARG (unreadable file), ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_system_set_atoms {#seams-system-set-atoms}

```c
seams_status seams_system_set_atoms(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len, const int *types, const int *ids)
```

Store points and a general cell. Lengths are size_t.

Parameters

- **sys**: Frame handle.
- **xyz**: n times 3 Cartesian coordinates, Angstrom. Every value must be finite.
- **n**: Atom count. Zero is rejected. A count that does not fit in the frame’s integer width is rejected.
- **box9**: Cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.
- **box_len**: Must be 9.
- **types**: Atom types, or NULL to store type 1 for every atom.
- **ids**: Atom ids, or NULL to store 1..n. Ids must be unique.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys, xyz, or box9; ::SEAMS_ERR_RANGE for n == 0 or box_len != 9; ::SEAMS_ERR_BOX for a non-finite cell or a non-positive edge; ::SEAMS_ERR_ARG for a non-finite coordinate or repeated ids; ::SEAMS_ERR_NOMEM; ::SEAMS_ERR_INTERNAL.

### seams_system_set_frame {#seams-system-set-frame}

```c
seams_status seams_system_set_frame(seams_system *sys, const double *xyz, int n, const double *box9, const int *types, const int *ids)
```

Store points and a general LAMMPS cell, with an int count.

Parameters

- **sys**: Frame handle.
- **xyz**: n times 3 coordinates, Angstrom.
- **n**: Atom count. A non-positive count is rejected.
- **box9**: Nine cell numbers: ox, oy, oz, lx, ly, lz, xy, xz, yz.
- **types**: Types, or NULL for type 1.
- **ids**: Ids, or NULL for 1..n.

Returns

Same codes as seams_system_set_atoms().

### seams_system_set_hydrogen_frame {#seams-system-set-hydrogen-frame}

```c
seams_status seams_system_set_hydrogen_frame( seams_system *sys, const double *heavy_xyz, size_t n_heavy, const int *heavy_types, const int *heavy_ids, const int *heavy_mols, const double *hydrogen_xyz, size_t n_hydrogen, const int *hydrogen_types, const int *hydrogen_ids, const int *hydrogen_mols, const double *box9, size_t box_len)
```

Store a heavy-atom cloud and a separate hydrogen cloud.

Parameters

- **sys**: Frame handle.
- **heavy_xyz**: Interleaved heavy coordinates, Angstrom.
- **n_heavy**: Heavy-atom count. Must be positive.
- **heavy_types**: One species per heavy atom.
- **heavy_ids**: One id per heavy atom, or NULL for 1..n_heavy.
- **heavy_mols**: One molID per heavy atom, or NULL to copy the atom id.
- **hydrogen_xyz**: Interleaved hydrogen coordinates, Angstrom.
- **n_hydrogen**: Hydrogen count. Zero allows the hydrogen pointers to be NULL.
- **hydrogen_types**: One species per hydrogen, required when n_hydrogen is positive.
- **hydrogen_ids**: One id per hydrogen, or NULL to continue after the heavies.
- **hydrogen_mols**: One molID per hydrogen, or NULL to copy the atom id.
- **box9**: ox, oy, oz, lx, ly, lz, xy, xz, yz.
- **box_len**: Length of box9. Must be 9.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_system_set_ortho {#seams-system-set-ortho}

```c
seams_status seams_system_set_ortho(seams_system *sys, const double *xyz, int n, double ox, double oy, double oz, double lx, double ly, double lz, const int *types, const int *ids)
```

Store points and an orthorhombic cell.

Parameters

- **sys**: Frame handle.
- **xyz**: n times 3 coordinates, Angstrom.
- **n**: Atom count.
- **ox**: Origin x, Angstrom.
- **oy**: Origin y, Angstrom.
- **oz**: Origin z, Angstrom.
- **lx**: Edge length x, Angstrom.
- **ly**: Edge length y, Angstrom.
- **lz**: Edge length z, Angstrom.
- **types**: Types, or NULL.
- **ids**: Ids, or NULL.

Returns

Same codes as seams_system_set_frame(). The tilt is zero.

### seams_system_set_positions {#seams-system-set-positions}

```c
seams_status seams_system_set_positions(seams_system *sys, const double *xyz, size_t n, const double *box9, size_t box_len)
```

Overwrite the coordinates, and optionally the cell, of the stored frame.

Parameters

- **sys**: Frame handle with a stored frame.
- **xyz**: n times 3 Cartesian coordinates, Angstrom. Every value must be finite.
- **n**: Atom count. Must equal the stored count.
- **box9**: New cell (ox, oy, oz, lx, ly, lz, xy, xz, yz), or NULL to keep the stored cell.
- **box_len**: 9 when box9 is not NULL; ignored otherwise.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys or xyz; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE when n differs from the stored count or box_len != 9; ::SEAMS_ERR_BOX for a non-finite cell or a non-positive edge; ::SEAMS_ERR_ARG for a non-finite coordinate; ::SEAMS_ERR_NOMEM; ::SEAMS_ERR_INTERNAL.

### seams_system_stage_positions {#seams-system-stage-positions}

```c
seams_status seams_system_stage_positions(seams_system *sys, int enable)
```

Trade memory for a single pass in seams_system_set_positions().

Parameters

- **sys**: Frame handle, with or without a stored frame.
- **enable**: Nonzero to keep the second set, 0 to free it.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys.

### seams_system_symmetry {#seams-system-symmetry}

```c
seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd, int *labels, size_t labels_len)
```

Per-atom ice shell label.

Parameters

- **sys**: Frame handle.
- **shell_rmsd**: Positive scale-free residual ceiling, or a non-positive value for the library ceiling.
- **labels**: Caller buffer.
- **labels_len**: Must be at least the frame count.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_INTERNAL, ::SEAMS_ERR_NOMEM.

### seams_system_types {#seams-system-types}

```c
seams_status seams_system_types(const seams_system *sys, int *types, size_t len)
```

Copy the stored atom types out.

Parameters

- **sys**: Frame handle.
- **types**: Caller buffer for n ints, in stored order.
- **len**: Length of types. Must be at least n.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

### seams_tetrahedral_order {#seams-tetrahedral-order}

```c
seams_status seams_tetrahedral_order( seams_system *sys, const int *centres, size_t n_centres, const int *partners, size_t n_partners, int convention, double cutoff, double *q, size_t q_len)
```

Errington-Debenedetti tetrahedral order for arbitrary selections.

Parameters

- **sys**: Frame handle.
- **centres**: Centre indices, or NULL when n_centres is zero.
- **n_centres**: Centre count. Zero means every atom.
- **partners**: Partner indices, or NULL when n_partners is zero.
- **n_partners**: Partner count. Zero means every atom.
- **convention**: One ::seams_tetrahedral_neighbours value. Any other value returns ::SEAMS_ERR_ARG.
- **cutoff**: Positive finite cutoff in WITHIN_CUTOFF mode; ignored in NEAREST_FOUR mode.
- **q**: Caller buffer, one value per selected centre.
- **q_len**: Length of q.

Returns

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

### seams_traj_unwrap_centres {#seams-traj-unwrap-centres}

```c
seams_status seams_traj_unwrap_centres( const double *xyz, size_t nframes, size_t natoms, const double *boxes, size_t nboxes, const int *mols, size_t mols_len, const double *masses, size_t masses_len, int *mol_ids_out, size_t mol_ids_cap, double *centres_out, size_t centres_len, size_t *n_mols)
```

Assemble per-molecule centres and unwrap them across frames.

Parameters

- **xyz**: Wrapped atom coordinates, nframes * 3 * natoms doubles.
- **nframes**: Number of frames, at least 1.
- **natoms**: Number of atoms per frame, at least 1.
- **boxes**: Cell array of 9 * nboxes doubles.
- **nboxes**: Cell count, either 1 or nframes.
- **mols**: Molecule id per atom, length mols_len.
- **mols_len**: Must equal natoms.
- **masses**: Optional positive mass per atom, or NULL for equal weights.
- **masses_len**: Must equal natoms when masses is non-NULL, or 0 when NULL.
- **mol_ids_out**: Receives unique molecule ids, or NULL.
- **mol_ids_cap**: Capacity of mol_ids_out in ints when non-NULL.
- **centres_out**: Receives nframes * 3 * n_mols unwrapped coordinates, or NULL.
- **centres_len**: Capacity of centres_out in doubles when non-NULL.
- **n_mols**: Receives the molecule count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_traj_unwrap_positions {#seams-traj-unwrap-positions}

```c
seams_status seams_traj_unwrap_positions(const double *xyz, size_t nframes, size_t natoms, const double *boxes, size_t nboxes, double *unwrapped, size_t unwrapped_len)
```

Unwrap periodic coordinates across consecutive frames.

Parameters

- **xyz**: Wrapped coordinates, nframes * 3 * natoms doubles in Angstrom.
- **nframes**: Number of frames, at least 1.
- **natoms**: Number of atoms per frame, at least 1.
- **boxes**: Cell array of 9 * nboxes doubles.
- **nboxes**: Cell count, either 1 or nframes.
- **unwrapped**: Output buffer of unwrapped_len doubles.
- **unwrapped_len**: Capacity in doubles, at least nframes * 3 * natoms.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

### seams_version {#seams-version}

```c
int seams_version(void)
```

Packed release number.

Returns

SEAMS_VERSION_NUMBER, (major << 16) | (minor << 8) | patch.

### seams_version_major {#seams-version-major}

```c
int seams_version_major(void)
```

Major version of the library at run time.

Returns

SEAMS_VERSION_MAJOR of the build.

### seams_version_minor {#seams-version-minor}

```c
int seams_version_minor(void)
```

Minor version of the library at run time.

Returns

SEAMS_VERSION_MINOR of the build.

### seams_version_patch {#seams-version-patch}

```c
int seams_version_patch(void)
```

Patch version of the library at run time.

Returns

SEAMS_VERSION_PATCH of the build.

### seams_version_string {#seams-version-string}

```c
const char *seams_version_string(void)
```

Version text of the built library.

Returns

Pointer to a static string. Never NULL.

### seams_water_vectors {#seams-water-vectors}

```c
seams_status seams_water_vectors(const seams_system *sys, int oxygen_type, int hydrogen_type, int *mol_ids, double *oh_mid, double *bisector, double *hh, double *hh_unit, size_t cap, size_t *n_waters)
```

Water orientation vectors: O to H-midpoint and H-to-H under PBC.

Parameters

- **sys**: Frame handle.
- **oxygen_type**: Species integer of water oxygen.
- **hydrogen_type**: Species integer of water hydrogen, distinct from oxygen_type.
- **mol_ids**: Receives molecule ids, length cap, or NULL.
- **oh_mid**: Receives O-to-H-midpoint vectors in Angstrom, 3 * cap doubles, or NULL.
- **bisector**: Receives unit bisector vectors, 3 * cap doubles, or NULL.
- **hh**: Receives H1-to-H2 vectors in Angstrom, 3 * cap doubles, or NULL.
- **hh_unit**: Receives unit H1-to-H2 vectors, 3 * cap doubles, or NULL.
- **cap**: Molecule capacity of every non-NULL output array.
- **n_waters**: Receives the water count.

Returns

::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.

## C API {#c-api}

Each native call and the `seams.h` or `seams_lang.h` function it uses. The C column is the stable ABI.

| Native | C |
| --- | --- |
| seams_abi | seams_abi |
| seams_brine_analyse | seams_brine_analyse |
| seams_brine_destroy | seams_brine_destroy |
| seams_brine_fronts | seams_brine_fronts |
| seams_brine_pockets | seams_brine_pockets |
| seams_brine_profile | seams_brine_profile |
| seams_brine_rejection | seams_brine_rejection |
| seams_cages | seams_cages |
| seams_cages_counts | seams_cages_counts |
| seams_cages_destroy | seams_cages_destroy |
| seams_cages_find | seams_cages_find |
| seams_cages_labels | seams_cages_labels |
| seams_cages_signature | seams_cages_signature |
| seams_chill | seams_chill |
| seams_chill_plus | seams_chill_plus |
| seams_chill_plus_n | seams_chill_plus_n |
| seams_contact_pairs | seams_contact_pairs |
| seams_coordination | seams_coordination |
| seams_density_axis | seams_density_axis |
| seams_density_fractional | seams_density_fractional |
| seams_directed_hbonds | seams_directed_hbonds |
| seams_domains | seams_domains |
| seams_fingerprint | seams_fingerprint |
| seams_fingerprint_atom | seams_fingerprint_atom |
| seams_fingerprint_destroy | seams_fingerprint_destroy |
| seams_fingerprint_key | seams_fingerprint_key |
| seams_fingerprint_rings | seams_fingerprint_rings |
| seams_get_num_threads | seams_get_num_threads |
| seams_handle_error | seams_handle_error |
| seams_hbond_edges | seams_hbond_edges |
| seams_hbond_tcf_create | seams_hbond_tcf_create |
| seams_hbond_tcf_curve | seams_hbond_tcf_curve |
| seams_hbond_tcf_destroy | seams_hbond_tcf_destroy |
| seams_hbond_tcf_fit | seams_hbond_tcf_fit |
| seams_hbond_tcf_push | seams_hbond_tcf_push |
| seams_hbonds | seams_hbonds |
| seams_host_guest | seams_host_guest |
| seams_label_name | seams_label_name |
| seams_last_error | seams_last_error |
| seams_last_error_message | seams_last_error_message |
| seams_lsi | seams_lsi |
| seams_molecule_centres | seams_molecule_centres |
| seams_molecule_dipoles | seams_molecule_dipoles |
| seams_mols | seams_mols |
| seams_monolayer | seams_monolayer |
| seams_msd_create | seams_msd_create |
| seams_msd_curve | seams_msd_curve |
| seams_msd_destroy | seams_msd_destroy |
| seams_msd_fit | seams_msd_fit |
| seams_msd_push | seams_msd_push |
| seams_neighbours_build | seams_neighbours_build |
| seams_neighbours_copy | seams_neighbours_copy |
| seams_neighbours_count | seams_neighbours_count |
| seams_neighbours_csr | seams_neighbours_csr |
| seams_neighbours_destroy | seams_neighbours_destroy |
| seams_neighbours_knn | seams_neighbours_knn |
| seams_neighbours_nrows | seams_neighbours_nrows |
| seams_orientation_order | seams_orientation_order |
| seams_path_nframes | seams_path_nframes |
| seams_pore_frame | seams_pore_frame |
| seams_pore_profiles | seams_pore_profiles |
| seams_prisms | seams_prisms |
| seams_q4q6 | seams_q4q6 |
| seams_q4q6_n | seams_q4q6_n |
| seams_q6bar | seams_q6bar |
| seams_ql | seams_ql |
| seams_quenched_analysis | seams_quenched_analysis |
| seams_rdf | seams_rdf |
| seams_reorient_create | seams_reorient_create |
| seams_reorient_curve | seams_reorient_curve |
| seams_reorient_destroy | seams_reorient_destroy |
| seams_reorient_fit | seams_reorient_fit |
| seams_reorient_push | seams_reorient_push |
| seams_ring_symmetry_name | seams_ring_symmetry_name |
| seams_rings_chordless | seams_rings_chordless |
| seams_rings_count | seams_rings_count |
| seams_rings_destroy | seams_rings_destroy |
| seams_rings_find | seams_rings_find |
| seams_rings_fraction | seams_rings_fraction |
| seams_rings_from_csr | seams_rings_from_csr |
| seams_rings_histogram | seams_rings_histogram |
| seams_rings_histogram_n | seams_rings_histogram_n |
| seams_rings_labels | seams_rings_labels |
| seams_rings_labels_n | seams_rings_labels_n |
| seams_rings_nrings | seams_rings_nrings |
| seams_rings_vertices | seams_rings_vertices |
| seams_rings_winding | seams_rings_winding |
| seams_running_cn | seams_running_cn |
| seams_set_mols | seams_set_mols |
| seams_set_num_threads | seams_set_num_threads |
| seams_signature_cages | seams_signature_cages |
| seams_signature_cages_edge | seams_signature_cages_edge |
| seams_site_table | seams_site_table |
| seams_status_string | seams_status_string |
| seams_steinhardt | seams_steinhardt |
| seams_system_box | seams_system_box |
| seams_system_count | seams_system_count |
| seams_system_create | seams_system_create |
| seams_system_destroy | seams_system_destroy |
| seams_system_ids | seams_system_ids |
| seams_system_natoms | seams_system_natoms |
| seams_system_positions | seams_system_positions |
| seams_system_read | seams_system_read |
| seams_system_set_atoms | seams_system_set_atoms |
| seams_system_set_frame | seams_system_set_frame |
| seams_system_set_hydrogen_frame | seams_system_set_hydrogen_frame |
| seams_system_set_ortho | seams_system_set_ortho |
| seams_system_set_positions | seams_system_set_positions |
| seams_system_stage_positions | seams_system_stage_positions |
| seams_system_symmetry | seams_system_symmetry |
| seams_system_types | seams_system_types |
| seams_tetrahedral_order | seams_tetrahedral_order |
| seams_traj_unwrap_centres | seams_traj_unwrap_centres |
| seams_traj_unwrap_positions | seams_traj_unwrap_positions |
| seams_version | seams_version |
| seams_version_major | seams_version_major |
| seams_version_minor | seams_version_minor |
| seams_version_patch | seams_version_patch |
| seams_version_string | seams_version_string |
| seams_water_vectors | seams_water_vectors |
