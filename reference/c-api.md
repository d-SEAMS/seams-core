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

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

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

::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

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

Version text, for example “2.14.0”.

Returns

Pointer to a static string. Never NULL.

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
| seams_domains | seams_domains |
| seams_fingerprint | seams_fingerprint |
| seams_fingerprint_atom | seams_fingerprint_atom |
| seams_fingerprint_destroy | seams_fingerprint_destroy |
| seams_fingerprint_key | seams_fingerprint_key |
| seams_fingerprint_rings | seams_fingerprint_rings |
| seams_get_num_threads | seams_get_num_threads |
| seams_handle_error | seams_handle_error |
| seams_hbonds | seams_hbonds |
| seams_label_name | seams_label_name |
| seams_last_error | seams_last_error |
| seams_last_error_message | seams_last_error_message |
| seams_monolayer | seams_monolayer |
| seams_neighbours_build | seams_neighbours_build |
| seams_neighbours_copy | seams_neighbours_copy |
| seams_neighbours_count | seams_neighbours_count |
| seams_neighbours_csr | seams_neighbours_csr |
| seams_neighbours_destroy | seams_neighbours_destroy |
| seams_neighbours_knn | seams_neighbours_knn |
| seams_neighbours_nrows | seams_neighbours_nrows |
| seams_path_nframes | seams_path_nframes |
| seams_prisms | seams_prisms |
| seams_q4q6 | seams_q4q6 |
| seams_q4q6_n | seams_q4q6_n |
| seams_q6bar | seams_q6bar |
| seams_ql | seams_ql |
| seams_rdf | seams_rdf |
| seams_ring_symmetry_name | seams_ring_symmetry_name |
| seams_rings_chordless | seams_rings_chordless |
| seams_rings_count | seams_rings_count |
| seams_rings_destroy | seams_rings_destroy |
| seams_rings_find | seams_rings_find |
| seams_rings_fraction | seams_rings_fraction |
| seams_rings_histogram | seams_rings_histogram |
| seams_rings_histogram_n | seams_rings_histogram_n |
| seams_rings_labels | seams_rings_labels |
| seams_rings_labels_n | seams_rings_labels_n |
| seams_rings_nrings | seams_rings_nrings |
| seams_rings_vertices | seams_rings_vertices |
| seams_running_cn | seams_running_cn |
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
| seams_system_set_ortho | seams_system_set_ortho |
| seams_system_set_positions | seams_system_set_positions |
| seams_system_stage_positions | seams_system_stage_positions |
| seams_system_symmetry | seams_system_symmetry |
| seams_system_types | seams_system_types |
| seams_version | seams_version |
| seams_version_major | seams_version_major |
| seams_version_minor | seams_version_minor |
| seams_version_patch | seams_version_patch |
| seams_version_string | seams_version_string |
