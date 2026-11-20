Further analyses
================


Brine fronts, topology, and the remaining exports.

The stable C interface provides 89 exported functions through ``include/seams.h``.

.. _capi-seams_brine_analyse:

seams_brine_analyse
-------------------


.. code-block:: c

   seams_status seams_brine_analyse(seams_system *sys, const int *ice, size_t n_ice, int water_type, int ion_type, int n_bins, int axis, seams_brine **out)

Planar ice/brine front from an ice assignment.

``sys``
    Frame handle.

``ice``
    Per-atom ice flags, length n_ice.

``n_ice``
    Must equal the frame count.

``water_type``
    Species kept as water, or 0.

``ion_type``
    Species kept as ions, or 0 for none.

``n_bins``
    Bin count along axis, or 0 for 20.

``axis``
    0, 1, or 2.

``out``
    Receives a new handle. On failure \*out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

ice[i] is non-zero when atom i is ice. water_type 0 keeps every atom whose type is not ion_type. ion_type 0 means the frame has no ions. n_bins 0 selects 20 bins. axis is 0, 1, or 2. The ice fraction cut is 0.5, and the bulk ice and liquid cuts are 0.85 and 0.15.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees \*out with seams_brine_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_brine \*front = NULL; if (seams_brine_analyse(sys, ice, n, 1, 2, 0, 2, &front) == SEAMS_OK) { seams_brine_destroy(front); } @endcode

.. _capi-seams_brine_destroy:

seams_brine_destroy
-------------------


.. code-block:: c

   void seams_brine_destroy(seams_brine *front)

Free a brine-front handle.

``front``
    Handle from seams_brine_analyse(), or NULL. NULL is ignored.

@par Ownership Releases the handle and its profile. @par Thread safety The caller must not use front after this returns. @par Example @code seams_brine_destroy(front); @endcode

.. _capi-seams_brine_fronts:

seams_brine_fronts
------------------


.. code-block:: c

   seams_status seams_brine_fronts(const seams_brine *front, double *z_enter, double *z_exit, double *thickness, int *found)

The two planar crossings of the ice slab.

``front``
    Brine handle.

``z_enter``
    Where the ice fraction crosses the cut on entering the ice along the positive axis, or NULL.

``z_exit``
    The crossing on leaving the ice, or NULL.

``thickness``
    Ice run along the axis, or NULL.

``found``
    Receives 1 when both crossings exist, otherwise 0, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or every output pointer is NULL.

Lengths are Angstrom, in the frame's cell.

@par Ownership The caller owns front; the call borrows it. z_enter, z_exit, thickness and found are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on front, and not alongside a non-const call on front. @par Example @code double thick = 0; int found = 0; seams_brine_fronts(front, NULL, NULL, &thick, &found); @endcode

.. _capi-seams_brine_rejection:

seams_brine_rejection
---------------------


.. code-block:: c

   seams_status seams_brine_rejection(const seams_brine *front, double *c_ice, double *c_liquid, double *k, int *defined)

Molal rejection coefficient k = c_ice / c_liquid.

``front``
    Brine handle.

``c_ice``
    Bulk ice molality, or NULL.

``c_liquid``
    Bulk liquid molality, or NULL.

``k``
    The ratio, or NULL. Unchanged when the ratio is undefined.

``defined``
    Receives 1 when the liquid bulk contains salt or the ice bulk does not, otherwise 0, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or every output pointer is NULL.

Molality is formula units per kilogram of water.

@par Ownership The caller owns front; the call borrows it. c_ice, c_liquid, k and defined are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on front, and not alongside a non-const call on front. @par Example @code double k = 0; int defined = 0; seams_brine_rejection(front, NULL, NULL, &k, &defined); @endcode

.. _capi-seams_brine_profile:

seams_brine_profile
-------------------


.. code-block:: c

   seams_status seams_brine_profile(const seams_brine *front, double *z_mid, double *molality, double *ice_fraction, size_t cap, size_t *n_out)

Concentration profile along the front axis.

``front``
    Brine handle.

``z_mid``
    Bin centers, or NULL to skip the column.

``molality``
    Bin molalities, or NULL to skip the column.

``ice_fraction``
    Bin ice fractions, or NULL to skip the column.

``cap``
    Capacity of each non-NULL column.

``n_out``
    Receives the bin count.

Returns: ::SEAMS_OK when every non-NULL column was filled, or when every column pointer was NULL and n_out was written. ::SEAMS_ERR_RANGE when a non-NULL column is short. ::SEAMS_ERR_NULL when front or n_out is NULL.

z_mid is Angstrom. molality is mol/kg. ice_fraction is dimensionless.

@par Ownership The caller owns the columns. @par Thread safety Read-only: may run alongside other const calls on front, and not alongside a non-const call on front. @par Example @code size_t bins = 0; seams_brine_profile(front, NULL, NULL, NULL, 0, &bins); @endcode

.. _capi-seams_brine_pockets:

seams_brine_pockets
-------------------


.. code-block:: c

   seams_status seams_brine_pockets(const seams_brine *front, int *n_pockets, int *reservoir_waters)

Brine pockets trapped in the ice, excluding the liquid reservoir.

``front``
    Brine handle.

``n_pockets``
    Receives the pocket count, or NULL.

``reservoir_waters``
    Receives the water count of the liquid reservoir, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL when front is NULL or both output pointers are NULL.

@par Ownership The caller owns front; the call borrows it. n_pockets and reservoir_waters are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on front, and not alongside a non-const call on front. @par Example @code int pockets = 0; seams_brine_pockets(front, &pockets, NULL); @endcode

.. _capi-seams_rdf:

seams_rdf
---------


.. code-block:: c

   seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *g, int *count, size_t cap, size_t *n_out, double *volume)

Partial radial distribution.

``sys``
    Frame handle.

``rmax``
    Positive finite outer radius, Angstrom.

``bins``
    Bin count. Zero selects about rmax / 0.1.

``type_i``
    First species. Zero and negative keep every atom of a single-species frame, or every atom as one species otherwise.

``type_j``
    Second species. Zero and negative reuse the type_i species.

``r``
    Bin centers, or NULL to query.

``g``
    g(r) per bin, or NULL to query.

``count``
    Pair count per bin, or NULL to query.

``cap``
    Capacity of each non-NULL column.

``n_out``
    Receives the bin count.

``volume``
    Receives the cell volume, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

r is Angstrom. g is dimensionless. volume is cubic Angstrom.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns the columns. @par Thread safety Exclusive use of sys. @par Example @code size_t bins = 0; seams_rdf(sys, 8.0, 80, 0, 0, NULL, NULL, NULL, 0, &bins, NULL); @endcode

.. _capi-seams_running_cn:

seams_running_cn
----------------


.. code-block:: c

   seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *cn, size_t cap, size_t *n_out)

Running coordination number along the partial RDF.

``sys``
    Frame handle.

``rmax``
    Positive finite outer radius, Angstrom.

``bins``
    Bin count. Zero selects about rmax / 0.1.

``type_i``
    Central species, as in seams_rdf().

``type_j``
    Counted species, as in seams_rdf().

``r``
    Bin centers, Angstrom, or NULL to query.

``cn``
    Running coordination number, or NULL to query.

``cap``
    Capacity of each non-NULL column.

``n_out``
    Receives the bin count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns the columns. @par Thread safety Exclusive use of sys. @par Example @code double r[80], cn[80]; size_t bins = 0; seams_running_cn(sys, 8.0, 80, 0, 0, r, cn, 80, &bins); @endcode

.. _capi-seams_coordination:

seams_coordination
------------------


.. code-block:: c

   seams_status seams_coordination(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *cn)

Coordination number integrated to rmax.

``sys``
    Frame handle.

``rmax``
    Positive finite outer radius, Angstrom.

``bins``
    Bin count. Zero selects about rmax / 0.1.

``type_i``
    Central species, as in seams_rdf().

``type_j``
    Counted species, as in seams_rdf().

``cn``
    Receives the coordination number.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys; the call borrows it. cn is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Exclusive use of sys. @par Example @code double cn = 0; seams_coordination(sys, 3.5, 0, 1, 1, &cn); @endcode

.. _capi-seams_hbonds:

seams_hbonds
------------


.. code-block:: c

   seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int *n_bonds)

Hydrogen-bond count on one frame.

``sys``
    Frame handle.

``cutoff``
    Positive finite donor-acceptor cutoff, Angstrom.

``oxygen_type``
    Heavy-atom species. Zero keeps every atom that is not hydrogen_type.

``hydrogen_type``
    Positive hydrogen species.

``dist_cutoff``
    Acceptor-H distance, Angstrom.

``angle_cutoff``
    O-O-H angle, degrees.

``n_bonds``
    Receives the number of hydrogen bonds.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys; the call borrows it. n_bonds is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Exclusive use of sys. @par Example @code int n_bonds = 0; seams_hbonds(sys, 3.5, 1, 2, 2.4, 30.0, &n_bonds); @endcode

.. _capi-seams_hbond_edges:

seams_hbond_edges
-----------------


.. code-block:: c

   seams_status seams_hbond_edges(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int fields, seams_hbond_edge *edges, size_t cap, size_t *n_out)

Hydrogen-bond edges for one frame.

``sys``
    Frame handle.

``cutoff``
    Donor-acceptor cutoff, Angstrom.

``oxygen_type``
    Heavy-atom species filter as in seams_hbonds().

``hydrogen_type``
    Hydrogen species (positive).

``dist_cutoff``
    Acceptor-H distance, Angstrom.

``angle_cutoff``
    O-O-H angle, degrees.

``fields``
    Bitmask: SEAMS_HBOND_EDGE_WITH_H and/or SEAMS_HBOND_EDGE_WITH_GEOMETRY.

``edges``
    Caller buffer with caller-set size on each element, or NULL to query.

``cap``
    Capacity of @p edges. Ignored when @p edges is NULL.

``n_out``
    Receives the edge count.

Returns: ::SEAMS_OK, or the same error classes as seams_hbonds(), plus ::SEAMS_ERR_ARG when fields is unknown or size is invalid, and ::SEAMS_ERR_RANGE when cap is smaller than \*n_out.

The pair list is stable and deterministic: donor, then acceptor, then hydrogen index. donor and acceptor are frame indices over the stored atoms. When edges is non-NULL, every element must start with size set by the caller to sizeof(seams_hbond_edge); the library validates size and writes only up to that prefix per element.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys, edges, and n_out. The call borrows sys and writes caller buffers only. @par Thread safety Exclusive use of sys.

.. _capi-seams_directed_hbonds:

seams_directed_hbonds
---------------------


.. code-block:: c

   seams_status seams_directed_hbonds( seams_system *sys, double pair_cutoff, double dist_cutoff, double angle_cutoff, const int *donor_types, size_t n_donor_types, const int *acceptor_types, size_t n_acceptor_types, int hydrogen_type, int exclude_same_mol, seams_directed_hbond *edges, size_t cap, size_t *n_out)

Directed hydrogen-bond edges for selected donor and acceptor types.

``sys``
    Frame handle.

``pair_cutoff``
    Positive finite donor-acceptor cutoff, Angstrom.

``dist_cutoff``
    Positive finite acceptor-hydrogen cutoff, Angstrom.

``angle_cutoff``
    Positive finite acceptor-first angle, degrees.

``donor_types``
    Species that may donate. NULL only when the count is 0.

``n_donor_types``
    Length of donor_types. Zero writes no edges.

``acceptor_types``
    Species that may accept. NULL only when the count is 0.

``n_acceptor_types``
    Length of acceptor_types. Zero writes no edges.

``hydrogen_type``
    Positive hydrogen species, distinct from both type sets.

``exclude_same_mol``
    1 drops a shared molID, 0 keeps it.

``edges``
    Caller buffer, or NULL to query the count.

``cap``
    Capacity of edges. Ignored when edges is NULL.

``n_out``
    Receives the edge count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

Each edge is one accepting hydrogen. Rows are sorted by donor, then hydrogen, then acceptor, as frame indices. The angle is acceptor-first. The donor-acceptor vector is the minimum image of acceptor minus donor. The acceptor-hydrogen vector is the minimum image of acceptor minus the hydrogen. A linear donor-hydrogen-acceptor geometry has angle 0 degrees. The acceptor-hydrogen distance is kept only when it is strictly below dist_cutoff. The angle is kept when it is less than or equal to angle_cutoff. The pair list keeps a donor and acceptor only when their minimum image is strictly below pair_cutoff, on an orthorhombic or triclinic cell. Each hydrogen is assigned to the nearest donor-type heavy atom that shares its molID, and the lower index wins a tie. exclude_same_mol of 1 drops a pair that shares a molID. An empty donor set or acceptor set writes zero edges. Atoms whose type equals hydrogen_type are hydrogens and are never donors or acceptors. Water-water and water-framework counts are two calls with different type sets. When edges is non-NULL, every element must start with size set to sizeof(seams_directed_hbond). @par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys, the type arrays, edges, and n_out. The call borrows them and writes caller buffers only. @par Thread safety Exclusive use of sys. @par Example @code const int donors[] = {8}; const int acceptors[] = {10}; size_t n_out = 0; seams_directed_hbonds(sys, 3.5, 2.42, 30.0, donors, 1, acceptors, 1, 2, 1, NULL, 0, &n_out); @endcode

.. _capi-seams_quenched_analysis:

seams_quenched_analysis
-----------------------


.. code-block:: c

   seams_status seams_quenched_analysis( seams_system *sys, seams_quench_callback callback, void *user_data, double cutoff, int type, int graph, int max_depth, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int *raw_labels, size_t raw_labels_len, int *quenched_labels, size_t quenched_labels_len, int counts_raw[9], int counts_quenched[9], int *rings_raw, int *rings_quenched, size_t ring_counts_len, int *hbond_raw, int *hbond_quenched, int *hbond_edge_delta, int *ring_count_delta, int *class_change_count)

Run CHILL+, ring census, and H-bonds on raw and quenched frames.

``sys``
    Frame handle.

``callback``
    Quench callback that writes output coordinates.

``user_data``
    Opaque callback payload.

``cutoff``
    Positive finite neighbour cutoff, Angstrom.

``type``
    Species filter for neighbour, CHILL+, and rings.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

``max_depth``
    Primitive-ring depth, from 3 through 12.

``oxygen_type``
    Heavy-atom species filter for H-bonds.

``hydrogen_type``
    Hydrogen species filter for H-bonds.

``dist_cutoff``
    Acceptor-H distance, Angstrom.

``angle_cutoff``
    O-O-H angle, degrees.

``raw_labels``
    Raw-frame CHILL+ labels, length at least atom count.

``raw_labels_len``
    Length of @p raw_labels.

``quenched_labels``
    Quenched-frame CHILL+ labels, length at least atom count.

``quenched_labels_len``
    Length of @p quenched_labels.

``counts_raw``
    Raw-frame CHILL+ population counts.

``counts_quenched``
    Quenched-frame CHILL+ population counts.

``rings_raw``
    Raw-frame primitive-ring histogram.

``rings_quenched``
    Quenched-frame primitive-ring histogram.

``ring_counts_len``
    Length of @p rings_raw and @p rings_quenched.

``hbond_raw``
    Receives raw-frame H-bond edge count.

``hbond_quenched``
    Receives quenched-frame H-bond edge count.

``hbond_edge_delta``
    Receives symmetric-difference edge count between raw and quenched H-bond edge sets.

``ring_count_delta``
    Receives sum over ring sizes of absolute raw-minus- quenched count differences.

``class_change_count``
    Receives the atom count whose CHILL+ class changed.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_BOX, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

The callback receives one frame and returns quenched coordinates. This call then runs CHILL+, primitive rings, and hydrogen-bond analysis on both raw and quenched frames and emits their per-frame differences. @par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns every output buffer, callback state, and sys. The call borrows callback pointers and writes caller buffers only. @par Thread safety Exclusive use of sys.

.. _capi-seams_density_axis:

seams_density_axis
------------------


.. code-block:: c

   seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis, double *z, double *rho, size_t cap, size_t *n_out)

Number density along one cell axis.

``sys``
    Frame handle.

``type``
    Species filter. Zero keeps every atom.

``bins``
    Bin count. Zero selects about one bin per 0.1 Angstrom.

``axis``
    0, 1, or 2.

``z``
    Bin centers, or NULL to query.

``rho``
    Density per bin, or NULL to query.

``cap``
    Capacity of each non-NULL column.

``n_out``
    Receives the bin count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

z is Angstrom. rho is number density, atoms per cubic Angstrom.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns the columns. @par Thread safety Exclusive use of sys. @par Example @code double z[50], rho[50]; size_t bins = 0; seams_density_axis(sys, 0, 50, 2, z, rho, 50, &bins); @endcode

.. _capi-seams_fingerprint:

seams_fingerprint
-----------------


.. code-block:: c

   seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k, int graph, int hops, int color_types, seams_topo **out)

Topological fingerprint of the bonded graph.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter.

``k``
    Neighbour count for the mutual graph. Zero selects 4.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

``hops``
    Neighbourhood depth. Zero selects 2.

``color_types``
    Non-zero colors vertices by LAMMPS type.

``out``
    Receives a ::seams_topo. On failure \*out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees \*out with seams_fingerprint_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_topo \*fp = NULL; seams_fingerprint(sys, 3.5, 0, 4, SEAMS_GRAPH_CUTOFF, 2, 0, &fp); seams_fingerprint_destroy(fp); @endcode

.. _capi-seams_fingerprint_destroy:

seams_fingerprint_destroy
-------------------------


.. code-block:: c

   void seams_fingerprint_destroy(seams_topo *fp)

Free a fingerprint.

``fp``
    Handle from seams_fingerprint(), or NULL. NULL is ignored.

@par Ownership Releases the handle and its key strings. @par Thread safety The caller must not use fp after this returns. @par Example @code seams_fingerprint_destroy(fp); @endcode

.. _capi-seams_fingerprint_key:

seams_fingerprint_key
---------------------


.. code-block:: c

   seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap, size_t *n_out)

Frame key text.

``fp``
    Fingerprint handle.

``buf``
    Caller buffer of cap bytes, or NULL to query.

``cap``
    Capacity of buf, including the terminating NUL.

``n_out``
    Receives the key length without the terminating NUL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when buf is short.

@par Ownership The caller owns buf. The copy is NUL-terminated. @par Thread safety Read-only: may run alongside other const calls on fp, and not alongside a non-const call on fp. @par Example @code size_t len = 0; seams_fingerprint_key(fp, NULL, 0, &len); char \*key = malloc(len + 1); seams_fingerprint_key(fp, key, len + 1, &len); free(key); @endcode

.. _capi-seams_fingerprint_atom:

seams_fingerprint_atom
----------------------


.. code-block:: c

   seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf, size_t cap, size_t *n_out)

Local key of one atom. Same buffer rule as seams_fingerprint_key().

``fp``
    Fingerprint handle.

``index``
    Atom index in the frame.

``buf``
    Caller buffer of cap bytes, or NULL to query.

``cap``
    Capacity of buf, including the terminating NUL.

``n_out``
    Receives the key length without the terminating NUL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE for a bad index or a short buffer.

@par Ownership The caller owns buf. @par Thread safety Read-only: may run alongside other const calls on fp, and not alongside a non-const call on fp. @par Example @code char key[256]; size_t len = 0; seams_fingerprint_atom(fp, 0, key, sizeof key, &len); @endcode

.. _capi-seams_fingerprint_rings:

seams_fingerprint_rings
-----------------------


.. code-block:: c

   seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts, size_t counts_len)

Primitive-ring census stored on the fingerprint.

``fp``
    Fingerprint handle.

``counts``
    Caller buffer. counts[k] receives the number of rings of size k. Sizes past the buffer are dropped.

``counts_len``
    Length of counts. Must be positive.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when counts_len is 0.

@par Ownership The caller owns counts. @par Thread safety Read-only: may run alongside other const calls on fp, and not alongside a non-const call on fp. @par Example @code int census[13]; seams_fingerprint_rings(fp, census, 13); @endcode

.. _capi-seams_contact_pairs:

seams_contact_pairs
-------------------


.. code-block:: c

   seams_status seams_contact_pairs(seams_system *sys, const char *site_spec, size_t *n_pairs, int *n_cation, int *n_anion)

Mutual nearest unlike contacts of a site table.

``sys``
    Frame handle.

``site_spec``
    Text such as "1=cationHead,2=anion". Must be non-NULL.

``n_pairs``
    Receives the contact count, or NULL.

``n_cation``
    Receives the cation molecule count, or NULL.

``n_anion``
    Receives the anion molecule count, or NULL. At least one of the three outputs must be non-NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

Cation and anion counts are molecules in the ion cloud.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys; the call borrows it. site_spec is read during the call only. n_pairs, n_cation and n_anion are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Exclusive use of sys. @par Example @code size_t pairs = 0; seams_contact_pairs(sys, "1=cationHead,2=anion", &pairs, NULL, NULL); @endcode

.. _capi-seams_domains:

seams_domains
-------------


.. code-block:: c

   seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec, int subset, int *n_subset, int *largest, double *percolation)

Largest domain of a polar or apolar site subset.

``sys``
    Frame handle.

``cutoff``
    Positive finite cluster cutoff, Angstrom.

``site_spec``
    Site table text. Must be non-NULL.

``subset``
    0 is apolar, 1 is polar.

``n_subset``
    Receives the number of sites in the subset, or NULL.

``largest``
    Receives the size of the largest domain, or NULL.

``percolation``
    Receives the largest domain size over the subset size, or NULL. At least one output must be non-NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys; the call borrows it. site_spec is read during the call only. n_subset, largest and percolation are caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Exclusive use of sys. @par Example @code int largest = 0; seams_domains(sys, 5.0, "1=cationHead,2=anion,3=tail", 0, NULL, &largest, NULL); @endcode

.. _capi-seams_site_table:

seams_site_table
----------------


.. code-block:: c

   seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds, size_t len)

Per-atom site kind from a site-table spec.

``sys``
    Frame handle.

``spec``
    Site table text such as "1=cationHead,2=anion". Must be non-NULL.

``kinds``
    One ::seams_site_kind per atom.

``len``
    Length of kinds. Must be at least the frame count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns kinds. @par Thread safety Exclusive use of sys. @par Example @code seams_site_table(sys, "1=cationHead,2=anion", kinds, n); @endcode

