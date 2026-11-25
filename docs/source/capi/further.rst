Further analyses
================


Brine fronts, topology, and the remaining exports.

The stable C interface provides 125 exported functions through ``include/seams.h``.

.. _capi-seams_tetrahedral_order:

seams_tetrahedral_order
-----------------------


.. code-block:: c

   seams_status seams_tetrahedral_order( seams_system *sys, const int *centres, size_t n_centres, const int *partners, size_t n_partners, int convention, double cutoff, double *q, size_t q_len)

Errington-Debenedetti tetrahedral order for arbitrary selections.

``sys``
    Frame handle.

``centres``
    Centre indices, or NULL when n_centres is zero.

``n_centres``
    Centre count. Zero means every atom.

``partners``
    Partner indices, or NULL when n_partners is zero.

``n_partners``
    Partner count. Zero means every atom.

``convention``
    One ::seams_tetrahedral_neighbours value. Any other value returns ::SEAMS_ERR_ARG.

``cutoff``
    Positive finite cutoff in WITHIN_CUTOFF mode; ignored in NEAREST_FOUR mode.

``q``
    Caller buffer, one value per selected centre.

``q_len``
    Length of q.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

centres and partners contain zero-based frame indices. A NULL pointer with zero length selects every atom. A centre is excluded from its own partner shell. q = 1 - 9/(2 N (N-1)) sum_{j<k} (cos psi_jk + 1/3)^2 over the N shell neighbours: a perfect tetrahedron is 1 and an ideal gas averages 0. NEAREST_FOUR uses the four nearest selected partners, N = 4 and the Errington-Debenedetti factor 3/8, and ignores cutoff. WITHIN_CUTOFF uses every selected partner at or inside cutoff. In both modes a centre with fewer than four neighbours receives quiet NaN.

@par Handle access Non-const only because the call clears the handle error slot. @par Ownership The caller owns every selection and output buffer. @par Thread safety Exclusive use of sys.

.. _capi-seams_lsi:

seams_lsi
---------


.. code-block:: c

   seams_status seams_lsi(seams_system *sys, const int *centres, size_t n_centres, const int *partners, size_t n_partners, double cutoff, double *lsi, size_t lsi_len)

Local Structure Index for arbitrary selections.

``sys``
    Frame handle.

``centres``
    Centre indices, or NULL when n_centres is zero.

``n_centres``
    Centre count. Zero means every atom.

``partners``
    Partner indices, or NULL when n_partners is zero.

``n_partners``
    Partner count. Zero means every atom.

``cutoff``
    Positive finite shell cutoff.

``lsi``
    Caller buffer, one value per selected centre.

``lsi_len``
    Length of lsi.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

Shiratani-Sasai index. For each centre, selected-partner distances are sorted, r_1 < ... < r_N <= cutoff < r_{N+1}. The result is the population variance of the N gaps r_{i+1} - r_i, i = 1..N; the last gap reaches the first partner beyond cutoff. A centre with no partner inside cutoff, or none beyond it, receives quiet NaN.

@par Handle access Non-const only because the call clears the handle error slot. @par Ownership The caller owns every selection and output buffer. @par Thread safety Exclusive use of sys.

.. _capi-seams_switch_rational:

seams_switch_rational
---------------------


.. code-block:: c

   seams_status seams_switch_rational(double r, double r0, double d0, int n, int m, double d_on, double d_max, double *s, double *dsdr)

PLUMED RATIONAL switching function with a C2 taper to zero.

``r``
    Distance, Angstrom. Not negative.

``r0``
    Positive length scale, Angstrom.

``d0``
    Plateau radius, Angstrom, at least zero.

``n``
    Numerator exponent, at least one.

``m``
    Denominator exponent, greater than n. Zero selects 2 \* n.

``d_on``
    Taper start, Angstrom, from d0 and below d_max.

``d_max``
    Radius where s reaches zero, Angstrom.

``s``
    Receives s(r).

``dsdr``
    Receives ds/dr in 1/Angstrom, the right derivative at r = d0, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

For x = (r - d0) / r0, s = (1 - x^n) / (1 - x^m) for r > d0 and 1 for r <= d0, times a taper that is 1 below d_on, 0 at and beyond d_max, and 1 - u^3 (10 - 15u + 6u^2) for u = (r - d_on) / (d_max - d_on) between them. s and its first two derivatives vanish at d_max. On r >= 0, s is C^2 in r except at r = d0 > 0 with n < 3, where it is only C^(n-1). At r = d0, dsdr is the right derivative: -1 / r0 for n = 1 and 0 otherwise. @par Ownership The caller owns s and dsdr. The call keeps no state. @par Thread safety Safe from any thread. @par Example @code double sw = 0, dsw = 0; seams_switch_rational(3.0, 0.5, 2.5, 6, 0, 3.3, 3.5, &sw, &dsw); @endcode

.. _capi-seams_cv_create:

seams_cv_create
---------------


.. code-block:: c

   seams_status seams_cv_create(double r0, double d0, int n, int m, double d_on, double d_max, int max_ring, seams_cv **out)

Create a smooth collective-variable workspace.

``r0``
    Switching length scale, Angstrom. See seams_switch_rational().

``d0``
    Switching plateau radius, Angstrom.

``n``
    Switching numerator exponent.

``m``
    Switching denominator exponent; zero selects 2 \* n.

``d_on``
    Taper start, Angstrom.

``d_max``
    Radius where every bond weight reaches zero, Angstrom.

``max_ring``
    Largest candidate ring, from 3 through 10.

``out``
    Receives a new handle. On failure \*out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM.

@par Ownership The caller frees \*out with seams_cv_destroy(). The handle holds no frame until seams_cv_rebuild(). @par Thread safety The new handle is private to the caller. @par Example @code seams_cv \*cv = NULL; if (seams_cv_create(0.5, 2.5, 6, 0, 3.3, 3.5, 6, &cv) == SEAMS_OK) { seams_cv_destroy(cv); } @endcode

.. _capi-seams_cv_destroy:

seams_cv_destroy
----------------


.. code-block:: c

   void seams_cv_destroy(seams_cv *cv)

Free a collective-variable workspace.

``cv``
    Handle from seams_cv_create(), or NULL. NULL is ignored.

@par Ownership Releases the handle, its topology and its state clauses. @par Thread safety The caller must not use cv after this returns. @par Example @code seams_cv \*cv = NULL; seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv); seams_cv_destroy(cv); @endcode

.. _capi-seams_cv_rebuild:

seams_cv_rebuild
----------------


.. code-block:: c

   seams_status seams_cv_rebuild(seams_cv *cv, const seams_system *sys, int type, const int *atoms, size_t n_atoms, double skin)

List the graph and fix the candidate rings for the next frames.

``cv``
    Workspace handle.

``sys``
    Frame handle that holds a frame.

``type``
    Species to select. Zero and negative select every atom.

``atoms``
    Frame indices to select, or NULL for every atom of type.

``n_atoms``
    Length of atoms. Ignored when atoms is NULL.

``skin``
    Verlet skin, Angstrom, at least zero. d_max + skin must be below half the narrowest width of the cell.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_BOX, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

Bonds are the selected pairs within d_max + skin, from linkcell. The candidate rings, their shortcut paths and the hexagon pairs are fixed until the next rebuild. A value computed from this topology equals the value from a fresh rebuild while no selected atom has moved more than skin / 2 from its position here and the cell is unchanged; seams_cv_needs_rebuild() reports when that no longer holds, and seams_cv_value() then refuses the frame. @par Handle access Writes the topology into cv. Reads sys. @par Ownership The caller owns atoms. The workspace keeps a copy of the selection and of the coordinates and cell of this frame. @par Thread safety Exclusive use of cv. Read-only use of sys. @par Example @code seams_cv \*cv = NULL; seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv); seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1); seams_cv_destroy(cv); @endcode

.. _capi-seams_cv_needs_rebuild:

seams_cv_needs_rebuild
----------------------


.. code-block:: c

   seams_status seams_cv_needs_rebuild(const seams_cv *cv, const seams_system *sys, int *stale)

Whether a frame needs seams_cv_rebuild() before the next value.

``cv``
    Workspace handle.

``sys``
    Frame handle that holds a frame.

``stale``
    Receives 1 when cv was never rebuilt, the atom count or cell changed, or a selected atom moved more than skin / 2 since the rebuild, and 0 otherwise.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE.

@par Ownership The caller owns stale. @par Thread safety Read-only: may run alongside other const calls on cv and sys. @par Example @code seams_cv \*cv = NULL; seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv); seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1); int stale = 0; seams_cv_needs_rebuild(cv, sys, &stale); seams_cv_destroy(cv); @endcode

.. _capi-seams_cv_topology:

seams_cv_topology
-----------------


.. code-block:: c

   seams_status seams_cv_topology(const seams_cv *cv, size_t *n_edges, size_t *n_rings, size_t *n_pairs)

Sizes of the fixed topology.

``cv``
    Workspace handle.

``n_edges``
    Receives the number of listed bonds, or NULL.

``n_rings``
    Receives the number of candidate rings, or NULL.

``n_pairs``
    Receives the number of candidate hexagon pairs, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE before the first rebuild.

@par Ownership The caller owns the three outputs. @par Thread safety Read-only: may run alongside other const calls on cv. @par Example @code seams_cv \*cv = NULL; seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv); seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1); size_t edges = 0, rings = 0, pairs = 0; seams_cv_topology(cv, &edges, &rings, &pairs); seams_cv_destroy(cv); @endcode

.. _capi-seams_cv_value:

seams_cv_value
--------------


.. code-block:: c

   seams_status seams_cv_value(seams_cv *cv, const seams_system *sys, int kind, int arg, double *value, double *grad, size_t grad_len)

Value and gradient of one smooth collective variable.

``cv``
    Workspace handle after seams_cv_rebuild().

``sys``
    Frame with the atom count and cell of the rebuild.

``kind``
    A ::seams_cv_kind code.

``arg``
    Ring size, weight choice, or state index, as above.

``value``
    Receives the value.

``grad``
    Caller buffer of 3 \* natoms doubles for d value / d x in 1/Angstrom times the unit of value, or NULL to skip it.

``grad_len``
    Length of grad. Ignored when grad is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE when cv needs a rebuild for sys, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

::SEAMS_CV_RING_COUNT sums the ring weights of size arg. ::SEAMS_CV_RING_ATOMS sums 1 - prod (1 - w) over the rings of size arg through each selected atom. ::SEAMS_CV_HEX_BASAL_PAIRS sums w_a w_b B over vertex-disjoint hexagon pairs, B the soft form of the hexagonal-cage basal test of either ring against the other; arg is 0. ::SEAMS_CV_GYRATION_RG2 and ::SEAMS_CV_GYRATION_VOLUME use the gyration tensor G of the selected atoms, unwrapped along bonds, with weight 1 when arg is -1 and the ring participation of size arg otherwise; the volume is (4 pi / 3) 5^(3/2) sqrt(det G), the volume of the uniform ellipsoid with that tensor. ::SEAMS_CV_STATE returns state arg. @par Handle access Non-const because the evaluation keeps bond-weight scratch on cv. The topology is not changed. @par Ownership The caller owns value and grad. @par Thread safety Exclusive use of cv. Read-only use of sys. @par Example @code seams_cv \*cv = NULL; seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv); seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1); double n6 = 0; seams_cv_value(cv, sys, SEAMS_CV_RING_COUNT, 6, &n6, NULL, 0); seams_cv_destroy(cv); @endcode

.. _capi-seams_cv_state_clause:

seams_cv_state_clause
---------------------


.. code-block:: c

   seams_status seams_cv_state_clause(seams_cv *cv, int state, int kind, int arg, double threshold, double width, int above)

Add a soft threshold factor to a state.

``cv``
    Workspace handle.

``state``
    Index of the state to extend, or the current number of states to open a new one.

``kind``
    A ::seams_cv_kind code.

``arg``
    Argument of kind, as for seams_cv_value(); for ::SEAMS_CV_STATE an index below state.

``threshold``
    Threshold, in the unit of the variable.

``width``
    Positive logistic width, in the unit of the variable. Zero only for ::SEAMS_CV_STATE.

``above``
    Nonzero for a factor that rises through the threshold.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM.

A state is the product of its factors. A factor on a collective variable v is sigma((v - threshold) / width) when above is nonzero and 1 - sigma((v - threshold) / width) otherwise, sigma(z) = 1 / (1 + e^-z). A factor of kind ::SEAMS_CV_STATE reads an earlier state chi; with width zero it is chi itself, or 1 - chi when above is zero, so states compose by product (and), complement (not), and 1 - (1 - a)(1 - b) (or). @par Handle access Appends to the state list of cv. @par Ownership The workspace keeps the clause until it is destroyed. @par Thread safety Exclusive use of cv. @par Example @code seams_cv \*cv = NULL; seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv); seams_cv_state_clause(cv, 0, SEAMS_CV_RING_COUNT, 6, 40.0, 4.0, 1); seams_cv_destroy(cv); @endcode

.. _capi-seams_cv_states:

seams_cv_states
---------------


.. code-block:: c

   seams_status seams_cv_states(seams_cv *cv, const seams_system *sys, double *chi, size_t chi_len, size_t *n_states, double *grad, size_t grad_len)

Soft indicator of every state, with gradients on request.

``cv``
    Workspace handle after seams_cv_rebuild().

``sys``
    Frame with the atom count and cell of the rebuild.

``chi``
    Caller buffer of one double per state, each in [0, 1], or NULL to query n_states.

``chi_len``
    Length of chi. Ignored when chi is NULL.

``n_states``
    Receives the number of states.

``grad``
    Caller buffer of n_states rows of 3 \* natoms doubles, d chi / d x in 1/Angstrom, or NULL to skip it.

``grad_len``
    Length of grad. Ignored when grad is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE when cv needs a rebuild for sys, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const because the evaluation keeps bond-weight scratch on cv. @par Ownership The caller owns chi, n_states and grad. @par Thread safety Exclusive use of cv. Read-only use of sys. @par Example @code seams_cv \*cv = NULL; seams_cv_create(0.3, 2.8, 6, 0, 2.9, 3.0, 6, &cv); seams_cv_rebuild(cv, sys, 0, NULL, 0, 0.1); seams_cv_state_clause(cv, 0, SEAMS_CV_RING_COUNT, 6, 40.0, 4.0, 1); seams_cv_state_clause(cv, 1, SEAMS_CV_STATE, 0, 0.0, 0.0, 0); double chi[2]; size_t n_states = 0; seams_cv_states(cv, sys, chi, 2, &n_states, NULL, 0); seams_cv_destroy(cv); @endcode

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

Each edge is one accepting hydrogen. Rows are sorted by donor, then hydrogen, then acceptor, as frame indices. The angle is acceptor-first. The donor-acceptor vector is the minimum image of acceptor minus donor. The acceptor-hydrogen vector is the minimum image of acceptor minus the hydrogen. A linear donor-hydrogen-acceptor geometry has angle 0 degrees. The acceptor-hydrogen distance is kept only when it is strictly below dist_cutoff. The angle is kept when it is less than or equal to angle_cutoff. The pair list keeps a donor and acceptor only when their minimum image is strictly below pair_cutoff, on an orthorhombic or triclinic cell. For a triclinic cell, every face separation must be at least twice pair_cutoff; a thinner cell can miss valid pairs. Each hydrogen is assigned to the nearest heavy atom that shares its molID, then kept only when that atom's type is in donor_types. The lower index wins a distance tie. exclude_same_mol of 1 drops a pair that shares a molID. An empty donor set or acceptor set writes zero edges. Atoms whose type equals hydrogen_type are hydrogens and are never donors or acceptors. Water-water and water-framework counts are two calls with different type sets. When edges is non-NULL, every element must start with size set to sizeof(seams_directed_hbond). @par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller owns sys, the type arrays, edges, and n_out. The call borrows them and writes caller buffers only. @par Thread safety Exclusive use of sys. @par Example @code const int donors[] = {8}; const int acceptors[] = {10}; size_t n_out = 0; seams_directed_hbonds(sys, 3.5, 2.42, 30.0, donors, 1, acceptors, 1, 2, 1, NULL, 0, &n_out); @endcode

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

.. _capi-seams_mols:

seams_mols
----------


.. code-block:: c

   seams_status seams_mols(const seams_system *sys, int *mols, size_t len)

Copy the stored molecule ids out.

``sys``
    Frame handle.

``mols``
    Caller buffer for n ints, in stored order.

``len``
    Length of mols. Must be at least n.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE for a short buffer.

@par Ownership The caller owns mols. @par Thread safety Read-only: may run alongside other const calls on sys. @par Example @code int mols[8] = {0}; seams_mols(sys, mols, 8); @endcode

.. _capi-seams_set_mols:

seams_set_mols
--------------


.. code-block:: c

   seams_status seams_set_mols(seams_system *sys, const int *mols, size_t len)

Replace the stored molecule ids on a loaded frame.

``sys``
    Frame handle.

``mols``
    Caller buffer of n molecule ids.

``len``
    Length of mols. Must equal the stored atom count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE when no frame is stored; ::SEAMS_ERR_RANGE when len differs from the atom count.

@par Handle access Non-const: updates stored molecule ids and clears the handle error slot. @par Ownership The library copies mols and keeps no pointer to caller memory. @par Thread safety Exclusive use of sys. @par Example @code const int mols[8] = {1, 1, 2, 2, 3, 3, 4, 4}; seams_set_mols(sys, mols, 8); @endcode

.. _capi-seams_molecule_centres:

seams_molecule_centres
----------------------


.. code-block:: c

   seams_status seams_molecule_centres(const seams_system *sys, const int *mols, size_t mols_len, const double *masses, size_t masses_len, int *mol_ids_out, double *centres_out, size_t cap, size_t *n_mols)

Intact periodic centres of mass or geometry for molecules in a frame.

``sys``
    Frame handle.

``mols``
    Optional molecule id per atom, or NULL to use stored molecule ids.

``mols_len``
    Length of mols when mols is non-NULL; must equal n, or 0 when NULL.

``masses``
    Optional positive weight per atom, or NULL for equal weights.

``masses_len``
    Length of masses when masses is non-NULL; must equal n, or 0 when NULL.

``mol_ids_out``
    Receives one molecule id per molecule, or NULL.

``centres_out``
    Receives 3 \* n_mols coordinates in Angstrom, or NULL.

``cap``
    Molecule capacity of non-NULL output buffers.

``n_mols``
    Receives the molecule count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.

Atoms sharing a molecule id are assembled onto the first atom of that molecule through the minimum image of the stored orthorhombic or triclinic cell before averaging. Pass NULL for both mol_ids_out and centres_out to query the molecule count in n_mols.

@par Ownership The caller owns all buffers. @par Thread safety Read-only on sys. @par Example @code size_t nm = 0; seams_molecule_centres(sys, NULL, 0, NULL, 0, NULL, NULL, 0, &nm); @endcode

.. _capi-seams_traj_unwrap_positions:

seams_traj_unwrap_positions
---------------------------


.. code-block:: c

   seams_status seams_traj_unwrap_positions(const double *xyz, size_t nframes, size_t natoms, const double *boxes, size_t nboxes, double *unwrapped, size_t unwrapped_len)

Unwrap periodic coordinates across consecutive frames.

``xyz``
    Wrapped coordinates, nframes \* 3 \* natoms doubles in Angstrom.

``nframes``
    Number of frames, at least 1.

``natoms``
    Number of atoms per frame, at least 1.

``boxes``
    Cell array of 9 \* nboxes doubles.

``nboxes``
    Cell count, either 1 or nframes.

``unwrapped``
    Output buffer of unwrapped_len doubles.

``unwrapped_len``
    Capacity in doubles, at least nframes \* 3 \* natoms.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

Each step adds the minimum-image displacement under an orthorhombic or triclinic nine-number cell (ox, oy, oz, lx, ly, lz, xy, xz, yz). Pass nboxes equal to 1 for a constant cell or nframes for a deforming cell.

@par Ownership The caller owns xyz, boxes, and unwrapped. @par Thread safety Pure aside from the thread-local error string. @par Example @code const double box[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0}; const double wrapped[6] = {9.5, 1.0, 1.0, 0.5, 1.0, 1.0}; double unwrapped[6] = {0}; seams_traj_unwrap_positions(wrapped, 2, 1, box, 1, unwrapped, 6); @endcode

.. _capi-seams_traj_unwrap_centres:

seams_traj_unwrap_centres
-------------------------


.. code-block:: c

   seams_status seams_traj_unwrap_centres( const double *xyz, size_t nframes, size_t natoms, const double *boxes, size_t nboxes, const int *mols, size_t mols_len, const double *masses, size_t masses_len, int *mol_ids_out, size_t mol_ids_cap, double *centres_out, size_t centres_len, size_t *n_mols)

Assemble per-molecule centres and unwrap them across frames.

``xyz``
    Wrapped atom coordinates, nframes \* 3 \* natoms doubles.

``nframes``
    Number of frames, at least 1.

``natoms``
    Number of atoms per frame, at least 1.

``boxes``
    Cell array of 9 \* nboxes doubles.

``nboxes``
    Cell count, either 1 or nframes.

``mols``
    Molecule id per atom, length mols_len.

``mols_len``
    Must equal natoms.

``masses``
    Optional positive mass per atom, or NULL for equal weights.

``masses_len``
    Must equal natoms when masses is non-NULL, or 0 when NULL.

``mol_ids_out``
    Receives unique molecule ids, or NULL.

``mol_ids_cap``
    Capacity of mol_ids_out in ints when non-NULL.

``centres_out``
    Receives nframes \* 3 \* n_mols unwrapped coordinates, or NULL.

``centres_len``
    Capacity of centres_out in doubles when non-NULL.

``n_mols``
    Receives the molecule count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

Pass NULL for both mol_ids_out and centres_out to query n_mols first.

@par Ownership The caller owns all buffers. @par Thread safety Pure aside from the thread-local error string. @par Example @code const double box[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0}; const double xyz[12] = {9.8, 0, 0, 0.2, 0, 0, 0.4, 0, 0, 0.8, 0, 0}; const int mols[2] = {1, 1}; size_t nm = 0; double out[6] = {0}; seams_traj_unwrap_centres(xyz, 2, 2, box, 1, mols, 2, NULL, 0, NULL, 0, out, 6, &nm); @endcode

.. _capi-seams_water_vectors:

seams_water_vectors
-------------------


.. code-block:: c

   seams_status seams_water_vectors(const seams_system *sys, int oxygen_type, int hydrogen_type, int *mol_ids, double *oh_mid, double *bisector, double *hh, double *hh_unit, size_t cap, size_t *n_waters)

Water orientation vectors: O to H-midpoint and H-to-H under PBC.

``sys``
    Frame handle.

``oxygen_type``
    Species integer of water oxygen.

``hydrogen_type``
    Species integer of water hydrogen, distinct from oxygen_type.

``mol_ids``
    Receives molecule ids, length cap, or NULL.

``oh_mid``
    Receives O-to-H-midpoint vectors in Angstrom, 3 \* cap doubles, or NULL.

``bisector``
    Receives unit bisector vectors, 3 \* cap doubles, or NULL.

``hh``
    Receives H1-to-H2 vectors in Angstrom, 3 \* cap doubles, or NULL.

``hh_unit``
    Receives unit H1-to-H2 vectors, 3 \* cap doubles, or NULL.

``cap``
    Molecule capacity of every non-NULL output array.

``n_waters``
    Receives the water count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.

Each water pairs one oxygen_type atom with two hydrogen_type atoms via shared molecule ids when available, or via minimum-image covalent distance otherwise. Pass NULL for all five output arrays to query n_waters.

@par Ownership The caller owns all output arrays. @par Thread safety Read-only on sys. @par Example @code size_t nw = 0; seams_water_vectors(sys, 1, 2, NULL, NULL, NULL, NULL, NULL, 0, &nw); @endcode

.. _capi-seams_molecule_dipoles:

seams_molecule_dipoles
----------------------


.. code-block:: c

   seams_status seams_molecule_dipoles( const seams_system *sys, const int *mols, size_t mols_len, const double *charges, size_t charges_len, int *mol_ids_out, double *dipoles, double *unit_vectors, size_t cap, size_t *n_mols)

Generic per-molecule dipole or geometric orientation vectors.

``sys``
    Frame handle.

``mols``
    Optional molecule id per atom, or NULL to use stored molecule ids.

``mols_len``
    Length of mols when non-NULL, or 0 when NULL.

``charges``
    Optional partial charge per atom, or NULL for geometric mode.

``charges_len``
    Length of charges when non-NULL, or 0 when NULL.

``mol_ids_out``
    Receives molecule ids, length cap, or NULL.

``dipoles``
    Receives 3 \* cap raw dipole vectors, or NULL.

``unit_vectors``
    Receives 3 \* cap normalized orientation vectors, or NULL.

``cap``
    Molecule capacity of non-NULL output arrays.

``n_mols``
    Receives the molecule count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_BOX; ::SEAMS_ERR_NOMEM.

With charges, computes sum_a q_a (r_a - r_centre) in elementary-charge Angstrom on each PBC-assembled molecule. Without charges, computes the mean minimum-image vector from each molecule's first atom to its remaining atoms.

@par Ownership The caller owns all buffers. @par Thread safety Read-only on sys. @par Example @code size_t nm = 0; seams_molecule_dipoles(sys, NULL, 0, NULL, 0, NULL, NULL, NULL, 0, &nm); @endcode

.. _capi-seams_orientation_order:

seams_orientation_order
-----------------------


.. code-block:: c

   seams_status seams_orientation_order(const double *vectors, size_t n_vectors, const double axis[3], double *cos_theta, double *p2, size_t cap, double *mean_p1, double *mean_p2)

Project molecular vectors onto an axis for P1 and P2 order parameters.

``vectors``
    Input 3-vectors, 3 \* n_vectors doubles (or NULL when n_vectors is 0).

``n_vectors``
    Number of 3-vectors.

``axis``
    Non-zero reference axis of 3 finite doubles.

``cos_theta``
    Optional output buffer of cap doubles for cos(theta), or NULL.

``p2``
    Optional output buffer of cap doubles for P2(cos(theta)), or NULL.

``cap``
    Capacity of non-NULL per-vector buffers; must be at least n_vectors.

``mean_p1``
    Receives ensemble mean of cos(theta), or NULL.

``mean_p2``
    Receives ensemble mean of P2(cos(theta)), or NULL.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

Computes cos_theta = u . e and P2 = 0.5 \* (3 cos_theta^2 - 1) for each non-zero vector, where e is the normalized axis. Zero-length vectors receive NaN in per-molecule outputs and are omitted from mean_p1 and mean_p2.

@par Ownership The caller owns all arrays. @par Thread safety Pure aside from the thread-local error string. @par Example @code const double vecs[6] = {0, 0, 1, 0, 0, 2}; const double axis[3] = {0, 0, 1}; double p1 = 0, p2 = 0; seams_orientation_order(vecs, 2, axis, NULL, NULL, 0, &p1, &p2); @endcode

.. _capi-seams_msd_create:

seams_msd_create
----------------


.. code-block:: c

   seams_status seams_msd_create(size_t natoms, size_t max_lag, size_t origin_stride, double dt, const double axis[3], seams_msd **out)

Create a streaming multi-origin MSD accumulator with axial/lateral split.

``natoms``
    Particle or molecule count per frame, at least 1. Above PTRDIFF_MAX / (3 \* sizeof(double)) the coordinates cannot be allocated, which is ::SEAMS_ERR_RANGE.

``max_lag``
    Maximum lag in frames, or 0 to retain all lags. Above PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated, which is ::SEAMS_ERR_RANGE.

``origin_stride``
    Frame stride between time origins, at least 1.

``dt``
    Positive finite frame spacing in physical time units.

``axis``
    Non-zero 3-vector defining the axial direction.

``out``
    Receives the new handle; set to NULL on failure.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

@par Ownership The caller frees \*out with seams_msd_destroy(). @par Thread safety Pure aside from handle registration and the thread-local error string. @par Example @code const double axis[3] = {0, 0, 1}; seams_msd \*msd = NULL; seams_msd_create(4, 16, 1, 0.1, axis, &msd); seams_msd_destroy(msd); @endcode

.. _capi-seams_msd_destroy:

seams_msd_destroy
-----------------


.. code-block:: c

   void seams_msd_destroy(seams_msd *msd)

Free an MSD accumulator.

``msd``
    Handle from seams_msd_create(), or NULL.

@par Ownership Releases the handle and its ring buffer. @par Thread safety Exclusive use of msd. @par Example @code seams_msd_destroy(NULL); @endcode

.. _capi-seams_msd_push:

seams_msd_push
--------------


.. code-block:: c

   seams_status seams_msd_push(seams_msd *msd, const double *xyz, size_t natoms, const double *box9, size_t box_len)

Push one frame into a streaming MSD accumulator.

``msd``
    MSD handle.

``xyz``
    Coordinates of 3 \* natoms doubles in Angstrom.

``natoms``
    Particle count; must match the count passed to seams_msd_create().

``box9``
    Optional 9-number periodic cell, or NULL for unwrapped input.

``box_len``
    Must be 9 when box9 is non-NULL, or 0 when box9 is NULL.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_BOX; ::SEAMS_ERR_ARG; ::SEAMS_ERR_STATE; ::SEAMS_ERR_NOMEM.

Pass box9 as NULL (with box_len 0) when xyz is already unwrapped, or pass a nine-number cell (with box_len 9) to unwrap incrementally across frames.

@par Ownership The library copies coordinates into active origin slots as needed. @par Thread safety Exclusive use of msd. @par Example @code const double axis[3] = {0, 0, 1}; const double f0[3] = {0, 0, 0}; seams_msd \*msd = NULL; seams_msd_create(1, 8, 1, 1.0, axis, &msd); seams_msd_push(msd, f0, 1, NULL, 0); seams_msd_destroy(msd); @endcode

.. _capi-seams_msd_curve:

seams_msd_curve
---------------


.. code-block:: c

   seams_status seams_msd_curve(const seams_msd *msd, double *lag_time, double *total, double *axial, double *lateral, size_t *n_origins, size_t cap, size_t *n_lags)

Read accumulated lag times and total, axial, and lateral MSD curves.

``msd``
    MSD handle.

``lag_time``
    Receives lag times k \* dt, or NULL.

``total``
    Receives 3D total MSD, or NULL.

``axial``
    Receives 1D axial MSD along axis, or NULL.

``lateral``
    Receives 2D lateral MSD perpendicular to axis, or NULL.

``n_origins``
    Receives time-origin counts per lag, or NULL.

``cap``
    Capacity of non-NULL output arrays.

``n_lags``
    Receives the available lag count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_NOMEM.

Any of lag_time, total, axial, lateral, and n_origins may be NULL to skip that column; pass all five as NULL to query n_lags.

@par Ownership The caller owns all output buffers. @par Thread safety Read-only on msd. @par Example @code const double axis[3] = {0, 0, 1}; const double f0[3] = {0, 0, 0}, f1[3] = {0, 0, 1}; seams_msd \*msd = NULL; size_t nlags = 0; seams_msd_create(1, 4, 1, 1.0, axis, &msd); seams_msd_push(msd, f0, 1, NULL, 0); seams_msd_push(msd, f1, 1, NULL, 0); seams_msd_curve(msd, NULL, NULL, NULL, NULL, NULL, 0, &nlags); seams_msd_destroy(msd); @endcode

.. _capi-seams_msd_fit:

seams_msd_fit
-------------


.. code-block:: c

   seams_status seams_msd_fit(const seams_msd *msd, double t_start, double t_end, double *d_total, double *d_axial, double *d_lateral, double *alpha)

Linear Einstein diffusion fit on an MSD accumulator.

``msd``
    MSD handle with at least 2 lag points.

``t_start``
    Window start time, >= 0.

``t_end``
    Window end time, > t_start (or both 0 for default window).

``d_total``
    Receives 3D diffusion coefficient, or NULL.

``d_axial``
    Receives 1D axial diffusion coefficient, or NULL.

``d_lateral``
    Receives 2D lateral diffusion coefficient, or NULL.

``alpha``
    Receives the log-log anomalous exponent, or NULL. Receives NaN when the window holds fewer than two lags with t > 0 and MSD_total > 1e-20, so the log-log fit was skipped; test it with isnan(). The diffusion coefficients are still the linear fit.

Returns: ::SEAMS_OK when the linear fit succeeded, including a skipped log-log fit (alpha NaN); ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

Fits MSD_total = 6 D_total t + c, MSD_axial = 2 D_axial t + c, and MSD_lateral = 4 D_lateral t + c over [t_start, t_end]. Pass t_start = 0 and t_end = 0 for the default [0.1 t_max, 0.5 t_max] window.

@par Ownership Writes scalar outputs owned by the caller. @par Thread safety Read-only on msd. @par Example @code const double axis[3] = {0, 0, 1}; const double f0[3] = {0, 0, 0}, f1[3] = {0, 0, 1}, f2[3] = {0, 0, 2}; seams_msd \*msd = NULL; double dtot = 0, dax = 0, dlat = 0, alpha = 0; seams_msd_create(1, 4, 1, 1.0, axis, &msd); seams_msd_push(msd, f0, 1, NULL, 0); seams_msd_push(msd, f1, 1, NULL, 0); seams_msd_push(msd, f2, 1, NULL, 0); seams_msd_fit(msd, 0.0, 0.0, &dtot, &dax, &dlat, &alpha); seams_msd_destroy(msd); @endcode

.. _capi-seams_hbond_tcf_create:

seams_hbond_tcf_create
----------------------


.. code-block:: c

   seams_status seams_hbond_tcf_create(size_t max_lag, size_t origin_stride, double dt, int directed, seams_hbond_tcf **out)

Create a streaming H-bond lifetime accumulator from per-frame edges.

``max_lag``
    Maximum lag in frames, or 0 to retain all lags. Above PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated, which is ::SEAMS_ERR_RANGE.

``origin_stride``
    Frame stride between time origins, at least 1.

``dt``
    Positive finite frame spacing.

``directed``
    Non-zero treats (donor, acceptor) as directed; zero canonicalizes pairs.

``out``
    Receives the new handle; set to NULL on failure.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

@par Ownership The caller frees \*out with seams_hbond_tcf_destroy(). @par Thread safety Pure aside from handle registration and the thread-local error string. @par Example @code seams_hbond_tcf \*tcf = NULL; seams_hbond_tcf_create(16, 1, 0.1, 1, &tcf); seams_hbond_tcf_destroy(tcf); @endcode

.. _capi-seams_hbond_tcf_destroy:

seams_hbond_tcf_destroy
-----------------------


.. code-block:: c

   void seams_hbond_tcf_destroy(seams_hbond_tcf *tcf)

Free an H-bond lifetime accumulator.

``tcf``
    Handle from seams_hbond_tcf_create(), or NULL.

@par Ownership Releases the handle and its active origin sets. @par Thread safety Exclusive use of tcf. @par Example @code seams_hbond_tcf_destroy(NULL); @endcode

.. _capi-seams_hbond_tcf_push:

seams_hbond_tcf_push
--------------------


.. code-block:: c

   seams_status seams_hbond_tcf_push(seams_hbond_tcf *tcf, const int *donors, const int *acceptors, size_t n_edges)

Push one frame's hydrogen-bond edge set into the lifetime accumulator.

``tcf``
    H-bond TCF handle.

``donors``
    Donor ids of length n_edges (may be NULL when n_edges is 0).

``acceptors``
    Acceptor ids of length n_edges (may be NULL when n_edges is 0).

``n_edges``
    Number of edges in this frame.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_NOMEM.

Duplicate edges within the same frame are deduplicated before updating intermittent C_HB(t) and continuous survival S_HB(t).

@par Ownership The library copies edge ids and keeps no caller pointer. @par Thread safety Exclusive use of tcf. @par Example @code const int don[1] = {1}, acc[1] = {2}; seams_hbond_tcf \*tcf = NULL; seams_hbond_tcf_create(8, 1, 1.0, 1, &tcf); seams_hbond_tcf_push(tcf, don, acc, 1); seams_hbond_tcf_destroy(tcf); @endcode

.. _capi-seams_hbond_tcf_curve:

seams_hbond_tcf_curve
---------------------


.. code-block:: c

   seams_status seams_hbond_tcf_curve(const seams_hbond_tcf *tcf, double *lag_time, double *c_intermittent, double *s_continuous, size_t *n_origins, size_t cap, size_t *n_lags)

Read intermittent C_HB(t) and continuous S_HB(t) lifetime curves.

``tcf``
    H-bond TCF handle.

``lag_time``
    Receives lag times k \* dt, or NULL.

``c_intermittent``
    Receives intermittent autocorrelation C_HB, or NULL.

``s_continuous``
    Receives unbroken survival probability S_HB, or NULL.

``n_origins``
    Receives origin counts per lag, or NULL.

``cap``
    Capacity of non-NULL output arrays.

``n_lags``
    Receives the available lag count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_NOMEM.

@par Ownership The caller owns all output buffers. @par Thread safety Read-only on tcf. @par Example @code const int don[1] = {1}, acc[1] = {2}; seams_hbond_tcf \*tcf = NULL; size_t nlags = 0; seams_hbond_tcf_create(4, 1, 1.0, 1, &tcf); seams_hbond_tcf_push(tcf, don, acc, 1); seams_hbond_tcf_push(tcf, don, acc, 1); seams_hbond_tcf_curve(tcf, NULL, NULL, NULL, NULL, 0, &nlags); seams_hbond_tcf_destroy(tcf); @endcode

.. _capi-seams_hbond_tcf_fit:

seams_hbond_tcf_fit
-------------------


.. code-block:: c

   seams_status seams_hbond_tcf_fit(const seams_hbond_tcf *tcf, double t_start, double t_end, double *tau_c_int, double *tau_c_exp, double *tau_s_int, double *tau_s_exp)

Integrated and single-exponential lifetimes for C_HB and S_HB.

``tcf``
    H-bond TCF handle with at least 2 lag points.

``t_start``
    Fit window start time, >= 0.

``t_end``
    Fit window end time, > t_start (or both 0 for full range).

``tau_c_int``
    Receives trapezoidal integral of C_HB, or NULL.

``tau_c_exp``
    Receives exponential decay time of C_HB, or NULL.

``tau_s_int``
    Receives trapezoidal integral of S_HB, or NULL.

``tau_s_exp``
    Receives exponential decay time of S_HB, or NULL.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

Pass t_start = 0 and t_end = 0 to fit across the full positive decay range.

@par Ownership Writes scalar outputs owned by the caller. @par Thread safety Read-only on tcf. @par Example @code const int don[2] = {1, 2}, acc[2] = {3, 4}; seams_hbond_tcf \*tcf = NULL; double tc = 0, ts = 0; seams_hbond_tcf_create(4, 1, 1.0, 1, &tcf); seams_hbond_tcf_push(tcf, don, acc, 2); seams_hbond_tcf_push(tcf, don, acc, 1); seams_hbond_tcf_fit(tcf, 0.0, 0.0, &tc, NULL, &ts, NULL); seams_hbond_tcf_destroy(tcf); @endcode

.. _capi-seams_reorient_create:

seams_reorient_create
---------------------


.. code-block:: c

   seams_status seams_reorient_create(size_t n_vectors, size_t max_lag, size_t origin_stride, double dt, seams_reorient **out)

Create a streaming orientational correlation accumulator for C1 and C2.

``n_vectors``
    Number of molecular vectors per frame, at least 1. Above PTRDIFF_MAX / (3 \* sizeof(double)) the vectors cannot be allocated, which is ::SEAMS_ERR_RANGE.

``max_lag``
    Maximum lag in frames, or 0 to retain all lags. Above PTRDIFF_MAX / sizeof(double) - 1 the lag bins cannot be allocated, which is ::SEAMS_ERR_RANGE.

``origin_stride``
    Frame stride between time origins, at least 1.

``dt``
    Positive finite frame spacing.

``out``
    Receives the new handle; set to NULL on failure.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

@par Ownership The caller frees \*out with seams_reorient_destroy(). @par Thread safety Pure aside from handle registration and the thread-local error string. @par Example @code seams_reorient \*rot = NULL; seams_reorient_create(4, 16, 1, 0.1, &rot); seams_reorient_destroy(rot); @endcode

.. _capi-seams_reorient_destroy:

seams_reorient_destroy
----------------------


.. code-block:: c

   void seams_reorient_destroy(seams_reorient *reorient)

Free a reorientation accumulator.

``reorient``
    Handle from seams_reorient_create(), or NULL.

@par Ownership Releases the handle and its ring buffer. @par Thread safety Exclusive use of reorient. @par Example @code seams_reorient_destroy(NULL); @endcode

.. _capi-seams_reorient_push:

seams_reorient_push
-------------------


.. code-block:: c

   seams_status seams_reorient_push(seams_reorient *reorient, const double *vectors, size_t n_vectors)

Push one frame of molecular 3-vectors into a reorientation accumulator.

``reorient``
    Reorientation handle.

``vectors``
    Array of 3 \* n_vectors finite doubles.

``n_vectors``
    Vector count; must match seams_reorient_create().

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

Non-zero vectors are normalized per molecule; zero-norm vectors are skipped in correlation pairs.

@par Ownership The library copies normalized vectors into active origin slots. @par Thread safety Exclusive use of reorient. @par Example @code const double v0[3] = {0, 0, 1}; seams_reorient \*rot = NULL; seams_reorient_create(1, 4, 1, 1.0, &rot); seams_reorient_push(rot, v0, 1); seams_reorient_destroy(rot); @endcode

.. _capi-seams_reorient_curve:

seams_reorient_curve
--------------------


.. code-block:: c

   seams_status seams_reorient_curve(const seams_reorient *reorient, double *lag_time, double *c1, double *c2, size_t *n_origins, size_t cap, size_t *n_lags)

Read accumulated C1(t) = <P1(u(0).u(t))> and C2(t) = <P2(u(0).u(t))>.

``reorient``
    Reorientation handle.

``lag_time``
    Receives lag times k \* dt, or NULL.

``c1``
    Receives first Legendre correlation C1, or NULL.

``c2``
    Receives second Legendre correlation C2, or NULL.

``n_origins``
    Receives origin counts per lag, or NULL.

``cap``
    Capacity of non-NULL output arrays.

``n_lags``
    Receives the available lag count.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_NOMEM.

@par Ownership The caller owns all output buffers. @par Thread safety Read-only on reorient. @par Example @code const double v0[3] = {0, 0, 1}, v1[3] = {0, 1, 0}; seams_reorient \*rot = NULL; size_t nlags = 0; seams_reorient_create(1, 4, 1, 1.0, &rot); seams_reorient_push(rot, v0, 1); seams_reorient_push(rot, v1, 1); seams_reorient_curve(rot, NULL, NULL, NULL, NULL, 0, &nlags); seams_reorient_destroy(rot); @endcode

.. _capi-seams_reorient_fit:

seams_reorient_fit
------------------


.. code-block:: c

   seams_status seams_reorient_fit(const seams_reorient *reorient, double t_start, double t_end, double *tau1_int, double *tau1_exp, double *tau2_int, double *tau2_exp)

Integrated and single-exponential reorientation times tau_1 and tau_2.

``reorient``
    Reorientation handle with at least 2 lag points.

``t_start``
    Fit window start time, >= 0.

``t_end``
    Fit window end time, > t_start (or both 0 for full range).

``tau1_int``
    Receives trapezoidal integral of C1, or NULL.

``tau1_exp``
    Receives exponential decay time of C1, or NULL.

``tau2_int``
    Receives trapezoidal integral of C2, or NULL.

``tau2_exp``
    Receives exponential decay time of C2, or NULL.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL; ::SEAMS_ERR_STATE; ::SEAMS_ERR_RANGE; ::SEAMS_ERR_ARG; ::SEAMS_ERR_NOMEM.

@par Ownership Writes scalar outputs owned by the caller. @par Thread safety Read-only on reorient. @par Example @code const double v0[3] = {0, 0, 1}, v1[3] = {0.6, 0, 0.8}; seams_reorient \*rot = NULL; double t1 = 0, t2 = 0; seams_reorient_create(1, 4, 1, 1.0, &rot); seams_reorient_push(rot, v0, 1); seams_reorient_push(rot, v1, 1); seams_reorient_fit(rot, 0.0, 0.0, &t1, NULL, &t2, NULL); seams_reorient_destroy(rot); @endcode

