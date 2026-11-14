Further analyses
================


Brine fronts, topology, and the remaining exports.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 69 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

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

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

ice[i] is non-zero when atom i is ice. water_type 0 keeps every atom whose type is not ion_type. ion_type 0 means the frame has no ions. n_bins 0 selects 20 bins. axis is 0, 1, or 2. The ice fraction cut is 0.5, and the bulk ice and liquid cuts are 0.85 and 0.15.

Ownership: the caller frees \*out with seams_brine_destroy(). Thread safety: exclusive use of sys.

.. _capi-seams_brine_destroy:

seams_brine_destroy
-------------------


.. code-block:: c

   void seams_brine_destroy(seams_brine *front)

Free a brine-front handle.

``front``
    Handle from seams_brine_analyse(), or NULL.

Thread safety: the caller must not use front after this returns.

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

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL. At least one output pointer must be non-NULL.

Lengths are Angstrom, in the frame's cell. Thread safety: exclusive use of front.

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

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL. At least one output pointer must be non-NULL.

Molality is formula units per kilogram of water. Thread safety: exclusive use of front.

.. _capi-seams_brine_profile:

seams_brine_profile
-------------------


.. code-block:: c

   seams_status seams_brine_profile(const seams_brine *front, double *z_mid, double *molality, double *ice_fraction, size_t cap, size_t *n_out)

Concentration profile along the front axis.

``front``
    Brine handle.

``z_mid``
    Bin centres, or NULL to skip the column.

``molality``
    Bin molalities, or NULL to skip the column.

``ice_fraction``
    Bin ice fractions, or NULL to skip the column.

``cap``
    Capacity of each non-NULL column.

``n_out``
    Receives the bin count.

Returns: SEAMS_OK when every non-NULL column was filled, or when every column pointer was NULL and n_out was written. SEAMS_ERR_RANGE when a non-NULL column is short. SEAMS_ERR_NULL when front or n_out is NULL.

z_mid is Angstrom. molality is mol/kg. ice_fraction is dimensionless. Ownership: the caller owns the columns. Thread safety: exclusive use of front.

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

Returns: SEAMS_OK, or SEAMS_ERR_NULL when front is NULL. At least one output pointer must be non-NULL.

Thread safety: exclusive use of front.

.. _capi-seams_rdf:

seams_rdf
---------


.. code-block:: c

   seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *g, int *count, size_t cap, size_t *n_out, double *volume)

Partial radial distribution.

``bins``
    Bin count. Zero selects about rmax / 0.1.

``r,``
    g, count Caller columns, or NULL to query.

``n_out``
    Receives the bin count.

``volume``
    Receives the cell volume, or NULL.

r is Angstrom. g is dimensionless. volume is cubic Angstrom. Thread safety: exclusive use of sys.

.. _capi-seams_running_cn:

seams_running_cn
----------------


.. code-block:: c

   seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *cn, size_t cap, size_t *n_out)

Running coordination number along the partial RDF.

Thread safety: exclusive use of sys.

.. _capi-seams_coordination:

seams_coordination
------------------


.. code-block:: c

   seams_status seams_coordination(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *cn)

Coordination number integrated to rmax.

Thread safety: exclusive use of sys.

.. _capi-seams_hbonds:

seams_hbonds
------------


.. code-block:: c

   seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int *n_bonds)

Hydrogen-bond count on one frame.

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

Thread safety: exclusive use of sys.

.. _capi-seams_density_axis:

seams_density_axis
------------------


.. code-block:: c

   seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis, double *z, double *rho, size_t cap, size_t *n_out)

Number density along one cell axis.

``bins``
    Bin count. Zero selects about one bin per 0.1 Angstrom.

``axis``
    0, 1, or 2.

z is Angstrom. rho is number density. Thread safety: exclusive use of sys.

.. _capi-seams_fingerprint:

seams_fingerprint
-----------------


.. code-block:: c

   seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k, int graph, int hops, int colour_types, seams_topo **out)

Topological fingerprint of the bonded graph.

``hops``
    Neighbourhood depth. Zero selects 2.

``colour_types``
    Non-zero colours vertices by LAMMPS type.

``out``
    Receives a seams_topo. Free it with seams_fingerprint_destroy().

Thread safety: exclusive use of sys.

.. _capi-seams_fingerprint_destroy:

seams_fingerprint_destroy
-------------------------


.. code-block:: c

   void seams_fingerprint_destroy(seams_topo *fp)

@brief Free a fingerprint. NULL is ignored.

.. _capi-seams_fingerprint_key:

seams_fingerprint_key
---------------------


.. code-block:: c

   seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap, size_t *n_out)

Frame key text. A null buffer queries \*n_out, the length without the terminating NUL. A short buffer returns SEAMS_ERR_RANGE.

.. _capi-seams_fingerprint_atom:

seams_fingerprint_atom
----------------------


.. code-block:: c

   seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf, size_t cap, size_t *n_out)

@brief Local key of one atom. Same buffer rule as seams_fingerprint_key().

.. _capi-seams_fingerprint_rings:

seams_fingerprint_rings
-----------------------


.. code-block:: c

   seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts, size_t counts_len)

@brief Primitive-ring census stored on the fingerprint. counts[k] is size k.

.. _capi-seams_contact_pairs:

seams_contact_pairs
-------------------


.. code-block:: c

   seams_status seams_contact_pairs(seams_system *sys, const char *site_spec, size_t *n_pairs, int *n_cation, int *n_anion)

Mutual nearest unlike contacts of a site table.

``site_spec``
    Text such as "1=cationHead,2=anion".

Cation and anion counts are molecules in the ion cloud. Thread safety: exclusive use of sys.

.. _capi-seams_domains:

seams_domains
-------------


.. code-block:: c

   seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec, int subset, int *n_subset, int *largest, double *percolation)

Largest domain of a polar or apolar site subset.

``subset``
    0 is apolar, 1 is polar.

Thread safety: exclusive use of sys.

.. _capi-seams_site_table:

seams_site_table
----------------


.. code-block:: c

   seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds, size_t len)

Per-atom site kind from a site-table spec.

``kinds``
    One seams_site_kind per atom.

Thread safety: exclusive use of sys.

