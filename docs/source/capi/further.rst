Further analyses
================


Brine fronts, topology, and the remaining exports.

Generated from the Doxygen XML of ``include/seams.h`` at ``capi/polish-2.14`` commit ``0127f41efc3f44cefdcf532bdea3bba36b05176d``. The C sources of that branch are not part of this docs tree.

.. _capi-seams_brine_analyse:

seams_brine_analyse
-------------------


.. code-block:: c

   seams_status seams_brine_analyse(seams_system *sys, const int *ice, size_t n_ice, int water_type, int ion_type, int n_bins, int axis, seams_brine **out)

Planar ice/brine front from an ice assignment.

ice[i] is non-zero when atom i is ice. water_type 0 keeps every atom whose type is not ion_type. ion_type 0 means the frame has no ions. n_bins 0 selects 20 bins. axis is 0, 1, or 2. The ice fraction cut is 0.5, and the bulk ice and liquid cuts are 0.85 and 0.15. sys Frame handle. ice Per-atom ice flags, length n_ice. n_ice Must equal the frame count. water_type Species kept as water, or 0. ion_type Species kept as ions, or 0 for none. n_bins Bin count along axis, or 0 for 20. axis 0, 1, or 2. out Receives a new handle. On failure \*out is NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller frees \*out with seams_brine_destroy(). Thread safety Exclusive use of sys. Example seams_brine\*front=NULL; if(seams_brine_analyse(sys,ice,n,1,2,0,2,&front)==SEAMS_OK){ seams_brine_destroy(front); }

.. _capi-seams_brine_destroy:

seams_brine_destroy
-------------------


.. code-block:: c

   void seams_brine_destroy(seams_brine *front)

Free a brine-front handle.

front Handle from seams_brine_analyse(), or NULL. NULL is ignored. Ownership Releases the handle and its profile. Thread safety The caller must not use front after this returns. Example seams_brine_destroy(front);

.. _capi-seams_brine_fronts:

seams_brine_fronts
------------------


.. code-block:: c

   seams_status seams_brine_fronts(const seams_brine *front, double *z_enter, double *z_exit, double *thickness, int *found)

The two planar crossings of the ice slab.

Lengths are Angstrom, in the frame's cell. front Brine handle. z_enter Where the ice fraction crosses the cut on entering the ice along the positive axis, or NULL. z_exit The crossing on leaving the ice, or NULL. thickness Ice run along the axis, or NULL. found Receives 1 when both crossings exist, otherwise 0, or NULL. SEAMS_OK, or SEAMS_ERR_NULL when front is NULL or every output pointer is NULL. Thread safety Exclusive use of front. Example doublethick=0; intfound=0; seams_brine_fronts(front,NULL,NULL,&thick,&found);

.. _capi-seams_brine_pockets:

seams_brine_pockets
-------------------


.. code-block:: c

   seams_status seams_brine_pockets(const seams_brine *front, int *n_pockets, int *reservoir_waters)

Brine pockets trapped in the ice, excluding the liquid reservoir.

front Brine handle. n_pockets Receives the pocket count, or NULL. reservoir_waters Receives the water count of the liquid reservoir, or NULL. SEAMS_OK, or SEAMS_ERR_NULL when front is NULL or both output pointers are NULL. Thread safety Exclusive use of front. Example intpockets=0; seams_brine_pockets(front,&pockets,NULL);

.. _capi-seams_brine_profile:

seams_brine_profile
-------------------


.. code-block:: c

   seams_status seams_brine_profile(const seams_brine *front, double *z_mid, double *molality, double *ice_fraction, size_t cap, size_t *n_out)

Concentration profile along the front axis.

z_mid is Angstrom. molality is mol/kg. ice_fraction is dimensionless. front Brine handle. z_mid Bin centres, or NULL to skip the column. molality Bin molalities, or NULL to skip the column. ice_fraction Bin ice fractions, or NULL to skip the column. cap Capacity of each non-NULL column. n_out Receives the bin count. SEAMS_OK when every non-NULL column was filled, or when every column pointer was NULL and n_out was written. SEAMS_ERR_RANGE when a non-NULL column is short. SEAMS_ERR_NULL when front or n_out is NULL. Ownership The caller owns the columns. Thread safety Exclusive use of front. Example size_tbins=0; seams_brine_profile(front,NULL,NULL,NULL,0,&bins);

.. _capi-seams_brine_rejection:

seams_brine_rejection
---------------------


.. code-block:: c

   seams_status seams_brine_rejection(const seams_brine *front, double *c_ice, double *c_liquid, double *k, int *defined)

Molal rejection coefficient k = c_ice / c_liquid.

Molality is formula units per kilogram of water. front Brine handle. c_ice Bulk ice molality, or NULL. c_liquid Bulk liquid molality, or NULL. k The ratio, or NULL. Unchanged when the ratio is undefined. defined Receives 1 when the liquid bulk contains salt or the ice bulk does not, otherwise 0, or NULL. SEAMS_OK, or SEAMS_ERR_NULL when front is NULL or every output pointer is NULL. Thread safety Exclusive use of front. Example doublek=0; intdefined=0; seams_brine_rejection(front,NULL,NULL,&k,&defined);

.. _capi-seams_contact_pairs:

seams_contact_pairs
-------------------


.. code-block:: c

   seams_status seams_contact_pairs(seams_system *sys, const char *site_spec, size_t *n_pairs, int *n_cation, int *n_anion)

Mutual nearest unlike contacts of a site table.

Cation and anion counts are molecules in the ion cloud. sys Frame handle. site_spec Text such as "1=cationHead,2=anion". Must be non-NULL. n_pairs Receives the contact count, or NULL. n_cation Receives the cation molecule count, or NULL. n_anion Receives the anion molecule count, or NULL. At least one of the three outputs must be non-NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Thread safety Exclusive use of sys. Example size_tpairs=0; seams_contact_pairs(sys,"1=cationHead,2=anion",&pairs,NULL,NULL);

.. _capi-seams_coordination:

seams_coordination
------------------


.. code-block:: c

   seams_status seams_coordination(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *cn)

Coordination number integrated to rmax.

sys Frame handle. rmax Positive finite outer radius, Angstrom. bins Bin count. Zero selects about rmax / 0.1. type_i Central species, as in seams_rdf(). type_j Counted species, as in seams_rdf(). cn Receives the coordination number. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Thread safety Exclusive use of sys. Example doublecn=0; seams_coordination(sys,3.5,0,1,1,&cn);

.. _capi-seams_density_axis:

seams_density_axis
------------------


.. code-block:: c

   seams_status seams_density_axis(seams_system *sys, int type, int bins, int axis, double *z, double *rho, size_t cap, size_t *n_out)

Number density along one cell axis.

z is Angstrom. rho is number density, atoms per cubic Angstrom. sys Frame handle. type Species filter. Zero keeps every atom. bins Bin count. Zero selects about one bin per 0.1 Angstrom. axis 0, 1, or 2. z Bin centres, or NULL to query. rho Density per bin, or NULL to query. cap Capacity of each non-NULL column. n_out Receives the bin count. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns the columns. Thread safety Exclusive use of sys. Example doublez[50],rho[50]; size_tbins=0; seams_density_axis(sys,0,50,2,z,rho,50,&bins);

.. _capi-seams_domains:

seams_domains
-------------


.. code-block:: c

   seams_status seams_domains(seams_system *sys, double cutoff, const char *site_spec, int subset, int *n_subset, int *largest, double *percolation)

Largest domain of a polar or apolar site subset.

sys Frame handle. cutoff Positive finite cluster cutoff, Angstrom. site_spec Site table text. Must be non-NULL. subset 0 is apolar, 1 is polar. n_subset Receives the number of sites in the subset, or NULL. largest Receives the size of the largest domain, or NULL. percolation Receives the largest domain size over the subset size, or NULL. At least one output must be non-NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Thread safety Exclusive use of sys. Example intlargest=0; seams_domains(sys,5.0,"1=cationHead,2=anion,3=tail",0,NULL,&largest,NULL);

.. _capi-seams_fingerprint:

seams_fingerprint
-----------------


.. code-block:: c

   seams_status seams_fingerprint(seams_system *sys, double cutoff, int type, int k, int graph, int hops, int colour_types, seams_topo **out)

Topological fingerprint of the bonded graph.

sys Frame handle. cutoff Positive finite cutoff, Angstrom. type Species filter. k Neighbour count for the mutual graph. Zero selects 4. graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL. hops Neighbourhood depth. Zero selects 2. colour_types Non-zero colours vertices by LAMMPS type. out Receives a seams_topo. On failure \*out is NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller frees \*out with seams_fingerprint_destroy(). Thread safety Exclusive use of sys. Example seams_topo\*fp=NULL; seams_fingerprint(sys,3.5,0,4,SEAMS_GRAPH_CUTOFF,2,0,&fp); seams_fingerprint_destroy(fp);

.. _capi-seams_fingerprint_atom:

seams_fingerprint_atom
----------------------


.. code-block:: c

   seams_status seams_fingerprint_atom(const seams_topo *fp, int index, char *buf, size_t cap, size_t *n_out)

Local key of one atom. Same buffer rule as seams_fingerprint_key().

fp Fingerprint handle. index Atom index in the frame. buf Caller buffer of cap bytes, or NULL to query. cap Capacity of buf, including the terminating NUL. n_out Receives the key length without the terminating NUL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE for a bad index or a short buffer. Ownership The caller owns buf. Thread safety Exclusive use of fp. Example charkey[256]; size_tlen=0; seams_fingerprint_atom(fp,0,key,sizeofkey,&len);

.. _capi-seams_fingerprint_destroy:

seams_fingerprint_destroy
-------------------------


.. code-block:: c

   void seams_fingerprint_destroy(seams_topo *fp)

Free a fingerprint.

fp Handle from seams_fingerprint(), or NULL. NULL is ignored. Ownership Releases the handle and its key strings. Thread safety The caller must not use fp after this returns. Example seams_fingerprint_destroy(fp);

.. _capi-seams_fingerprint_key:

seams_fingerprint_key
---------------------


.. code-block:: c

   seams_status seams_fingerprint_key(const seams_topo *fp, char *buf, size_t cap, size_t *n_out)

Frame key text.

fp Fingerprint handle. buf Caller buffer of cap bytes, or NULL to query. cap Capacity of buf, including the terminating NUL. n_out Receives the key length without the terminating NUL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when buf is short. Ownership The caller owns buf. The copy is NUL-terminated. Thread safety Exclusive use of fp. Example size_tlen=0; seams_fingerprint_key(fp,NULL,0,&len); char\*key=malloc(len+1); seams_fingerprint_key(fp,key,len+1,&len);

.. _capi-seams_fingerprint_rings:

seams_fingerprint_rings
-----------------------


.. code-block:: c

   seams_status seams_fingerprint_rings(const seams_topo *fp, int *counts, size_t counts_len)

Primitive-ring census stored on the fingerprint.

fp Fingerprint handle. counts Caller buffer. counts[k] receives the number of rings of size k. Sizes past the buffer are dropped. counts_len Length of counts. Must be positive. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when counts_len is 0. Ownership The caller owns counts. Thread safety Exclusive use of fp. Example intcensus[13]; seams_fingerprint_rings(fp,census,13);

.. _capi-seams_hbonds:

seams_hbonds
------------


.. code-block:: c

   seams_status seams_hbonds(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, double dist_cutoff, double angle_cutoff, int *n_bonds)

Hydrogen-bond count on one frame.

sys Frame handle. cutoff Positive finite donor-acceptor cutoff, Angstrom. oxygen_type Heavy-atom species. Zero keeps every atom that is not hydrogen_type. hydrogen_type Positive hydrogen species. dist_cutoff Acceptor-H distance, Angstrom. angle_cutoff O-O-H angle, degrees. n_bonds Receives the number of hydrogen bonds. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Thread safety Exclusive use of sys. Example intn_bonds=0; seams_hbonds(sys,3.5,1,2,2.4,30.0,&n_bonds);

.. _capi-seams_rdf:

seams_rdf
---------


.. code-block:: c

   seams_status seams_rdf(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *g, int *count, size_t cap, size_t *n_out, double *volume)

Partial radial distribution.

r is Angstrom. g is dimensionless. volume is cubic Angstrom. sys Frame handle. rmax Positive finite outer radius, Angstrom. bins Bin count. Zero selects about rmax / 0.1. type_i First species. Zero and negative keep every atom of a single-species frame, or every atom as one species otherwise. type_j Second species. Zero and negative reuse the type_i species. r Bin centres, or NULL to query. g g(r) per bin, or NULL to query. count Pair count per bin, or NULL to query. cap Capacity of each non-NULL column. n_out Receives the bin count. volume Receives the cell volume, or NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns the columns. Thread safety Exclusive use of sys. Example size_tbins=0; seams_rdf(sys,8.0,80,0,0,NULL,NULL,NULL,0,&bins,NULL);

.. _capi-seams_running_cn:

seams_running_cn
----------------


.. code-block:: c

   seams_status seams_running_cn(seams_system *sys, double rmax, int bins, int type_i, int type_j, double *r, double *cn, size_t cap, size_t *n_out)

Running coordination number along the partial RDF.

sys Frame handle. rmax Positive finite outer radius, Angstrom. bins Bin count. Zero selects about rmax / 0.1. type_i Central species, as in seams_rdf(). type_j Counted species, as in seams_rdf(). r Bin centres, Angstrom, or NULL to query. cn Running coordination number, or NULL to query. cap Capacity of each non-NULL column. n_out Receives the bin count. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns the columns. Thread safety Exclusive use of sys. Example doubler[80],cn[80]; size_tbins=0; seams_running_cn(sys,8.0,80,0,0,r,cn,80,&bins);

.. _capi-seams_site_table:

seams_site_table
----------------


.. code-block:: c

   seams_status seams_site_table(seams_system *sys, const char *spec, int *kinds, size_t len)

Per-atom site kind from a site-table spec.

sys Frame handle. spec Site table text such as "1=cationHead,2=anion". Must be non-NULL. kinds One seams_site_kind per atom. len Length of kinds. Must be at least the frame count. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns kinds. Thread safety Exclusive use of sys. Example seams_site_table(sys,"1=cationHead,2=anion",kinds,n);

