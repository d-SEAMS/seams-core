Cages
=====


Hexagonal and cubic cage membership.

Generated from the Doxygen XML of ``include/seams.h`` at ``capi/polish-2.14`` commit ``0127f41efc3f44cefdcf532bdea3bba36b05176d``. The C sources of that branch are not part of this docs tree.

``meson test docs_capi_examples`` compiles and runs this example.

.. literalinclude:: ../../../docs/examples/capi_cages.c
   :language: c

.. _capi-seams_cages_find:

seams_cages_find
----------------


.. code-block:: c

   seams_status seams_cages_find(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)

Hexagonal and double-diamond cage membership.

sys Frame handle. cutoff Positive finite cutoff, Angstrom. The seeded graph uses cutoff + 1.5 as the candidate shell. type Species filter. graph SEAMS_GRAPH_CUTOFF, SEAMS_GRAPH_MUTUAL, or SEAMS_GRAPH_SEEDED. k Neighbour count for the seeded graph. Zero selects 4. complete Non-zero runs ring-adjacent completion on the seeded labels. Ignored for the cutoff and mutual graphs. out Receives a new handle. On failure \*out is NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller frees \*out with seams_cages_destroy(). Thread safety Exclusive use of sys. Example seams_cage_set\*cages=NULL; if(seams_cages_find(sys,3.5,0,SEAMS_GRAPH_CUTOFF,0,0,&cages)==SEAMS_OK){ seams_cages_destroy(cages); }

.. _capi-seams_cages:

seams_cages
-----------


.. code-block:: c

   seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)

Hexagonal and double-diamond cages. Same search as seams_cages_find().

sys Frame handle. cutoff Positive finite cutoff, Angstrom. type Species filter. graph SEAMS_GRAPH_CUTOFF, SEAMS_GRAPH_MUTUAL, or SEAMS_GRAPH_SEEDED. k Neighbour count for the seeded graph. Zero selects 4. complete Non-zero runs ring-adjacent completion on seeded labels. out Receives a seams_cage_set. On failure \*out is NULL. Same codes as seams_cages_find(). Ownership The caller frees \*out with seams_cages_destroy(). Thread safety Exclusive use of sys. Example seams_cage_set\*cages=NULL; seams_cages(sys,3.5,0,SEAMS_GRAPH_SEEDED,4,1,&cages); seams_cages_destroy(cages);

.. _capi-seams_cages_destroy:

seams_cages_destroy
-------------------


.. code-block:: c

   void seams_cages_destroy(seams_cage_set *cages)

Free a cage-membership handle.

cages Handle from seams_cages_find() or seams_cages(), or NULL. NULL is ignored. Ownership Releases the handle and its per-atom labels. Thread safety The caller must not use cages after this returns. Example seams_cages_destroy(cages);

.. _capi-seams_cages_counts:

seams_cages_counts
------------------


.. code-block:: c

   seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex, int *n_cubic, int *n_water, int *n_both)

Disjoint cage population counts.

The four counts sum to the number of atoms the search kept. cages Cage handle. n_hex Hexagonal only, or NULL. n_cubic Double-diamond only, or NULL. n_water Neither cage, or NULL. n_both Both cages, or NULL. SEAMS_OK, or SEAMS_ERR_NULL when cages is NULL or every count pointer is NULL. Thread safety Exclusive use of cages. Example inthex=0,cubic=0; seams_cages_counts(cages,&hex,&cubic,NULL,NULL);

.. _capi-seams_cages_labels:

seams_cages_labels
------------------


.. code-block:: c

   seams_status seams_cages_labels(const seams_cage_set *cages, int *per_atom, size_t len)

Per-atom cage code.

cages Cage handle. per_atom Caller buffer of seams_cage_label values, one per atom in the frame. Atoms the species filter dropped are SEAMS_CAGE_WATER. len Length of per_atom. Must be at least the frame count stored on the handle. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE. Ownership The caller owns per_atom. Thread safety Exclusive use of cages. Example seams_cages_labels(cages,per_atom,n);

.. _capi-seams_cages_signature:

seams_cages_signature
---------------------


.. code-block:: c

   seams_status seams_cages_signature(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)

Count closed cages whose ring census matches a signature.

spec is a census such as "512", "51262", or "51264", or the names "hc" and "ddc". n_cages is the number of closed cages. n_atoms is the number of distinct vertices in those cages. sys Frame handle. cutoff Positive finite cutoff, Angstrom. type Species filter. graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL. spec Signature text. Must be non-NULL. n_cages Receives the cage count. n_atoms Receives the vertex count. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Thread safety Exclusive use of sys. Example size_tn_cages=0,n_atoms=0; seams_cages_signature(sys,3.5,0,SEAMS_GRAPH_CUTOFF,"512",&n_cages,&n_atoms);

.. _capi-seams_signature_cages:

seams_signature_cages
---------------------


.. code-block:: c

   seams_status seams_signature_cages(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)

Closed cages of one ring census. Same search as seams_cages_signature().

sys Frame handle. cutoff Positive finite cutoff, Angstrom. type Species filter. graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL. spec Signature text such as "512" or "hc". Must be non-NULL. n_cages Receives the cage count. n_atoms Receives the vertex count. Same codes as seams_cages_signature(). Thread safety Exclusive use of sys. Example size_tn_cages=0,n_atoms=0; seams_signature_cages(sys,3.5,0,SEAMS_GRAPH_CUTOFF,"hc",&n_cages,&n_atoms);

.. _capi-seams_system_symmetry:

seams_system_symmetry
---------------------


.. code-block:: c

   seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd, int *labels, size_t labels_len)

Per-atom ice shell label.

labels[i] is 0 for ice Ih, 1 for ice Ic, 2 for stacking disorder, or -1 when the shell does not match. A non-positive shell_rmsd selects the library ceiling (the scale-free second-shell residual, 0.15). sys Frame handle. shell_rmsd Positive scale-free residual ceiling, or a non-positive value for the library ceiling. labels Caller buffer. labels_len Must be at least the frame count. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_INTERNAL, SEAMS_ERR_NOMEM. Ownership The caller owns labels. Thread safety Exclusive use of sys. Example seams_system_symmetry(sys,0.0,labels,n);

.. _capi-seams_prisms:

seams_prisms
------------


.. code-block:: c

   seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect, int *imperfect, size_t len)

Perfect and imperfect prism counts by ring size.

Rings come from the mutual four-nearest graph inside cutoff + 1.5. sys Frame handle. cutoff Positive finite cutoff, Angstrom. type Species filter. perfect Caller buffer. perfect[k] receives the perfect k-prisms. imperfect Caller buffer. imperfect[k] receives the deformed k-prisms. len Length of each buffer. At least 7, so that sizes 3 through 6 are filled. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns both buffers. Thread safety Exclusive use of sys. Example intperfect[7],imperfect[7]; seams_prisms(sys,3.5,0,perfect,imperfect,7);

.. _capi-seams_monolayer:

seams_monolayer
---------------


.. code-block:: c

   seams_status seams_monolayer(seams_system *sys, double cutoff, int type, int max_depth, double sheet_area, int *counts, size_t counts_len)

Primitive-ring histogram of a monolayer.

sys Frame handle. cutoff Positive finite cutoff, Angstrom. type Species filter. max_depth Largest ring size, at least 3. sheet_area Sheet area, square Angstrom. Must be positive. counts counts[k] is the number of rings of size k. counts_len Must be at least max_depth + 1. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns counts. Thread safety Exclusive use of sys. Example intcounts[7]; seams_monolayer(sys,3.5,0,6,900.0,counts,7);

