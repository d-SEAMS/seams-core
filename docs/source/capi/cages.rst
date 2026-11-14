Cages
=====


Hexagonal, cubic, and signature cages.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 69 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

.. literalinclude:: ../../../docs/examples/capi_cages.c
   :language: c

.. _capi-seams_cages_find:

seams_cages_find
----------------


.. code-block:: c

   seams_status seams_cages_find(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)

Hexagonal and double-diamond cage membership.

``sys``
    Frame handle.

``cutoff``
    Positive cutoff, Angstrom. The seeded graph uses cutoff + 1.5 as the candidate shell.

``type``
    Species filter.

``graph``
    SEAMS_GRAPH_CUTOFF, SEAMS_GRAPH_MUTUAL, or SEAMS_GRAPH_SEEDED.

``k``
    Neighbour count for the seeded graph. Zero selects 4.

``complete``
    Non-zero runs ring-adjacent completion on the seeded labels. Ignored for the cutoff and mutual graphs.

``out``
    Receives a new handle. On failure \*out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: the caller frees \*out with seams_cages_destroy(). Thread safety: exclusive use of sys.

.. _capi-seams_cages_destroy:

seams_cages_destroy
-------------------


.. code-block:: c

   void seams_cages_destroy(seams_cage_set *cages)

Free a cage-membership handle.

``cages``
    Handle from seams_cages_find(), or NULL.

Thread safety: the caller must not use cages after this returns.

.. _capi-seams_cages_counts:

seams_cages_counts
------------------


.. code-block:: c

   seams_status seams_cages_counts(const seams_cage_set *cages, int *n_hex, int *n_cubic, int *n_water, int *n_both)

Disjoint cage population counts.

``cages``
    Cage handle.

``n_hex``
    Hexagonal only, or NULL.

``n_cubic``
    Double-diamond only, or NULL.

``n_water``
    Neither cage, or NULL.

``n_both``
    Both cages, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL when cages is NULL. At least one count pointer must be non-NULL.

The four counts sum to the number of atoms the search kept. Thread safety: exclusive use of cages.

.. _capi-seams_cages_labels:

seams_cages_labels
------------------


.. code-block:: c

   seams_status seams_cages_labels(const seams_cage_set *cages, int *per_atom, size_t len)

Per-atom cage code.

``cages``
    Cage handle.

``per_atom``
    Caller buffer of seams_cage_label values, one per atom in the frame. Atoms the species filter dropped are SEAMS_CAGE_WATER.

``len``
    Length of per_atom. Must be at least the frame count stored on the handle.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: the caller owns per_atom. Thread safety: exclusive use of cages.

.. _capi-seams_cages_signature:

seams_cages_signature
---------------------


.. code-block:: c

   seams_status seams_cages_signature(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)

Count closed cages whose ring census matches a signature.

``sys``
    Frame handle.

``cutoff``
    Positive cutoff, Angstrom.

``type``
    Species filter.

``graph``
    SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.

``spec``
    Signature text. Must be non-NULL.

``n_cages``
    Receives the cage count.

``n_atoms``
    Receives the vertex count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

spec is a census such as "512", "51262", or "51264", or the names "hc" and "ddc". n_cages is the number of closed cages. n_atoms is the number of distinct vertices in those cages.

Thread safety: exclusive use of sys.

.. _capi-seams_system_symmetry:

seams_system_symmetry
---------------------


.. code-block:: c

   seams_status seams_system_symmetry(seams_system *sys, double shell_rmsd, int *labels, size_t labels_len)

Per-atom ice shell label.

``sys``
    Frame handle.

``shell_rmsd``
    Positive scale-free residual ceiling, or a non-positive value for the library ceiling.

``labels``
    Caller buffer.

``labels_len``
    Must be at least the frame count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_INTERNAL, SEAMS_ERR_NOMEM.

labels[i] is 0 for ice Ih, 1 for ice Ic, 2 for stacking disorder, or -1 when the shell does not match. A non-positive shell_rmsd selects the library ceiling (the scale-free second-shell residual, 0.15).

Ownership: the caller owns labels. Thread safety: exclusive use of sys.

.. _capi-seams_cages:

seams_cages
-----------


.. code-block:: c

   seams_status seams_cages(seams_system *sys, double cutoff, int type, int graph, int k, int complete, seams_cage_set **out)

Hexagonal and double-diamond cages. Same search as seams_cages_find().

``out``
    Receives a seams_cage_set. The caller frees it with seams_cages_destroy().

Thread safety: exclusive use of sys.

.. _capi-seams_signature_cages:

seams_signature_cages
---------------------


.. code-block:: c

   seams_status seams_signature_cages(seams_system *sys, double cutoff, int type, int graph, const char *spec, size_t *n_cages, size_t *n_atoms)

Closed cages of one ring census. Same search as seams_cages_signature().

Thread safety: exclusive use of sys.

.. _capi-seams_prisms:

seams_prisms
------------


.. code-block:: c

   seams_status seams_prisms(seams_system *sys, double cutoff, int type, int *perfect, int *imperfect, size_t len)

Perfect and imperfect prism counts by ring size.

``perfect,``
    imperfect Length at least 7. Index is the ring size. Sizes 3 through 6 are filled.

Thread safety: exclusive use of sys.

.. _capi-seams_monolayer:

seams_monolayer
---------------


.. code-block:: c

   seams_status seams_monolayer(seams_system *sys, double cutoff, int type, int max_depth, double sheet_area, int *counts, size_t counts_len)

Primitive-ring histogram of a monolayer.

``sheet_area``
    Sheet area, square Angstrom. Must be positive.

``counts``
    counts[k] is the number of rings of size k. counts_len must be at least max_depth + 1.

Thread safety: exclusive use of sys.

