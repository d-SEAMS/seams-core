Rings
=====


Primitive rings and their histogram.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 69 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

.. literalinclude:: ../../../docs/examples/capi_rings.c
   :language: c

.. _capi-seams_ring_symmetry_name:

seams_ring_symmetry_name
------------------------


.. code-block:: c

   const char *seams_ring_symmetry_name(int symmetry)

Name of a ring-symmetry code.

``symmetry``
    A seams_ring_symmetry value.

Returns: Static string. Unknown codes return "unknown".

Thread safety: pure.

.. _capi-seams_rings_find:

seams_rings_find
----------------


.. code-block:: c

   seams_status seams_rings_find(seams_system *sys, double cutoff, int type, int graph, int max_depth, seams_rings **out)

Primitive rings (Franzblau) up to a maximum size.

``sys``
    Frame handle.

``cutoff``
    Positive cutoff, Angstrom.

``type``
    Species filter.

``graph``
    SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.

``max_depth``
    Largest ring size, from 3 through 12.

``out``
    Receives a new handle. On failure \*out is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Ownership: the caller frees \*out with seams_rings_destroy(). Thread safety: exclusive use of sys.

.. _capi-seams_rings_destroy:

seams_rings_destroy
-------------------


.. code-block:: c

   void seams_rings_destroy(seams_rings *rings)

Free a ring set.

``rings``
    Handle from seams_rings_find(), or NULL.

Thread safety: the caller must not use rings after this returns.

.. _capi-seams_rings_count:

seams_rings_count
-----------------


.. code-block:: c

   int seams_rings_count(const seams_rings *rings)

Number of primitive rings.

``rings``
    Ring handle, or NULL.

Returns: The ring count, or 0 when rings is NULL.

Thread safety: exclusive use of rings.

.. _capi-seams_rings_histogram:

seams_rings_histogram
---------------------


.. code-block:: c

   seams_status seams_rings_histogram(const seams_rings *rings, int *counts, int counts_len)

Count primitive rings by size.

``rings``
    Ring handle.

``counts``
    Caller buffer. counts[k] receives the number of rings of size k. Sizes outside the buffer are ignored.

``counts_len``
    Length of counts. Must be positive.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: the caller owns counts. Thread safety: exclusive use of rings.

.. _capi-seams_rings_labels:

seams_rings_labels
------------------


.. code-block:: c

   seams_status seams_rings_labels(const seams_rings *rings, int *sizes, int *symmetry)

Size and symmetry of each ring.

``rings``
    Ring handle.

``sizes``
    One int per ring, or NULL.

``symmetry``
    One seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL. Each non-NULL buffer has length seams_rings_count(rings).

Returns: SEAMS_OK, or SEAMS_ERR_NULL.

Ownership: the caller owns both arrays. Thread safety: exclusive use of rings.

.. _capi-seams_rings_vertices:

seams_rings_vertices
--------------------


.. code-block:: c

   seams_status seams_rings_vertices(const seams_rings *rings, int index, int *atoms, size_t cap, size_t *n_out)

Atom indices of one ring.

``rings``
    Ring handle.

``index``
    Ring index, from 0 to seams_rings_count(rings) - 1.

``atoms``
    Caller buffer, or NULL to query.

``cap``
    Capacity of atoms. Ignored when atoms is NULL.

``n_out``
    Receives the vertex count.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE.

Ownership: the caller owns atoms. Indices are positions in the frame. Thread safety: exclusive use of rings.

.. _capi-seams_rings_fraction:

seams_rings_fraction
--------------------


.. code-block:: c

   seams_status seams_rings_fraction(const seams_rings *rings, int size, double *fraction)

Fraction of atoms that sit on a ring.

``rings``
    Ring handle.

``size``
    Ring size to count. Zero counts atoms that sit on any ring.

``fraction``
    Receives a value in [0, 1].

Returns: SEAMS_OK, or SEAMS_ERR_NULL.

Thread safety: exclusive use of rings.

