Rings
=====


Primitive rings and their histogram.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 85 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

.. literalinclude:: ../../../docs/examples/capi_rings.c
   :language: c

.. _capi-seams_ring_symmetry_name:

seams_ring_symmetry_name
------------------------


.. code-block:: c

   const char *seams_ring_symmetry_name(int symmetry)

Name of a ring-symmetry code.

``symmetry``
    A ::seams_ring_symmetry value.

Returns: Static string. Unknown codes return "unknown".

@par Ownership The caller does not free it. @par Thread safety Pure. @par Example @code printf("%s\\n", seams_ring_symmetry_name(SEAMS_RING_DDC)); @endcode

.. _capi-seams_rings_find:

seams_rings_find
----------------


.. code-block:: c

   seams_status seams_rings_find(seams_system *sys, double cutoff, int type, int graph, int max_depth, seams_rings **out)

Primitive rings (Franzblau) up to a maximum size.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

``max_depth``
    Largest ring size, from 3 through 12.

``out``
    Receives a new handle. On failure \*out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees \*out with seams_rings_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_rings \*rings = NULL; if (seams_rings_find(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, 6, &rings) == SEAMS_OK) { seams_rings_destroy(rings); } @endcode

.. _capi-seams_rings_destroy:

seams_rings_destroy
-------------------


.. code-block:: c

   void seams_rings_destroy(seams_rings *rings)

Free a ring set.

``rings``
    Handle from seams_rings_find(), or NULL. NULL is ignored.

@par Ownership Releases the handle and every ring it holds. @par Thread safety The caller must not use rings after this returns. @par Example @code seams_rings_destroy(rings); @endcode

.. _capi-seams_rings_count:

seams_rings_count
-----------------


.. code-block:: c

   int seams_rings_count(const seams_rings *rings)

Number of primitive rings, as an int.

``rings``
    Ring handle, or NULL.

Returns: The ring count, or 0 when rings is NULL.

ABI 1 query. seams_rings_nrings() reports errors.

@par Ownership The caller owns rings; the call borrows it. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int nr = seams_rings_count(rings); @endcode

.. _capi-seams_rings_histogram:

seams_rings_histogram
---------------------


.. code-block:: c

   seams_status seams_rings_histogram(const seams_rings *rings, int *counts, int counts_len)

Count primitive rings by size, with an int length.

``rings``
    Ring handle.

``counts``
    Caller buffer. counts[k] receives the number of rings of size k. Sizes outside the buffer are ignored.

``counts_len``
    Length of counts. Must be positive.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

ABI 1 entry. New callers use seams_rings_histogram_n().

@par Ownership The caller owns counts. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int hist[13]; seams_rings_histogram(rings, hist, 13); @endcode

.. _capi-seams_rings_labels:

seams_rings_labels
------------------


.. code-block:: c

   seams_status seams_rings_labels(const seams_rings *rings, int *sizes, int *symmetry)

Size and symmetry of each ring, without buffer lengths.

``rings``
    Ring handle.

``sizes``
    One int per ring, or NULL.

``symmetry``
    One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

ABI 1 entry. Each non-NULL buffer has length seams_rings_count(rings). New callers use seams_rings_labels_n(), which checks the lengths.

@par Ownership The caller owns both arrays. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. The first symmetry request fills a cache on the handle under a lock; later requests only read it. @par Example @code seams_rings_labels(rings, sizes, symmetry); @endcode

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

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE.

@par Ownership The caller owns atoms. Indices are positions in the frame. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int atoms[12]; size_t len = 0; seams_rings_vertices(rings, 0, atoms, 12, &len); @endcode

.. _capi-seams_rings_chordless:

seams_rings_chordless
---------------------


.. code-block:: c

   seams_status seams_rings_chordless(const seams_neighbours *nb, int size, int *rings, size_t rings_len, size_t *n_rings, int *per_atom, size_t per_atom_len)

Chordless cycles of one size.

``nb``
    Graph from seams_neighbours_build().

``size``
    Cycle length, from 3 through 12.

``rings``
    Caller buffer of size times the cycle count ints, or NULL to skip the list.

``rings_len``
    Length of rings in ints. Ignored when rings is NULL.

``n_rings``
    Receives the cycle count. Required when rings is not NULL.

``per_atom``
    Caller buffer of one int per atom, or NULL. Receives the number of cycles through each atom.

``per_atom_len``
    Length of per_atom. Ignored when per_atom is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL. rings, n_rings and per_atom may not all be NULL. A short rings buffer returns ::SEAMS_ERR_RANGE after n_rings and per_atom are written.

A cycle is chordless when no edge of the graph joins two of its vertices that are not next to each other on the cycle. Every primitive ring is chordless; a chordless cycle need not be primitive. Each cycle is written once as size atom indices: it starts at its lowest index and runs towards the smaller of that atom's two cycle neighbours. Cycles are sorted lexicographically, so the list does not depend on the thread count. @par Ownership The caller owns rings and per_atom. Indices are positions in the frame. @par Thread safety Exclusive use of nb. @par Example @code size_t n_six = 0; seams_rings_chordless(nb, 6, NULL, 0, &n_six, NULL, 0); @endcode

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

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_NOMEM.

@par Ownership The caller owns rings; the call borrows it. fraction is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code double f6 = 0; seams_rings_fraction(rings, 6, &f6); @endcode

.. _capi-seams_rings_nrings:

seams_rings_nrings
------------------


.. code-block:: c

   seams_status seams_rings_nrings(const seams_rings *rings, size_t *n)

Number of primitive rings.

``rings``
    Ring handle.

``n``
    Receives the ring count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL when rings or n is NULL.

@par Ownership The caller owns rings; the call borrows it. n is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code size_t nr = 0; seams_rings_nrings(rings, &nr); @endcode

.. _capi-seams_rings_labels_n:

seams_rings_labels_n
--------------------


.. code-block:: c

   seams_status seams_rings_labels_n(const seams_rings *rings, int *sizes, size_t sizes_len, int *symmetry, size_t symmetry_len)

Size and symmetry of each ring into buffers of known length.

``rings``
    Ring handle.

``sizes``
    One int per ring, or NULL.

``sizes_len``
    Length of sizes. Ignored when sizes is NULL.

``symmetry``
    One ::seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL.

``symmetry_len``
    Length of symmetry. Ignored when symmetry is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when a buffer is shorter than the ring count, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Ownership The caller owns both arrays. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code size_t nr = 0; seams_rings_nrings(rings, &nr); int \*sym = malloc(nr \* sizeof \*sym); seams_rings_labels_n(rings, NULL, 0, sym, nr); free(sym); @endcode

.. _capi-seams_rings_histogram_n:

seams_rings_histogram_n
-----------------------


.. code-block:: c

   seams_status seams_rings_histogram_n(const seams_rings *rings, int *counts, size_t counts_len)

Count primitive rings by size into a buffer of known length.

``rings``
    Ring handle.

``counts``
    Caller buffer. counts[k] receives the number of rings of size k. Sizes outside the buffer are ignored.

``counts_len``
    Length of counts. Must be positive.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_RANGE when counts_len is 0.

@par Ownership The caller owns counts. @par Thread safety Read-only: may run alongside other const calls on rings, and not alongside a non-const call on rings. @par Example @code int hist[13]; seams_rings_histogram_n(rings, hist, 13); @endcode

