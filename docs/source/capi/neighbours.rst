Neighbours
==========


Build a neighbour list and read it back.

The stable C interface provides 129 exported functions through ``include/seams.h``.

.. literalinclude:: ../../../docs/examples/capi_neighbours.c
   :language: c

.. _capi-seams_neighbours_build:

seams_neighbours_build
----------------------


.. code-block:: c

   seams_status seams_neighbours_build(seams_system *sys, double cutoff, int type, int graph, seams_neighbours **out)

Build a cutoff or mutual neighbour list.

``sys``
    Frame handle that already holds a frame.

``cutoff``
    Positive finite cutoff, Angstrom. The mutual graph keeps the four nearest neighbours inside cutoff + 1.5.

``type``
    Species to keep. Zero and negative keep every atom.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL. ::SEAMS_GRAPH_SEEDED is rejected here.

``out``
    Receives a new handle. On failure \*out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees \*out with seams_neighbours_destroy(). The list is a copy and does not borrow sys. @par Thread safety Exclusive use of sys. The new handle is private to the caller. @par Example @code seams_neighbours \*nb = NULL; if (seams_neighbours_build(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, &nb) == SEAMS_OK) { seams_neighbours_destroy(nb); } @endcode

.. _capi-seams_neighbours_destroy:

seams_neighbours_destroy
------------------------


.. code-block:: c

   void seams_neighbours_destroy(seams_neighbours *nb)

Free a neighbour list.

``nb``
    Handle from seams_neighbours_build() or seams_neighbours_knn(), or NULL. NULL is ignored.

@par Ownership Releases the handle and its CSR arrays. @par Thread safety The caller must not use nb after this returns. @par Example @code seams_neighbours_destroy(nb); @endcode

.. _capi-seams_neighbours_count:

seams_neighbours_count
----------------------


.. code-block:: c

   int seams_neighbours_count(const seams_neighbours *nb)

Number of neighbour rows, as an int.

``nb``
    Neighbour handle, or NULL.

Returns: The row count, or 0 when nb is NULL.

ABI 1 query. seams_neighbours_nrows() reports errors.

@par Ownership The caller owns nb; the call borrows it. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on nb, and not alongside a non-const call on nb. @par Example @code int rows = seams_neighbours_count(nb); @endcode

.. _capi-seams_neighbours_csr:

seams_neighbours_csr
--------------------


.. code-block:: c

   seams_status seams_neighbours_csr(const seams_neighbours *nb, const int **offsets, int *n_offsets, const int **cols, int *nnz)

Deprecated. use seams_neighbours_copy; these pointers are owned by the handle.

Borrow the neighbour list as CSR.

``nb``
    Neighbour handle.

``offsets``
    Receives the offset pointer.

``n_offsets``
    Receives the offset count, n + 1.

``cols``
    Receives the column pointer. May be NULL when the list has no pairs; \*nnz is then 0.

``nnz``
    Receives the column count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL.

Offsets has n + 1 entries. Columns omit the self entry.

@par Ownership The pointers are owned by nb and remain valid until seams_neighbours_destroy(nb). Use seams_neighbours_copy() for caller buffers. @par Thread safety Exclusive use of nb. The returned pointers must not be used concurrently with a call that frees nb. @par Example @code const int \*off = NULL, \*cols = NULL; int n_off = 0, nnz = 0; seams_neighbours_csr(nb, &off, &n_off, &cols, &nnz); @endcode

.. _capi-seams_neighbours_copy:

seams_neighbours_copy
---------------------


.. code-block:: c

   seams_status seams_neighbours_copy(const seams_neighbours *nb, int *offsets, size_t offsets_cap, size_t *n_offsets, int *cols, size_t cols_cap, size_t *nnz)

Copy the neighbour CSR into caller buffers.

``nb``
    Neighbour handle.

``offsets``
    Caller buffer of offsets_cap ints, or NULL to query.

``offsets_cap``
    Capacity of offsets. Ignored when offsets is NULL.

``n_offsets``
    Receives the required offset count (n + 1).

``cols``
    Caller buffer of cols_cap ints, or NULL to query.

``cols_cap``
    Capacity of cols. Ignored when cols is NULL.

``nnz``
    Receives the required column count.

Returns: ::SEAMS_OK when every non-NULL buffer was filled, or when both buffers were NULL and the lengths were written. ::SEAMS_ERR_RANGE when a non-NULL buffer is short. ::SEAMS_ERR_NULL when nb, n_offsets, or nnz is NULL.

@par Ownership The caller owns the buffers. The handle is not modified. @par Thread safety Read-only: may run alongside other const calls on nb, and not alongside a non-const call on nb. @par Example @code size_t n_off = 0, nnz = 0; seams_neighbours_copy(nb, NULL, 0, &n_off, NULL, 0, &nnz); int \*off = malloc(n_off \* sizeof \*off), \*cols = malloc(nnz \* sizeof \*cols); seams_neighbours_copy(nb, off, n_off, &n_off, cols, nnz, &nnz); free(off); free(cols); @endcode

.. _capi-seams_neighbours_knn:

seams_neighbours_knn
--------------------


.. code-block:: c

   seams_status seams_neighbours_knn(seams_system *sys, int k, double cutoff, int type, int mutual, seams_neighbours **out)

k-nearest neighbour list.

``sys``
    Frame handle.

``k``
    Positive neighbour count.

``cutoff``
    Nominal positive cutoff, Angstrom. Candidates are taken inside cutoff + 1.5.

``type``
    Species filter. Zero and negative keep every atom.

``mutual``
    Non-zero keeps the mutual graph. Zero keeps the one-way k-nearest graph.

``out``
    Receives a neighbour handle. On failure \*out is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry, so seams_handle_error() reports this call. It reads the frame and keeps no scratch on the handle. @par Ownership The caller frees \*out with seams_neighbours_destroy(). @par Thread safety Exclusive use of sys. @par Example @code seams_neighbours \*nb = NULL; seams_neighbours_knn(sys, 4, 3.5, 0, 1, &nb); seams_neighbours_destroy(nb); @endcode

.. _capi-seams_neighbours_nrows:

seams_neighbours_nrows
----------------------


.. code-block:: c

   seams_status seams_neighbours_nrows(const seams_neighbours *nb, size_t *n)

Number of neighbour rows.

``nb``
    Neighbour handle.

``n``
    Receives the row count.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL when nb or n is NULL.

@par Ownership The caller owns nb; the call borrows it. n is caller storage that the call fills. The library keeps no pointer to any argument. @par Thread safety Read-only: may run alongside other const calls on nb, and not alongside a non-const call on nb. @par Example @code size_t rows = 0; seams_neighbours_nrows(nb, &rows); @endcode

