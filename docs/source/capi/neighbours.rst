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

Type 0 on the atoms (1, 1, 1) and (2, 1.2, 1.1) in an orthorhombic cell of edge 1e160 with origin (0, 0, 0) returns ::SEAMS_OK from the cutoff graph: an empty list at cutoff 0.5 and each other at cutoff 3.5, the same rows as that pair and type on an edge of 30 with that origin. The same atoms and type with origin (5, -2, 0.25) list each other at cutoff 0.5 on the edge of 1e160 and do not on an edge of 30. Type 0 on (-0.9, 0, 0) with (0.1, 0.2, 0.1), and on (-2, -1, -4) with (-1, -0.8, -3.9), returns ::SEAMS_OK and both neighbours at cutoff 0.5 on edge 1e160 with origin (0, 0, 0), and an empty list on edge 30. Type 0 on (0, 0, 0) and 0x1.6c2d4256ffcc4p+588 returns ::SEAMS_OK and both neighbours at cutoff 0.5 on that edge and origin.

When the orthorhombic cutoff graph bins a coordinate, the fraction is (coordinate - origin) times the binary64 reciprocal of the edge. When that fraction is outside [0, 1) it subtracts the floor. A result that is not strictly below 1, NaN included, becomes 0. The stored coordinate is the origin plus the fraction times the edge. Type 0 on (0, 0, 0) and (1e162, 0, 0), origin (0, 0, 0), edge 1e160, cutoff 0.5, returns ::SEAMS_OK and an empty list. The product is 99.99999999999999. The floor leaves 0.9999999999999858, and the stored x coordinate is 9.999999999999858e159. Type 0 on (0, 0, 0) and (1.5e160, 0, 0), that edge and origin, cutoff 0.5, returns ::SEAMS_OK and an empty list. The product is 1.5. The floor leaves 0.5, and the stored x coordinate is 5e159. Type 0 on (3.2581295480326344e+46, 0, 0) and (3.2581295480326338e+46, 0, 0), that edge and origin, cutoff 0.5, returns ::SEAMS_OK and both neighbours. Both stored coordinates are 3.2581295480326338e+46. Type 0 on (5, 0, 0) and (4096, 0, 0), origin (5, 0, 0), edge 1e160, cutoff 4091, returns ::SEAMS_OK and both neighbours. The wrap stores x = 4096 at 4095.9999999999995. Type 0 on (-1e308, 0, 0) and (1e308, 0, 0), origin (-1e308, 0, 0), edge 30, cutoff 0.5, returns ::SEAMS_OK and both neighbours. The product for the second atom is inf, inf minus its floor is NaN, and that NaN becomes 0, so both atoms are stored at the origin. Type 0 on (0, 0, 0) and (1, 0, 0), origin (0, 0, 0), edge 2^-1070, cutoff 0.5, returns ::SEAMS_OK and both neighbours. The reciprocal of that edge is inf, both fractions become 0, and both atoms are stored at the origin.

Type 0 on the mutual graph returns ::SEAMS_ERR_RANGE for the atoms (1, 1, 1) and (2, 1.2, 1.1) on the orthorhombic cell of edge 1e160 with origin (0, 0, 0), and on the triclinic cell of that edge with that origin and tilts 1, 0.5 and -0.25, at cutoffs 0.5 and 3.5. Type 99 on each of those frames returns ::SEAMS_OK and an empty list at those cutoffs. On that triclinic cell, type 0 on the cutoff graph returns ::SEAMS_OK and both neighbours for those atoms at cutoff 3.5, and ::SEAMS_ERR_RANGE at cutoff 0.5. Type 99 at that cutoff returns ::SEAMS_OK and an empty list.

Type 0 on (0, 0, 0) and (10 \* 2^52, 0, 0), edge 10, origin (0, 0, 0) and tilt xy = 1, returns ::SEAMS_ERR_RANGE from the cutoff graph at cutoffs 0.5 and 3.5. Type 0 on (5, 5, 5) and (45035996273704968, 5, 5) on that cell returns ::SEAMS_OK and an empty list at cutoff 0.5, and ::SEAMS_ERR_RANGE at cutoff 3.5. Type 0 on (0, 0, 0) and (2^53 - 1, 0, 0), edges 2, 2 and 2, origin (0, 0, 0) and tilt xy = 0.25, returns ::SEAMS_ERR_RANGE at cutoffs 0.5 and 1.0, with the message that the minimum-image wrap is not a lattice image, and ::SEAMS_OK with both neighbours at cutoff 3.5. Type 0 on (0, 0, 0) and (50 \* 2^52, 10 \* 2^51, 0), edges 10, 10 and 10, origin (0, 0, 0) and tilt xy = 100, returns ::SEAMS_OK and both neighbours at cutoffs 0.5 and 3.5. Type 0 on (0, 0, 0) and (10 \* 2^52, 0, 0) in an orthorhombic cell of edge 10 with origin (0, 0, 0) returns ::SEAMS_OK and both neighbours at cutoffs 0.5 and 3.5. Type 0 on (4, 0, 0) and (4 + 10 \* 2^52, 0, 0) in that cell returns ::SEAMS_OK and an empty list at both cutoffs. Type 0 on (0, 0, 0) and (0.5, 0, 0) in that cell returns ::SEAMS_OK and an empty list at cutoff 0.5.

For a pair the orthorhombic cutoff cell list visits, it keeps the neighbour only when the squared separation of the wrapped coordinates and a periodic cell shift is strictly below the squared cutoff, both computed in binary64. An overflowing square and a NaN square fail that test. Type 0 on (0, 0, 0) and (2e154, 0, 0), edge 1e160, origin (0, 0, 0), cutoff 3e154, returns ::SEAMS_OK and an empty list. Type 0 on (0.9e308, 0, 0) and (0.901e308, 0, 0), origin (1e308, 0, 0), edge 1e308, cutoff 0.5, returns ::SEAMS_OK and an empty list. Both stored x coordinates are inf, so each squared separation is NaN.

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

