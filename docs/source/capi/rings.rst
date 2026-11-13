Rings
=====


Primitive rings and their histogram.

Generated from the Doxygen XML of ``include/seams.h`` at ``capi/polish-2.14`` commit ``0127f41efc3f44cefdcf532bdea3bba36b05176d``. The C sources of that branch are not part of this docs tree.

``meson test docs_capi_examples`` compiles and runs this example.

.. literalinclude:: ../../../docs/examples/capi_rings.c
   :language: c

.. _capi-seams_rings_find:

seams_rings_find
----------------


.. code-block:: c

   seams_status seams_rings_find(seams_system *sys, double cutoff, int type, int graph, int max_depth, seams_rings **out)

Primitive rings (Franzblau) up to a maximum size.

sys Frame handle. cutoff Positive finite cutoff, Angstrom. type Species filter. graph SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL. max_depth Largest ring size, from 3 through 12. out Receives a new handle. On failure \*out is NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller frees \*out with seams_rings_destroy(). Thread safety Exclusive use of sys. Example seams_rings\*rings=NULL; if(seams_rings_find(sys,3.5,0,SEAMS_GRAPH_CUTOFF,6,&rings)==SEAMS_OK){ seams_rings_destroy(rings); }

.. _capi-seams_rings_destroy:

seams_rings_destroy
-------------------


.. code-block:: c

   void seams_rings_destroy(seams_rings *rings)

Free a ring set.

rings Handle from seams_rings_find(), or NULL. NULL is ignored. Ownership Releases the handle and every ring it holds. Thread safety The caller must not use rings after this returns. Example seams_rings_destroy(rings);

.. _capi-seams_rings_count:

seams_rings_count
-----------------


.. code-block:: c

   int seams_rings_count(const seams_rings *rings)

Number of primitive rings, as an int.

ABI 1 query. seams_rings_nrings() reports errors. rings Ring handle, or NULL. The ring count, or 0 when rings is NULL. Thread safety Exclusive use of rings. Example intnr=seams_rings_count(rings);

.. _capi-seams_rings_nrings:

seams_rings_nrings
------------------


.. code-block:: c

   seams_status seams_rings_nrings(const seams_rings *rings, size_t *n)

Number of primitive rings.

rings Ring handle. n Receives the ring count. SEAMS_OK, or SEAMS_ERR_NULL when rings or n is NULL. Thread safety Exclusive use of rings. Example size_tnr=0; seams_rings_nrings(rings,&nr);

.. _capi-seams_rings_histogram:

seams_rings_histogram
---------------------


.. code-block:: c

   seams_status seams_rings_histogram(const seams_rings *rings, int *counts, int counts_len)

Count primitive rings by size, with an int length.

ABI 1 entry. New callers use seams_rings_histogram_n(). rings Ring handle. counts Caller buffer. counts[k] receives the number of rings of size k. Sizes outside the buffer are ignored. counts_len Length of counts. Must be positive. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE. Ownership The caller owns counts. Thread safety Exclusive use of rings. Example inthist[13]; seams_rings_histogram(rings,hist,13);

.. _capi-seams_rings_histogram_n:

seams_rings_histogram_n
-----------------------


.. code-block:: c

   seams_status seams_rings_histogram_n(const seams_rings *rings, int *counts, size_t counts_len)

Count primitive rings by size into a buffer of known length.

rings Ring handle. counts Caller buffer. counts[k] receives the number of rings of size k. Sizes outside the buffer are ignored. counts_len Length of counts. Must be positive. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when counts_len is 0. Ownership The caller owns counts. Thread safety Exclusive use of rings. Example inthist[13]; seams_rings_histogram_n(rings,hist,13);

.. _capi-seams_rings_labels:

seams_rings_labels
------------------


.. code-block:: c

   seams_status seams_rings_labels(const seams_rings *rings, int *sizes, int *symmetry)

Size and symmetry of each ring, without buffer lengths.

ABI 1 entry. Each non-NULL buffer has length seams_rings_count(rings). New callers use seams_rings_labels_n(), which checks the lengths. rings Ring handle. sizes One int per ring, or NULL. symmetry One seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns both arrays. Thread safety Exclusive use of rings. The first symmetry request fills a cache on the handle. Example seams_rings_labels(rings,sizes,symmetry);

.. _capi-seams_rings_labels_n:

seams_rings_labels_n
--------------------


.. code-block:: c

   seams_status seams_rings_labels_n(const seams_rings *rings, int *sizes, size_t sizes_len, int *symmetry, size_t symmetry_len)

Size and symmetry of each ring into buffers of known length.

rings Ring handle. sizes One int per ring, or NULL. sizes_len Length of sizes. Ignored when sizes is NULL. symmetry One seams_ring_symmetry code per ring, or NULL. The two pointers may not both be NULL. symmetry_len Length of symmetry. Ignored when symmetry is NULL. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE when a buffer is shorter than the ring count, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. Ownership The caller owns both arrays. Thread safety Exclusive use of rings. Example size_tnr=0; seams_rings_nrings(rings,&nr); int\*sym=malloc(nr\*sizeof\*sym); seams_rings_labels_n(rings,NULL,0,sym,nr);

.. _capi-seams_rings_vertices:

seams_rings_vertices
--------------------


.. code-block:: c

   seams_status seams_rings_vertices(const seams_rings *rings, int index, int *atoms, size_t cap, size_t *n_out)

Atom indices of one ring.

rings Ring handle. index Ring index, from 0 to seams_rings_count(rings) - 1. atoms Caller buffer, or NULL to query. cap Capacity of atoms. Ignored when atoms is NULL. n_out Receives the vertex count. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_RANGE. Ownership The caller owns atoms. Indices are positions in the frame. Thread safety Exclusive use of rings. Example intatoms[12]; size_tlen=0; seams_rings_vertices(rings,0,atoms,12,&len);

.. _capi-seams_rings_fraction:

seams_rings_fraction
--------------------


.. code-block:: c

   seams_status seams_rings_fraction(const seams_rings *rings, int size, double *fraction)

Fraction of atoms that sit on a ring.

rings Ring handle. size Ring size to count. Zero counts atoms that sit on any ring. fraction Receives a value in [0, 1]. SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_NOMEM. Thread safety Exclusive use of rings. Example doublef6=0; seams_rings_fraction(rings,6,&f6);

.. _capi-seams_ring_symmetry_name:

seams_ring_symmetry_name
------------------------


.. code-block:: c

   const char * seams_ring_symmetry_name(int symmetry)

Name of a ring-symmetry code.

symmetry A seams_ring_symmetry value. Static string. Unknown codes return "unknown". Ownership The caller does not free it. Thread safety Pure. Example printf("%s\\n",seams_ring_symmetry_name(SEAMS_RING_DDC));

