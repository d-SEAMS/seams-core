Order parameters
================


CHILL, CHILL+, and Steinhardt parameters.

The stable C interface provides 85 exported functions through ``include/seams.h``.

.. literalinclude:: ../../../docs/examples/capi_order.c
   :language: c

.. _capi-seams_chill:

seams_chill
-----------


.. code-block:: c

   seams_status seams_chill(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, size_t labels_len, int counts9[9])

Original CHILL labels, published or OVITO edges.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter. Zero and negative keep every atom.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

``rule``
    ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.

``labels``
    Caller buffer of labels_len ints.

``labels_len``
    Must be at least the frame count.

``counts9``
    Nine ints that receive the population of each label, or NULL to skip the histogram.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

Published CHILL is Moore's four-neighbour rule with inclusive c_ij edges. OVITO mode uses the same thresholds and the open edges. Labels match ::seams_label.

@par Handle access Non-const: writes each atom's CHILL+ label into the stored frame and clears the handle's error slot on entry. Coordinates, cell, types, and ids are not changed. @par Ownership The caller owns labels and counts9. @par Thread safety Exclusive use of sys. @par Example @code int counts[9]; seams_chill(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n, counts); @endcode

.. _capi-seams_chill_plus:

seams_chill_plus
----------------


.. code-block:: c

   seams_status seams_chill_plus(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, int counts9[9])

CHILL+ labels, without a buffer length.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter. Zero and negative keep every atom.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

``rule``
    ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.

``labels``
    Caller buffer of one label per atom.

``counts9``
    Nine ints, or NULL.

Returns: Same codes as seams_chill_plus_n().

ABI 1 entry. labels must hold one int per atom in the frame. New callers use seams_chill_plus_n(), which checks the length.

@par Handle access Non-const: writes each atom's CHILL+ label into the stored frame and clears the handle's error slot on entry. Coordinates, cell, types, and ids are not changed. @par Ownership The caller owns labels and counts9. @par Thread safety Exclusive use of sys. @par Example @code seams_chill_plus(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, NULL); @endcode

.. _capi-seams_ql:

seams_ql
--------


.. code-block:: c

   seams_status seams_ql(seams_system *sys, double cutoff, int type, int order_l, int graph, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len)

Local Steinhardt q_l, and the Lechner-Dellago average on request.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter. Zero and negative keep every atom.

``order_l``
    Spherical-harmonic degree.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

``ql``
    Local q_l, or NULL.

``ql_len``
    Length of ql. Ignored when ql is NULL.

``qlbar``
    Averaged q_l, or NULL to skip the average.

``qlbar_len``
    Length of qlbar. Ignored when qlbar is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_ARG, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

order_l is 3, 4, 6, 8, or 12. Either array may be NULL. A non-NULL array has one double per atom and its length must be at least the frame count. The average is the neighbour mean of q_lm. The cutoff graph reads the cutoff shell. The mutual graph reads the four nearest neighbours inside cutoff + 1.5.

@par Handle access Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time. @par Ownership The caller owns both arrays. @par Thread safety Exclusive use of sys. Any q call on another handle may run beside it. @par Example @code seams_ql(sys, 3.5, 0, 6, SEAMS_GRAPH_CUTOFF, q6, n, q6bar, n); @endcode

.. _capi-seams_q4q6:

seams_q4q6
----------


.. code-block:: c

   seams_status seams_q4q6(seams_system *sys, double cutoff, int type, double *q4, double *q6)

Local q4 and q6 on the cutoff shell, without buffer lengths.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter.

``q4``
    Local q4, or NULL.

``q6``
    Local q6, or NULL. q4 and q6 may not both be NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

ABI 1 pair. Each non-NULL array has one double per atom. New callers use seams_q4q6_n(), which checks the lengths.

@par Handle access Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time. @par Ownership The caller owns q4 and q6. @par Thread safety Exclusive use of sys. Any q call on another handle may run beside it. @par Example @code seams_q4q6(sys, 3.5, 0, q4, q6); @endcode

.. _capi-seams_q6bar:

seams_q6bar
-----------


.. code-block:: c

   seams_status seams_q6bar(seams_system *sys, double cutoff, int type, double ice_cut, double *q6bar, size_t q6bar_len, double *mean, double *median, int *n_ice)

Lechner-Dellago q6-bar and the ice count above a cut.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter.

``ice_cut``
    Finite threshold on q6-bar.

``q6bar``
    Per-atom average, or NULL to keep only the summary.

``q6bar_len``
    Length of q6bar. Ignored when q6bar is NULL.

``mean``
    Receives the mean, or NULL.

``median``
    Receives the median, or NULL.

``n_ice``
    Receives the count above ice_cut, or NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL. At least one of q6bar, mean, median, and n_ice must be non-NULL.

The average is the same neighbour mean of q_6m as seams_ql() of order 6 on the cutoff graph. n_ice counts finite values strictly above ice_cut. The nucleation census uses 0.36. mean and median ignore non-finite values. An empty finite set stores 0.

@par Handle access Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time. @par Ownership The caller owns q6bar. @par Thread safety Exclusive use of sys. Any q call on another handle may run beside it. @par Example @code double mean = 0, median = 0; int n_ice = 0; seams_q6bar(sys, 3.5, 0, 0.36, NULL, 0, &mean, &median, &n_ice); @endcode

.. _capi-seams_steinhardt:

seams_steinhardt
----------------


.. code-block:: c

   seams_status seams_steinhardt(seams_system *sys, double cutoff, int type, int order_l, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len)

Local Steinhardt q_l on the cutoff shell, and the average on request.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter.

``order_l``
    3, 4, 6, 8, or 12.

``ql``
    Local q_l, or NULL.

``ql_len``
    Length of ql. Ignored when ql is NULL.

``qlbar``
    Averaged q_l, or NULL to skip the Lechner-Dellago average.

``qlbar_len``
    Length of qlbar. Ignored when qlbar is NULL.

Returns: Same codes as seams_ql().

This is seams_ql() with the cutoff graph.

@par Handle access Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time. @par Ownership The caller owns both arrays. @par Thread safety Exclusive use of sys. Any q call on another handle may run beside it. @par Example @code seams_steinhardt(sys, 3.5, 0, 6, q6, n, NULL, 0); @endcode

.. _capi-seams_chill_plus_n:

seams_chill_plus_n
------------------


.. code-block:: c

   seams_status seams_chill_plus_n(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, size_t labels_len, int counts9[9])

CHILL+ labels into a buffer of known length.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter. Zero and negative keep every atom.

``graph``
    ::SEAMS_GRAPH_CUTOFF or ::SEAMS_GRAPH_MUTUAL.

``rule``
    ::SEAMS_CHILL_PUBLISHED or ::SEAMS_CHILL_OVITO.

``labels``
    Caller buffer of labels_len ints.

``labels_len``
    Must be at least the frame count.

``counts9``
    Nine ints that receive the population of each label, or NULL.

Returns: ::SEAMS_OK; ::SEAMS_ERR_NULL for a NULL sys or labels; ::SEAMS_ERR_STATE when sys holds no frame; ::SEAMS_ERR_RANGE for a short labels buffer or a non-positive or non-finite cutoff; ::SEAMS_ERR_ARG for an unknown rule or graph; ::SEAMS_ERR_NOMEM; ::SEAMS_ERR_INTERNAL.

rule selects the q_lm shell. ::SEAMS_CHILL_PUBLISHED, the default, is Nguyen and Molinero equation 1: the four closest neighbours. ::SEAMS_CHILL_OVITO is every neighbour inside the cutoff, which is the OVITO Chill+ modifier. Labels match ::seams_label.

@par Handle access Non-const: writes each atom's CHILL+ label into the stored frame and clears the handle's error slot on entry. Coordinates, cell, types, and ids are not changed. @par Ownership The caller owns labels and counts9. @par Thread safety Exclusive use of sys. @par Example @code int counts[9]; seams_chill_plus_n(sys, 3.5, 0, SEAMS_GRAPH_CUTOFF, SEAMS_CHILL_PUBLISHED, labels, n, counts); @endcode

.. _capi-seams_q4q6_n:

seams_q4q6_n
------------


.. code-block:: c

   seams_status seams_q4q6_n(seams_system *sys, double cutoff, int type, double *q4, size_t q4_len, double *q6, size_t q6_len)

Local q4 and q6 on the cutoff shell into buffers of known length.

``sys``
    Frame handle.

``cutoff``
    Positive finite cutoff, Angstrom.

``type``
    Species filter.

``q4``
    Local q4, or NULL.

``q4_len``
    Length of q4. Ignored when q4 is NULL.

``q6``
    Local q6, or NULL. q4 and q6 may not both be NULL.

``q6_len``
    Length of q6. Ignored when q6 is NULL.

Returns: ::SEAMS_OK, or ::SEAMS_ERR_NULL, ::SEAMS_ERR_STATE, ::SEAMS_ERR_RANGE for a short buffer or a bad cutoff, ::SEAMS_ERR_NOMEM, ::SEAMS_ERR_INTERNAL.

@par Handle access Non-const only because the call clears the handle's error slot on entry; it reads the frame and nothing else on the handle. The q_lm work arrays belong to the calling thread, so calls on different handles may run at the same time. @par Ownership The caller owns q4 and q6. @par Thread safety Exclusive use of sys. Any q call on another handle may run beside it. @par Example @code seams_q4q6_n(sys, 3.5, 0, q4, n, q6, n); @endcode

