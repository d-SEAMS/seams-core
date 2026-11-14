Order parameters
================


CHILL, CHILL+, and Steinhardt parameters.

Generated from the comments on the ``SEAMS_API`` declarations in ``include/seams.h``. This tree has 69 exports. ``meson test docs_capi_symbols`` fails if a declaration is missing from these pages, and ``meson test docs_capi_examples`` compiles the example when this page has one.

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
    Positive cutoff, Angstrom.

``type``
    Species filter. Zero and negative keep every atom.

``graph``
    SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.

``rule``
    SEAMS_CHILL_PUBLISHED or SEAMS_CHILL_OVITO.

``labels``
    Caller buffer of labels_len ints.

``labels_len``
    Must be at least the frame count.

``counts9``
    Nine ints, or NULL to skip the histogram.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

Published CHILL is Moore's four-neighbour rule with inclusive c_ij edges. OVITO mode uses the same thresholds and the open edges. Labels match seams_label. counts9, when non-NULL, receives nine population counts in label order.

Ownership: the caller owns labels and counts9. Thread safety: exclusive use of sys.

.. _capi-seams_chill_plus:

seams_chill_plus
----------------


.. code-block:: c

   seams_status seams_chill_plus(seams_system *sys, double cutoff, int type, int graph, int rule, int *labels, int counts9[9])

CHILL+ labels.

``sys``
    Frame handle.

``cutoff``
    Positive cutoff, Angstrom.

``type``
    Species filter. Zero and negative keep every atom.

``graph``
    SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.

``rule``
    SEAMS_CHILL_PUBLISHED or SEAMS_CHILL_OVITO.

``labels``
    Caller buffer of one label per atom.

``counts9``
    Nine ints, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

rule selects the published Nguyen and Molinero table or the OVITO Chill+ modifier. counts9, when non-NULL, receives nine population counts. labels has one entry per atom in the frame; the length is the frame count.

Ownership: the caller owns labels and counts9. Thread safety: exclusive use of sys.

.. _capi-seams_ql:

seams_ql
--------


.. code-block:: c

   seams_status seams_ql(seams_system *sys, double cutoff, int type, int order_l, int graph, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len)

Local Steinhardt q_l, and the Lechner-Dellago average on request.

``sys``
    Frame handle.

``cutoff``
    Positive cutoff, Angstrom.

``type``
    Species filter. Zero and negative keep every atom.

``order_l``
    Spherical-harmonic degree.

``graph``
    SEAMS_GRAPH_CUTOFF or SEAMS_GRAPH_MUTUAL.

``ql``
    Local q_l, or NULL.

``ql_len``
    Length of ql. Ignored when ql is NULL.

``qlbar``
    Averaged q_l, or NULL to skip the average.

``qlbar_len``
    Length of qlbar. Ignored when qlbar is NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_ARG, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

order_l is 3, 4, 6, 8, or 12. Orders 3, 4, and 6 are q3, q4, and q6. Either array may be NULL. A non-NULL array has one double per atom and its length must be at least the frame count. The average is the neighbour mean of q_lm. The cutoff graph reads the cutoff shell. The mutual graph reads the four nearest neighbours inside cutoff + 1.5.

Ownership: the caller owns both arrays. Thread safety: exclusive use of sys. The spherical-harmonic scratch is process-global, so this call does not run beside seams_q4q6() or seams_q6bar().

.. _capi-seams_q4q6:

seams_q4q6
----------


.. code-block:: c

   seams_status seams_q4q6(seams_system *sys, double cutoff, int type, double *q4, double *q6)

Local q4 and q6 on the cutoff shell.

``sys``
    Frame handle.

``cutoff``
    Positive cutoff, Angstrom.

``type``
    Species filter.

``q4``
    Local q4, or NULL.

``q6``
    Local q6, or NULL.

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL.

ABI 1 pair. Either pointer may be NULL, and the two arrays may not both be NULL. Each non-NULL array has one double per atom. Degrees 4 and 6 use the solid-harmonic kernel. The Lechner-Dellago averages are seams_ql() with a non-NULL qlbar, and the nucleation summary is seams_q6bar().

Ownership: the caller owns q4 and q6. Thread safety: exclusive use of sys. Shares the process-global scratch with seams_ql() and seams_q6bar().

.. _capi-seams_q6bar:

seams_q6bar
-----------


.. code-block:: c

   seams_status seams_q6bar(seams_system *sys, double cutoff, int type, double ice_cut, double *q6bar, size_t q6bar_len, double *mean, double *median, int *n_ice)

Lechner-Dellago q6-bar and the ice count above a cut.

``sys``
    Frame handle.

``cutoff``
    Positive cutoff, Angstrom.

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

Returns: SEAMS_OK, or SEAMS_ERR_NULL, SEAMS_ERR_STATE, SEAMS_ERR_RANGE, SEAMS_ERR_NOMEM, SEAMS_ERR_INTERNAL. At least one of q6bar, mean, median, and n_ice must be non-NULL.

The average is the same neighbour mean of q_6m as seams_ql() of order 6 on the cutoff graph. n_ice counts finite values strictly above ice_cut. The nucleation census uses 0.36. mean and median ignore non-finite values. An empty finite set stores 0.

Ownership: the caller owns q6bar. Thread safety: exclusive use of sys. Shares the process-global scratch with seams_ql() and seams_q4q6().

.. _capi-seams_steinhardt:

seams_steinhardt
----------------


.. code-block:: c

   seams_status seams_steinhardt(seams_system *sys, double cutoff, int type, int order_l, double *ql, size_t ql_len, double *qlbar, size_t qlbar_len)

Local Steinhardt q_l on the cutoff shell, and the average on request.

order_l is 3, 4, 6, 8, or 12. qlbar NULL skips the Lechner-Dellago average. This is seams_ql() with the cutoff graph. Thread safety: exclusive use of sys. Shares the process-global scratch with seams_ql(), seams_q4q6(), and seams_q6bar().

