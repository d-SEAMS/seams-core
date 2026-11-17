===========================
Performance, 4096 to 110592
===========================

Neighbour lists and primitive rings from the committed tables other
than the release-candidate sitting and the 884736-molecule head-to-head.
A row is listed when its stored ``n`` is at least 6. The interval is
``q1 ms`` and ``q3 ms`` when the file stores them, ``iqr ms`` when that is
the interval it stores, and ``min ms`` when that is the second time it
stores. The version column is the version the row recorded. Several
watch files record 2.13.0. Two of the files store a variant
column, ``rings`` or ``bonds_and_rings``.

One cycless record also stores a 1000-molecule row. It is listed with
that record. Rows whose stored ``n`` is below 6 are omitted.

The release-candidate medians are on
`the sitting page <performance.html>`_.
The 884736-molecule ice and liquid medians are on
`the head-to-head <performance-h2h.html>`_.
What each sitting measures is
`the overview <performance-overview.html>`_.

.. include:: performance-range-tables.rst
