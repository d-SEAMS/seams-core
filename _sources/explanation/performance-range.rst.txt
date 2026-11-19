===========================
Performance, 4096 to 110592
===========================


Neighbour lists and primitive rings from the comparisons other
than the release-candidate sitting and the 884736-molecule head-to-head.
A row is listed when its stored ``n`` is at least 6. The interval is
``q1 ms`` and ``q3 ms`` when recorded, ``iqr ms`` when that is
the recorded interval, and ``min ms`` when that is the second time
recorded. The version column is the version the row recorded. Two of the tables store a variant
column, ``rings`` or ``bonds_and_rings``.

One cycless comparison also stores a 1000-molecule row. It is listed with
that table. Rows whose stored ``n`` is below 6 are omitted.

The release-candidate medians are on
`the sitting page <performance.rst>`_.
The 884736-molecule ice and liquid medians are on
`the head-to-head <performance-h2h.rst>`_.
What each comparison measures is on
`the overview <performance-overview.rst>`_.

.. include:: performance-range-tables.rst
