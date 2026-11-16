.. list-table::
   :header-rows: 0

   * - Rendered from ``benchmarks/pool-audit-v5.csv``.
       The CSV is the columnar copy of ``benchmarks/pool-audit-v5.txt``.
       Ice Ic, 884736 molecules, cutoff 3.5 angstrom. ``n`` is 8.
       Neighbour medians are recomputed from the sample lines.
       Ring ``total`` medians are recomputed the same way.
       ``cad887b2`` is the previous full ring call. ``9fcaab36`` is this sitting.

Neighbours against OVITO
^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: code, operation, threads, n, median_ms

   d-SEAMS, neighbours, 1, 8, 142.803
   OVITO, neighbours, 1, 8, 235.820
   d-SEAMS, neighbours, 8, 8, 46.160
   OVITO, neighbours, 8, 8, 52.130

Full cutoff ring call
^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: code, operation, threads, n, median_ms

   9fcaab36, total, 1, 8, 678.397
   9fcaab36, total, 2, 8, 408.562
   9fcaab36, total, 4, 8, 245.209
   9fcaab36, total, 8, 8, 163.807
   cad887b2, total, 1, 8, 671.204
   cad887b2, total, 2, 8, 425.875
   cad887b2, total, 4, 8, 261.232
   cad887b2, total, 8, 8, 177.018
   9fcaab36, link, 1, 8, 134.417
   9fcaab36, link, 2, 8, 91.491
   9fcaab36, link, 4, 8, 58.769
   9fcaab36, link, 8, 8, 45.987
   cad887b2, link, 1, 8, 139.027
   cad887b2, link, 2, 8, 89.277
   cad887b2, link, 4, 8, 54.601
   cad887b2, link, 8, 8, 36.755
   9fcaab36, search, 1, 8, 509.275
   9fcaab36, search, 2, 8, 297.029
   9fcaab36, search, 4, 8, 172.349
   9fcaab36, search, 8, 8, 112.440
   cad887b2, search, 1, 8, 496.326
   cad887b2, search, 2, 8, 320.089
   cad887b2, search, 4, 8, 199.153
   cad887b2, search, 8, 8, 134.539

