.. list-table::
   :header-rows: 0

   * - Ice Ic, 884736 molecules, cutoff 3.5 angstrom. ``n`` is 8.
       Neighbour medians and ring ``total`` medians are recomputed from the sample lines.

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

   cell-list, total, 1, 8, 678.397
   cell-list, total, 2, 8, 408.562
   cell-list, total, 4, 8, 245.209
   cell-list, total, 8, 8, 163.807
   linkcell, total, 1, 8, 671.204
   linkcell, total, 2, 8, 425.875
   linkcell, total, 4, 8, 261.232
   linkcell, total, 8, 8, 177.018
   cell-list, link, 1, 8, 134.417
   cell-list, link, 2, 8, 91.491
   cell-list, link, 4, 8, 58.769
   cell-list, link, 8, 8, 45.987
   linkcell, link, 1, 8, 139.027
   linkcell, link, 2, 8, 89.277
   linkcell, link, 4, 8, 54.601
   linkcell, link, 8, 8, 36.755
   cell-list, search, 1, 8, 509.275
   cell-list, search, 2, 8, 297.029
   cell-list, search, 4, 8, 172.349
   cell-list, search, 8, 8, 112.440
   linkcell, search, 1, 8, 496.326
   linkcell, search, 2, 8, 320.089
   linkcell, search, 4, 8, 199.153
   linkcell, search, 8, 8, 134.539

