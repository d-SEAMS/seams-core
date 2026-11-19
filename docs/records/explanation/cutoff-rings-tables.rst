.. list-table::
   :header-rows: 0

   * - Ice Ic at 884736 molecules. ``linkcell`` is the buffered cell-list builder.
       ``cell-list`` is the orthorhombic list inside the ring call.
       ``total`` medians are recomputed from ``samples_ms``.

Full cutoff ring call
^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: code, operation, threads, n, median_ms, iqr_ms

   cell-list, total, 1, 8, 654.896, 8.296
   cell-list, total, 2, 8, 418.741, 3.510
   cell-list, total, 4, 8, 252.877, 1.991
   cell-list, total, 8, 8, 173.215, 2.554
   linkcell, total, 1, 8, 653.639, 9.726
   linkcell, total, 2, 8, 521.980, 4.294
   linkcell, total, 4, 8, 323.101, 1.161
   linkcell, total, 8, 8, 228.632, 5.324
   cell-list, link, 1, 8, 132.279, 
   cell-list, link, 2, 8, 86.208, 
   cell-list, link, 4, 8, 52.179, 
   cell-list, link, 8, 8, 35.622, 
   linkcell, link, 1, 8, 131.875, 
   linkcell, link, 2, 8, 206.762, 
   linkcell, link, 4, 8, 137.551, 
   linkcell, link, 8, 8, 107.620, 
   cell-list, convert, 1, 8, 28.719, 
   cell-list, convert, 2, 8, 15.368, 
   cell-list, convert, 4, 8, 7.783, 
   cell-list, convert, 8, 8, 4.628, 
   linkcell, convert, 1, 8, 28.935, 
   linkcell, convert, 2, 8, 15.060, 
   linkcell, convert, 4, 8, 7.941, 
   linkcell, convert, 8, 8, 4.544, 
   cell-list, search, 1, 8, 492.962, 
   cell-list, search, 2, 8, 316.537, 
   cell-list, search, 4, 8, 192.863, 
   cell-list, search, 8, 8, 132.576, 
   linkcell, search, 1, 8, 493.166, 
   linkcell, search, 2, 8, 298.425, 
   linkcell, search, 4, 8, 177.345, 
   linkcell, search, 8, 8, 115.945, 

