Linkcell against the cell list
------------------------------

The rows are benchmarks/cutoff-rings.csv.
Ice Ic at 884736 molecules. Each row is the full cutoff ring call.
``n`` is 8.

.. csv-table::
   :header: code, operation, threads, median ms, iqr ms, n

   linkcell, total, 1, 653.639, 9.726, 8
   cell-list, total, 1, 654.896, 8.296, 8
   linkcell, total, 2, 521.980, 4.294, 8
   cell-list, total, 2, 418.741, 3.510, 8
   linkcell, total, 4, 323.101, 1.161, 8
   cell-list, total, 4, 252.877, 1.991, 8
   linkcell, total, 8, 228.632, 5.324, 8
   cell-list, total, 8, 173.215, 2.554, 8
