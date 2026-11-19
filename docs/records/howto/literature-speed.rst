.. list-table::
   :header-rows: 0

   * - ``t1`` and ``t8`` are the medians in the threads=1 and threads=8 columns.
       ``ratio`` is ``t1 / t8``.
       A PLUMED row was measured on one thread because ``PLUMED_NUM_THREADS`` was unset.
       Both of its columns are that measurement.

One thread against eight
^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: code, operation, molecules, t1_ms, t8_ms, ratio, note

   freud, q4, 110592, 276.580, 51.302, 5.39, 
   freud, q4, 884736, 2358.091, 431.376, 5.47, 
   freud, q6, 110592, 362.326, 68.760, 5.27, 
   freud, q6, 884736, 2973.372, 550.161, 5.40, 
   OVITO, CHILL+, 110592, 622.424, 88.685, 7.02, 
   OVITO, CHILL+, 884736, 4890.789, 666.288, 7.34, 
   OVITO, PTM, 110592, 2790.624, 370.341, 7.54, 
   OVITO, PTM, 884736, 21304.364, 2849.382, 7.48, 
   OVITO, rings, 110592, 583.334, 354.605, 1.65, 
   OVITO, rings, 884736, 4883.899, 2582.748, 1.89, 
   mdapy, CHILL+, 110592, 61.938, 12.379, 5.00, 
   mdapy, CHILL+, 884736, 744.487, 125.578, 5.93, 
   mdapy, q4, 110592, 64.687, 12.400, 5.22, 
   mdapy, q4, 884736, 807.101, 96.238, 8.39, 
   mdapy, q6, 110592, 79.232, 14.572, 5.44, 
   mdapy, q6, 884736, 921.574, 107.347, 8.59, 
   mdapy, rings, 110592, 2276.747, 375.967, 6.06, 
   mdapy, rings, 884736, 19221.718, 3361.913, 5.72, 
   vesin, neighbours, 110592, 51.194, 22.269, 2.30, 
   vesin, neighbours, 884736, 1112.702, 196.866, 5.65, 
   rings-statistics-matter, rings, 110592, 1151.862, 264.953, 4.35, 
   rings-statistics-matter, rings, 884736, 11504.451, 2665.750, 4.32, 
   PLUMED, q4, 110592, 1310.025, 1316.423, 1.00, one thread; PLUMED_NUM_THREADS unset
   PLUMED, q4, 884736, 10053.680, 10122.158, 0.99, one thread; PLUMED_NUM_THREADS unset
   PLUMED, q6, 110592, 1685.454, 1680.281, 1.00, one thread; PLUMED_NUM_THREADS unset
   PLUMED, q6, 884736, 12947.767, 12818.837, 1.01, one thread; PLUMED_NUM_THREADS unset

