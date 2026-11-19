.. list-table::
   :header-rows: 0

   * - ``ratio`` is the threads=1 median divided by the threads=8 median.
       ``flat`` means that ratio is between 0.95 and 1.05.
       A PLUMED row marked ``one thread; PLUMED_NUM_THREADS unset`` used one thread.
       Both of its columns were measured on one thread because that variable was unset.
       The large-cell PLUMED row set ``PLUMED_NUM_THREADS``, so it is not that label.

One thread against eight
^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: code, operation, molecules, t1_ms, t8_ms, ratio, with_threads

   ClavisCage, cages, 136, 69.332, 69.657, 0.995, flat
   ClavisCage, cages, 368, 301.514, 302.471, 0.997, flat
   ClavisCage, cages, 46, 51.122, 51.206, 0.998, flat
   OVITO, CHILL+, 110592, 622.424, 88.685, 7.018, falls
   OVITO, CHILL+, 884736, 4890.789, 666.288, 7.340, falls
   OVITO, PTM, 110592, 2790.624, 370.341, 7.535, falls
   OVITO, PTM, 884736, 21304.364, 2849.382, 7.477, falls
   OVITO, rings, 110592, 583.334, 354.605, 1.645, falls
   OVITO, rings, 884736, 4883.899, 2582.748, 1.891, falls
   PLUMED, q4, 110592, 1310.025, 1316.423, 0.995, one thread; PLUMED_NUM_THREADS unset
   PLUMED, q4, 884736, 10053.680, 10122.158, 0.993, one thread; PLUMED_NUM_THREADS unset
   PLUMED, q6, 110592, 1685.454, 1680.281, 1.003, one thread; PLUMED_NUM_THREADS unset
   PLUMED, q6, 884736, 12947.767, 12818.837, 1.010, one thread; PLUMED_NUM_THREADS unset
   d-SEAMS, CHILL+, 110592, 24.102, 5.697, 4.231, falls
   d-SEAMS, CHILL+, 884736, 263.592, 47.444, 5.556, falls
   d-SEAMS, PTM, 110592, 30.475, 6.976, 4.369, falls
   d-SEAMS, PTM, 884736, 287.530, 51.470, 5.586, falls
   d-SEAMS, cages, 136, 1.805, 1.623, 1.112, falls
   d-SEAMS, cages, 368, 4.829, 3.893, 1.240, falls
   d-SEAMS, cages, 46, 0.682, 0.583, 1.170, falls
   d-SEAMS, neighbours, 110592, 28.325, 6.453, 4.389, falls
   d-SEAMS, neighbours, 884736, 282.177, 60.408, 4.671, falls
   d-SEAMS, q4, 110592, 23.747, 6.094, 3.897, falls
   d-SEAMS, q4, 884736, 281.416, 53.025, 5.307, falls
   d-SEAMS, q6, 110592, 25.832, 5.820, 4.438, falls
   d-SEAMS, q6, 884736, 303.544, 57.752, 5.256, falls
   d-SEAMS, rings, 1000, 1.019, 1.087, 0.937, rises
   d-SEAMS, rings, 10648, 8.972, 2.858, 3.139, falls
   d-SEAMS, rings, 110592, 86.531, 32.665, 2.649, falls
   d-SEAMS, rings, 21952, 18.111, 5.736, 3.157, falls
   d-SEAMS, rings, 32768, 26.227, 8.816, 2.975, falls
   d-SEAMS, rings, 4096, 3.499, 1.131, 3.094, falls
   d-SEAMS, rings, 64000, 49.326, 16.147, 3.055, falls
   d-SEAMS, rings, 884736, 765.181, 301.039, 2.542, falls
   fastatomstruct, king, 1000, 24.865, 20.803, 1.195, falls
   fastatomstruct, king, 10648, 1813.620, 306.073, 5.925, falls
   fastatomstruct, king, 110592, 572149.958, 66413.885, 8.615, falls
   fastatomstruct, king, 21952, 16849.495, 2166.371, 7.778, falls
   fastatomstruct, king, 32768, 47798.094, 5787.638, 8.259, falls
   fastatomstruct, king, 4096, 289.722, 109.916, 2.636, falls
   fastatomstruct, king, 64000, 189511.480, 20791.598, 9.115, falls
   freud, q4, 110592, 276.580, 51.302, 5.391, falls
   freud, q4, 884736, 2358.091, 431.376, 5.466, falls
   freud, q6, 110592, 362.326, 68.760, 5.269, falls
   freud, q6, 884736, 2973.372, 550.161, 5.405, falls
   mdapy, CHILL+, 110592, 61.938, 12.379, 5.003, falls
   mdapy, CHILL+, 884736, 744.487, 125.578, 5.928, falls
   mdapy, q4, 110592, 64.687, 12.400, 5.217, falls
   mdapy, q4, 884736, 807.101, 96.238, 8.387, falls
   mdapy, q6, 110592, 79.232, 14.572, 5.437, falls
   mdapy, q6, 884736, 921.574, 107.347, 8.585, falls
   mdapy, rings, 110592, 2276.747, 375.967, 6.056, falls
   mdapy, rings, 884736, 19221.718, 3361.913, 5.717, falls
   rings-statistics-matter, rings, 110592, 1151.862, 264.953, 4.347, falls
   rings-statistics-matter, rings, 884736, 11504.451, 2665.750, 4.316, falls
   vesin, neighbours, 110592, 51.194, 22.269, 2.299, falls
   vesin, neighbours, 884736, 1112.702, 196.866, 5.652, falls

884736-molecule sitting, one thread against eight
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: phase, code, operation, t1_ms, t8_ms, ratio, with_threads

   ice, OVITO, chill+, 4927.647, 652.481, 7.552, falls
   ice, OVITO, coordination, 287.170, 58.378, 4.919, falls
   ice, OVITO, ptm-default, 11787.795, 1585.135, 7.436, falls
   ice, OVITO, ptm-ice, 21901.041, 2856.862, 7.666, falls
   ice, OVITO, q4, 393.774, 69.910, 5.633, falls
   ice, OVITO, q6, 476.740, 81.451, 5.853, falls
   ice, PLUMED, q6, 15780.364, 2864.697, 5.509, falls
   ice, d-SEAMS, chill+, 235.591, 45.688, 5.157, falls
   ice, d-SEAMS, classifier, 204.002, 48.499, 4.206, falls
   ice, d-SEAMS, neighbours, 260.183, 52.710, 4.936, falls
   ice, d-SEAMS, q4, 144.664, 31.895, 4.536, falls
   ice, d-SEAMS, q4+q6, 319.346, 85.776, 3.723, falls
   ice, d-SEAMS, q6, 150.095, 32.188, 4.663, falls
   ice, d-SEAMS, rings, 662.186, 246.607, 2.685, falls
   ice, freud, neighbours, 1033.405, 366.366, 2.821, falls
   ice, freud, q4, 2164.249, 457.446, 4.731, falls
   ice, freud, q6, 2870.108, 506.237, 5.669, falls
   ice, mdapy, chill+, 580.057, 93.155, 6.227, falls
   ice, mdapy, neighbours, 469.862, 88.871, 5.287, falls
   ice, mdapy, q4, 698.390, 111.528, 6.262, falls
   ice, mdapy, q6, 803.702, 134.511, 5.975, falls
   ice, mdapy, rings, 18678.338, 3447.121, 5.419, falls
   ice, pyscal3, neighbours, 708.702, 154.994, 4.572, falls
   ice, pyscal3, q6, 956.332, 185.698, 5.150, falls
   ice, rings-statistics-matter, rings, 11267.173, 2707.952, 4.161, falls
   ice, vesin, neighbours, 1083.309, 218.651, 4.955, falls
   liquid, OVITO, chill+, 6997.797, 1035.081, 6.761, falls
   liquid, OVITO, coordination, 472.467, 89.174, 5.298, falls
   liquid, OVITO, ptm-default, 14735.477, 2190.199, 6.728, falls
   liquid, OVITO, ptm-ice, 23780.638, 3494.725, 6.805, falls
   liquid, OVITO, q4, 579.641, 99.272, 5.839, falls
   liquid, OVITO, q6, 678.298, 112.396, 6.035, falls
   liquid, PLUMED, q6, 14772.141, 3248.878, 4.547, falls
   liquid, d-SEAMS, chill+, 388.618, 70.502, 5.512, falls
   liquid, d-SEAMS, classifier, 231.571, 63.954, 3.621, falls
   liquid, d-SEAMS, neighbours, 358.812, 70.649, 5.079, falls
   liquid, d-SEAMS, q4, 286.448, 51.916, 5.518, falls
   liquid, d-SEAMS, q4+q6, 602.804, 123.123, 4.896, falls
   liquid, d-SEAMS, q6, 292.794, 53.166, 5.507, falls
   liquid, d-SEAMS, rings, 2942.271, 677.074, 4.346, falls
   liquid, freud, neighbours, 1432.248, 562.501, 2.546, falls
   liquid, freud, q4, 2704.072, 621.275, 4.352, falls
   liquid, freud, q6, 3476.079, 693.877, 5.010, falls
   liquid, mdapy, chill+, 939.828, 168.668, 5.572, falls
   liquid, mdapy, neighbours, 771.802, 112.969, 6.832, falls
   liquid, mdapy, q4, 971.259, 148.271, 6.551, falls
   liquid, mdapy, q6, 1132.685, 178.443, 6.348, falls
   liquid, mdapy, rings, 27579.640, 5137.782, 5.368, falls
   liquid, pyscal3, neighbours, 5676.788, 4853.308, 1.170, falls
   liquid, pyscal3, q6, 6518.490, 4845.223, 1.345, falls
   liquid, rings-statistics-matter, rings, 46036.912, 7477.520, 6.157, falls
   liquid, vesin, neighbours, 1208.267, 249.099, 4.851, falls

signature search, one thread against eight
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: crystal, cell, molecules, t1_ms, t8_ms, ratio, with_threads

   sI, sI 2x2x2, 368, 1.671, 0.889, 1.879, falls
   sI, sI 4x4x4, 2944, 11.155, 5.296, 2.107, falls
   sI, sI 6x6x6, 9936, 34.513, 11.407, 3.026, falls
   sI, sI 8x8x8, 23552, 79.231, 26.431, 2.998, falls
   sII, sII 2x2x2, 1088, 4.815, 2.544, 1.893, falls
   sII, sII 4x4x4, 8704, 35.908, 11.174, 3.214, falls
   sII, sII 6x6x6, 29376, 121.799, 31.094, 3.917, falls
   sII, sII 8x8x8, 69632, 292.810, 73.109, 4.005, falls

pool-sized kernels, one thread against eight
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: kernel, molecules, t1_ms, t8_ms, ratio, with_threads

   chill, 110592, 24.159, 5.275, 4.580, falls
   chill, 32768, 7.429, 1.680, 4.422, falls
   chill, 4096, 0.908, 0.349, 2.600, falls
   chill, 884736, 232.528, 39.195, 5.933, falls
   hop, 110592, 26.382, 7.194, 3.667, falls
   hop, 32768, 7.778, 2.033, 3.827, falls
   hop, 4096, 0.979, 0.381, 2.565, falls
   hop, 884736, 256.551, 62.762, 4.088, falls
   q4-cutoff, 110592, 16.288, 3.592, 4.535, falls
   q4-cutoff, 32768, 5.263, 1.213, 4.338, falls
   q4-cutoff, 4096, 0.612, 0.244, 2.510, falls
   q4-cutoff, 884736, 154.053, 25.552, 6.029, falls
   q4-local, 110592, 20.345, 4.051, 5.022, falls
   q4-local, 32768, 5.905, 1.300, 4.542, falls
   q4-local, 4096, 0.667, 0.263, 2.538, falls
   q4-local, 884736, 201.293, 40.112, 5.018, falls
   q6-cutoff, 110592, 16.812, 3.852, 4.364, falls
   q6-cutoff, 32768, 5.427, 1.316, 4.123, falls
   q6-cutoff, 4096, 0.637, 0.242, 2.632, falls
   q6-cutoff, 884736, 158.820, 26.620, 5.966, falls
   q6-local, 110592, 23.795, 4.614, 5.158, falls
   q6-local, 32768, 6.898, 1.956, 3.526, falls
   q6-local, 4096, 0.772, 0.392, 1.970, falls
   q6-local, 884736, 239.786, 51.278, 4.676, falls
   rings, 110592, 62.012, 13.538, 4.581, falls
   rings, 32768, 19.361, 4.263, 4.541, falls
   rings, 4096, 2.467, 0.634, 3.892, falls
   rings, 884736, 495.117, 112.990, 4.382, falls
   symmetry, 110592, 24.393, 6.161, 3.960, falls
   symmetry, 32768, 7.109, 1.866, 3.810, falls
   symmetry, 4096, 0.913, 0.415, 2.202, falls
   symmetry, 884736, 214.114, 48.351, 4.428, falls

Pool audit at 4096 molecules, four threads against eight
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
   :header: kernel, t4_ms, t8_ms, t8_over_t4

   rings, 0.941, 0.634, 0.674
   q4-local, 0.267, 0.263, 0.983
   q6-local, 0.380, 0.392, 1.029
   q4-cutoff, 0.243, 0.244, 1.005
   q6-cutoff, 0.230, 0.242, 1.052
   chill, 0.348, 0.349, 1.002
   symmetry, 0.415, 0.415, 1.000
   hop, 0.411, 0.381, 0.929

