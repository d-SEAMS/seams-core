=======================
Linkcell strong scaling
=======================


An orthorhombic frame is cut by linkcell's weighted recursive coordinate
bisection of the periodic cell grid. The balance bound is the tightest a
cell-aligned cut can promise: one plus the part count times the heaviest
cell, over the total weight. The first call of an analysis weights cells
by how many points they hold. A later call of that analysis weights them
by the work measured on the first call.

The halo is painted in the same walk that visits the cells a rank owns.
A point in an owned cell stays owned. A point in a painted cell that
another rank owns is a ghost. The same cut is used inside one process.
Each thread owns a block of cells and writes the neighbour rows of the
atoms in those cells. The halo is painted in that owned walk, so the
threads do not share one cell mesh. The cut is reused on the next call
with the same grid. Neighbour rows of an unchanged frame are refilled
in place.

A skewed cell does not use that mesh. It keeps the slab, pencil, or
block cut.

Session
-------

One sitting, 12 November 2026. ``host_cores`` was read with
``OMP_NUM_THREADS`` unset and reported 8. The compiler is g++ 15.2.0. The
library is a Meson release build. Eight samples follow one warm-up. The
median of an even count is the mean of the two central samples. Speedup
is against one thread. The four operations are the neighbour rows,
CHILL+, the Steinhardt q6, and primitive rings of depth 6. The ice is
eight molecules in a cubic cell. No host name is recorded.

On 884736 molecules, eight threads build the neighbour rows in 25.7 ms,
5.72 times the one-thread time. CHILL+, q6 and the ring census use those
rows. OVITO 3.17.2 ``CreateBonds`` and freud 3.6.1 ``AABBQuery`` were timed
on the same ice in this sitting. OVITO stores one record per undirected
bond. freud stores both directions, which is the d-SEAMS pair count.

OpenMP, 110592 molecules
------------------------

.. table::

    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    | threads | neigh ms | speedup | chill ms | speedup | q6 ms | speedup | rings ms | speedup |
    +=========+==========+=========+==========+=========+=======+=========+==========+=========+
    |       1 |     15.2 |    1.00 |     47.4 |    1.00 |  43.5 |    1.00 |     83.4 |    1.00 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    |       2 |     11.6 |    1.32 |     34.8 |    1.36 |  26.4 |    1.65 |     64.0 |    1.30 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    |       4 |      7.4 |    2.05 |     28.3 |    1.67 |  16.7 |    2.60 |     34.0 |    2.46 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    |       8 |      5.6 |    2.71 |     23.8 |    1.99 |  10.7 |    4.07 |     20.2 |    4.13 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+

OpenMP, 884736 molecules
------------------------

.. table::

    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    | threads | neigh ms | speedup | chill ms | speedup | q6 ms | speedup | rings ms | speedup |
    +=========+==========+=========+==========+=========+=======+=========+==========+=========+
    |       1 |    146.9 |    1.00 |    442.6 |    1.00 | 497.9 |    1.00 |    666.9 |    1.00 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    |       2 |     74.2 |    1.98 |    285.8 |    1.55 | 318.5 |    1.56 |    406.4 |    1.64 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    |       4 |     41.4 |    3.55 |    224.2 |    1.97 | 225.7 |    2.21 |    236.4 |    2.82 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+
    |       8 |     25.7 |    5.72 |    157.8 |    2.81 | 181.0 |    2.75 |    158.9 |    4.20 |
    +---------+----------+---------+----------+---------+-------+---------+----------+---------+

Neighbours against OVITO and freud
----------------------------------

Same frames and the same sitting. Times are milliseconds.

.. table::

    +-----------+---------+---------+--------------+-------------+
    | molecules | threads | d-SEAMS | OVITO 3.17.2 | freud 3.6.1 |
    +===========+=========+=========+==============+=============+
    |    110592 |       1 |    15.2 |         43.9 |       101.7 |
    +-----------+---------+---------+--------------+-------------+
    |    110592 |       8 |     5.6 |         11.2 |       103.2 |
    +-----------+---------+---------+--------------+-------------+
    |    884736 |       1 |   146.9 |        302.9 |      1020.0 |
    +-----------+---------+---------+--------------+-------------+
    |    884736 |       8 |    25.7 |        104.1 |      1023.9 |
    +-----------+---------+---------+--------------+-------------+

MPI, one thread per rank
------------------------

These ranks were timed in the earlier sitting the same day, before each
thread owned a block. Speedup is against that sitting's one rank.

.. table::

    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    | molecules | ranks | neigh ms | speedup | chill ms | speedup | q6 ms | speedup | rings ms | speedup |
    +===========+=======+==========+=========+==========+=========+=======+=========+==========+=========+
    |    110592 |     1 |     15.4 |    1.00 |     33.2 |    1.00 |  39.8 |    1.00 |     75.1 |    1.00 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    |    110592 |     2 |     18.7 |    0.83 |     29.6 |    1.12 |  32.9 |    1.21 |     39.9 |    1.88 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    |    110592 |     4 |     10.3 |    1.49 |     17.3 |    1.92 |  24.0 |    1.66 |     21.4 |    3.51 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    |    110592 |     8 |      4.1 |    3.76 |      7.6 |    4.38 |   8.1 |    4.91 |     10.7 |    7.02 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    |    884736 |     1 |    166.4 |    1.00 |    340.2 |    1.00 | 556.7 |    1.00 |    618.4 |    1.00 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    |    884736 |     2 |    127.8 |    1.30 |    238.0 |    1.43 | 328.2 |    1.70 |    336.3 |    1.84 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    |    884736 |     4 |     62.2 |    2.68 |    115.9 |    2.94 | 156.7 |    3.55 |    168.0 |    3.68 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
    |    884736 |     8 |     26.4 |    6.30 |     59.3 |    5.74 |  67.8 |    8.22 |     84.6 |    7.31 |
    +-----------+-------+----------+---------+----------+---------+-------+---------+----------+---------+
