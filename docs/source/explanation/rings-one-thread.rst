===================
Rings on one thread
===================


``primitive::ringNetwork`` against OVITO 3.17.2 ``FindRingsModifier``.
The cutoff graph is built before the timer on both sides.
The cutoff is 3.5 Å and the largest ring has 6 members.
Eight timed calls follow two warmups. The median of those eight
is the mean of the two central samples.

Frames
------

These are the three sizes in the earlier one-thread campaign.

.. table::

    +----------+-----------+-------------------------------------------------+
    | name     | molecules | frame                                           |
    +==========+===========+=================================================+
    | tri-136  |       136 | ``genice_sH.lammpstrj``, guest-free structure H |
    +----------+-----------+-------------------------------------------------+
    | orth-4k  |      4096 | ``mW_cubic.lammpstrj``, cubic ice               |
    +----------+-----------+-------------------------------------------------+
    | orth-1e5 |    110592 | that cubic cell, replicated 3 by 3 by 3         |
    +----------+-----------+-------------------------------------------------+

The per-size census is the same in both codes.

.. table::

    +-----------+----+-----+--------+
    | molecules |  4 |   5 |      6 |
    +===========+====+=====+========+
    |       136 | 12 | 120 |     68 |
    +-----------+----+-----+--------+
    |      4096 |  0 |   0 |   8192 |
    +-----------+----+-----+--------+
    |    110592 |  0 |   0 | 221184 |
    +-----------+----+-----+--------+

This session
------------

Eight cores, one thread per core, Intel Xeon, family 6, model 207,
stepping 2. GCC 15.2.0, release build, OpenMP. OVITO 3.17.2.
The core count is 8. The record does not name the host.

.. table::

    +-----------+---------+--------------+------------+
    | molecules | threads | d-SEAMS (ms) | OVITO (ms) |
    +===========+=========+==============+============+
    |       136 |       1 |        0.212 |      0.550 |
    +-----------+---------+--------------+------------+
    |      4096 |       1 |        2.715 |     13.548 |
    +-----------+---------+--------------+------------+
    |    110592 |       1 |       82.343 |    474.104 |
    +-----------+---------+--------------+------------+
    |       136 |       8 |        0.090 |      0.830 |
    +-----------+---------+--------------+------------+
    |      4096 |       8 |        1.521 |     13.946 |
    +-----------+---------+--------------+------------+
    |    110592 |       8 |       27.246 |    298.562 |
    +-----------+---------+--------------+------------+

On one thread d-SEAMS is faster at each of the three sizes.
Samples are in ``markdown/rings-one-thread.json``.

Earlier campaign
----------------

That record used another 8-core Xeon and an older build of the search.
OVITO was faster there. The frames, the cutoff, the maximum ring size,
and the timing boundary (graph outside the timer) are the ones above.

.. table::

    +-----------+--------------+------------+
    | molecules | d-SEAMS (ms) | OVITO (ms) |
    +===========+==============+============+
    |       136 |        1.078 |      0.697 |
    +-----------+--------------+------------+
    |      4096 |       31.617 |     15.831 |
    +-----------+--------------+------------+
    |    110592 |      864.701 |    511.024 |
    +-----------+--------------+------------+

Search
------

The search packs adjacency into one column array. Rings of at most
seven members answer a shortcut from the neighbour rows; rings of size
eight to ten use an occupancy bitset for the bounded ball. These choices
do not change the rings returned.
