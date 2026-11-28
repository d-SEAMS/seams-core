===================
Rings on one thread
===================


``primitive::ringNetwork`` against OVITO 3.17.2 ``FindRingsModifier``.
The cutoff graph is built before the timer on both sides.
The cutoff is 3.5 Å and the largest ring has 6 members.
Eight timed calls follow two warmups. The median of those eight
is the mean of the two central samples.
The table rounds each median to 0.001 ms.

Frames
------

.. table::

    +----------+-----------+------------------------------------------------------------+
    | name     | molecules | frame                                                      |
    +==========+===========+============================================================+
    | tri-136  |       136 | ``input/traj/genice_sH.lammpstrj``, guest-free structure H |
    +----------+-----------+------------------------------------------------------------+
    | orth-4k  |      4096 | ``input/traj/mW_cubic.lammpstrj``, cubic ice               |
    +----------+-----------+------------------------------------------------------------+
    | orth-1e5 |    110592 | that cubic cell, replicated 3 by 3 by 3                    |
    +----------+-----------+------------------------------------------------------------+

The driver invokes the shared ``input/fetch.py`` helper when the cubic
trajectory is absent. The helper downloads the trajectory from
``d-SEAMS/seams-core-inputs`` and verifies its SHA-256 digest.

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

Measured times
--------------

The build is a release build with OpenMP on 8 cores. OVITO is 3.17.2.

.. table::

    +-----------+---------+--------------+------------+
    | molecules | threads | d-SEAMS (ms) | OVITO (ms) |
    +===========+=========+==============+============+
    |       136 |       1 |        0.183 |      0.534 |
    +-----------+---------+--------------+------------+
    |      4096 |       1 |        2.363 |     13.354 |
    +-----------+---------+--------------+------------+
    |    110592 |       1 |       68.773 |    415.508 |
    +-----------+---------+--------------+------------+
    |       136 |       8 |        0.187 |      0.541 |
    +-----------+---------+--------------+------------+
    |      4096 |       8 |        3.601 |     13.230 |
    +-----------+---------+--------------+------------+
    |    110592 |       8 |       19.874 |    358.581 |
    +-----------+---------+--------------+------------+

The one-thread medians are lower for d-SEAMS at each of the three
sizes, and so are the eight-thread medians.

Search
------

The timed call asks for rings of at most 6 members. A strict shortcut
on a ring of at most 7 members is a walk of one or two edges, read
from the neighbour rows, so this search never builds a distance ball.
The rings are the Franzblau set.
