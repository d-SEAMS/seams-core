===========================
Performance, 4096 to 110592
===========================


Neighbour lists and primitive rings against freud, OVITO, mdapy,
pyscal3, vesin, IRA, PTM, cycless, and every other code that has a
neighbour or ring row in ``benchmarks/*.csv``. IRA and PTM have no
neighbour or ring row in those files. The molecule field is
``molecules``. The interval is ``q1_ms`` and ``q3_ms`` when the file
stores them, and ``iqr_ms`` when that is the interval it stores.
``n`` is the sample count in the file. A row with no time, or a
verdict of absent or DNF, is written DNF.

For the classifier, polyhedral template matching, and IRA rows,
``ratio_vs_dseams`` is the ratio of medians. ``ratio_lo`` and ``ratio_hi``
give its interval.

.. include:: performance-range-tables.rst
