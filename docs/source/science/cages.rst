sII and sH hydrate cages
=========================

``input/traj/genice_sII.lammpstrj`` and ``input/traj/genice_sH.lammpstrj``
are guest-free GenIce hosts, 136 waters each. sH keeps the triclinic
tilt. Closed cages are counted with ``seams_cages_signature`` at 3.5
angstrom on the cutoff graph. The signatures are ``512`` for
:math:`5^{12}`, ``51264`` for :math:`5^{12}6^{4}`, ``sh`` for
:math:`4^{3}5^{6}6^{3}`, and ``51268`` for :math:`5^{12}6^{8}`.

.. code-block:: text

   cages sII 512 16
   cages sII 51264 8
   cages sH 512 10
   cages sH sh 4
   cages sH 51268 2

.. literalinclude:: ../../examples/science_tutorials.py
   :language: python
   :start-after: def signature
   :end-before: def read_xyz

.. literalinclude:: ../../examples/science_tutorials.c
   :language: c
   :start-after: static int cages_of
   :end-before: static int read_xyz

.. figure:: figures/cages.svg
   :alt: Cage counts for the GenIce sII and sH hosts.

   Cage counts for the GenIce sII and sH hosts.

What this tells you
-------------------

One sII cell holds 16 pentagonal dodecahedra and 8 of :math:`5^{12}6^{4}`.
Four sH formula units hold 12 of :math:`5^{12}`, 8 of
:math:`4^{3}5^{6}6^{3}`, and 4 of :math:`5^{12}6^{8}`. Those are the
cage numbers that separate the two hydrates. A count that dropped the
sheared sH edges would miss the cages whose bonds cross the xy tilt.
