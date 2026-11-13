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
   cages sH 512 12
   cages sH sh 8
   cages sH 51268 4

.. literalinclude:: ../../examples/science_tutorials.py
   :language: python
   :start-after: def signature
   :end-before: def read_xyz

.. literalinclude:: ../../examples/science_tutorials.c
   :language: c
   :start-after: static int cages_of
   :end-before: static int read_xyz

.. figure:: figures/sii.png
   :alt: sII cages. 5^12 polyhedra are orange and 5^12 6^4 polyhedra are blue.

   sII. Orange is :math:`5^{12}` and blue is :math:`5^{12}6^{4}`, the hexagonal and cubic colours. The legend is titled sII and the frame includes every cage.

.. figure:: figures/sh.png
   :alt: sH cages. 5^12 is orange, the medium cage is yellow, and 5^12 6^8 is green.

   sH. Orange is :math:`5^{12}`, yellow is :math:`4^{3}5^{6}6^{3}`, and green is :math:`5^{12}6^{8}`. Each face is the ring placed on the images nearest the cage centre. A face is omitted when an edge reaches 3.5 angstrom or the planarity residual exceeds 0.35 angstrom. The legend is titled sH and the frame includes every cage.

What this tells you
-------------------

One sII cell holds 16 pentagonal dodecahedra and 8 of :math:`5^{12}6^{4}`.
Four sH formula units hold 12 of :math:`5^{12}`, 8 of
:math:`4^{3}5^{6}6^{3}`, and 4 of :math:`5^{12}6^{8}`. Those are the
cage numbers that separate the two hydrates. A count that dropped the
sheared sH edges would miss the cages whose bonds cross the xy tilt.
