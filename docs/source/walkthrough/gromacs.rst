GROMACS
=======

``engines/gromacs/ice.gro`` stores the same 64-atom cell in nanometres.
The walkthrough converts the coordinates to angstrom and classifies them.

.. code-block:: text

   gromacs cubic 64 six-rings 128 q6 0.628539

.. literalinclude:: ../../examples/engine_walkthrough.py
   :language: python
   :start-after: def read_gro
   :end-before: def read_xyz

.. figure:: figures/gromacs.svg
   :alt: GROMACS cubic-ice census, 64 cubic molecules and 128 six-rings.

   Census of ``engines/gromacs/ice.gro``.

What this tells you
-------------------

After the nanometre-to-angstrom conversion the GROMACS frame is the
LAMMPS cell. CHILL+ reports 64 cubic molecules, 128 six-rings, and mean
q6 0.628539, so ``gmx-seams`` and the library are reading one structure.
