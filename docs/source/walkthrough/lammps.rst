LAMMPS
======

``engines/lammps/ice.data`` is the 64-atom cubic cell the in-run hook
reads. ``docs/examples/engine_walkthrough.py`` sends those coordinates
through ``libseams`` at 3.5 angstrom.

.. code-block:: text

   lammps cubic 64 six-rings 128 q6 0.628539

.. literalinclude:: ../../examples/engine_walkthrough.py
   :language: python
   :start-after: def read_lammps
   :end-before: def read_gro

.. figure:: figures/lammps.svg
   :alt: LAMMPS cubic-ice census, 64 cubic molecules and 128 six-rings.

   Census of ``engines/lammps/ice.data``.

What this tells you
-------------------

Every molecule in the cubic cell is cubic ice. The ring histogram holds
128 six-rings, two per molecule, and the mean local q6 is 0.628539.
That is the number the hook is expected to print when ``run_hook`` is
built against this data file.
