OpenMM
======

``engines/openmm/test_reporter.py`` builds a diamond cell with lattice
constant ``2.75 * 4 / sqrt(3)`` angstrom. The walkthrough classifies
the 2 by 2 by 2 repeat, the same 64 molecules as ``ice.data``. A
machine with OpenMM can also run ``DseamsReporter`` for one step; the
asserted census is this lattice.

.. code-block:: text

   openmm cubic 64 six-rings 128 q6 0.628539

.. literalinclude:: ../../examples/engine_walkthrough.py
   :language: python
   :start-after: def diamond
   :end-before: def figure

.. figure:: figures/openmm.svg
   :alt: OpenMM cubic-ice census, 64 cubic molecules and 128 six-rings.

   Census of the reporter's diamond cell.

What this tells you
-------------------

The reporter's cell is the same 64-molecule diamond as the LAMMPS and
GROMACS inputs. Mean q6 is 0.628539 and the six-ring count is 128, so
one OpenMM step that freezes those coordinates classifies cubic ice.
