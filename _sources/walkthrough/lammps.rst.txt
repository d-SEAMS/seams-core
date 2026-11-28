LAMMPS
======

``engines/lammps/in.dseams`` installs ``fix ext`` as
``external pf/callback 2 2`` and then ``run 8``. The hook in
``engines/lammps/run_hook.cpp`` classifies the cell on step 0 and every
second step after that, and adds the cubic count from each call.
``engines/lammps/ice.data`` is the cell the deck reads.

.. code-block:: text

   lammps calls 5 stride 2 ice 320

Five calls times 64 cubic molecules is 320. When ``run_hook`` is in the
build, the walkthrough runs that binary on the deck and requires the
same line.

.. figure:: figures/lammps.svg
   :alt: LAMMPS fix callback: 5 calls, stride 2, ice 320.

   The fix-external hook on ``in.dseams``.

What this tells you
-------------------

The callback fires five times across ``run 8`` because the stride is 2
and step 0 is included. Each call sees 64 cubic molecules, so the hook
accumulates ice 320. That total is the fix, not a single static frame.
