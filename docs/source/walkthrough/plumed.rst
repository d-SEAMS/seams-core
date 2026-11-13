PLUMED
======

``engines/plumed/ice.xyz`` is the cubic cell in nanometres, the frame
``plumed driver`` reads with ``engines/plumed/plumed.dat``. The
walkthrough converts it to angstrom and classifies it.

.. code-block:: text

   plumed cubic 64 six-rings 128 q6 0.628539

.. literalinclude:: ../../examples/engine_walkthrough.py
   :language: python
   :start-after: def read_xyz
   :end-before: def diamond

.. figure:: figures/plumed.svg
   :alt: PLUMED cubic-ice census, 64 cubic molecules and 128 six-rings.

   Census of ``engines/plumed/ice.xyz``.

What this tells you
-------------------

The collective variable ``SEAMS`` on this frame should report 64 cubic
molecules, 128 six-rings, and q6 0.628539. Those are the values printed
into ``COLVAR`` when the plugin is loaded.
