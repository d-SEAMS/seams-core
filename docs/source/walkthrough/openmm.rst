OpenMM
======

``DseamsReporter`` in ``adapters/python/dseams_insitu/openmm.py`` asks
the context for host positions every step. The walkthrough builds the
8-atom cell from ``engines/openmm/test_reporter.py`` and takes two
steps with the reporter's interval set to 1.

.. code-block:: text

   openmm calls 2 atoms 8 cubic 8

Both reports are entirely cubic. The docs job fails if OpenMM cannot
be imported.

.. figure:: figures/openmm.svg
   :alt: OpenMM reporter: two calls on an 8-atom cubic cell.

   Two reporter calls on the test cell.

What this tells you
-------------------

The reporter is not a file reader. It pulls the periodic box from the
OpenMM context twice, and both pulls classify every atom as cubic.
