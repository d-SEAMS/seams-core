GROMACS
=======

``engines/gromacs/seams_ice.cpp`` reads a trajectory in
``analyzeFrame`` and multiplies coordinates by 10, from nanometres to
angstrom. ``engines/gromacs/ice.gro`` stores the second oxygen at
0.1588 nm.

.. code-block:: text

   gromacs atoms 64 scaled_A 1.5880 ice 64

The walkthrough checks that 0.1588 nm becomes 1.5880 angstrom before
CHILL+, and that the converted cell is 64 cubic molecules. When
``gmx-seams`` is built it also has to print ``gromacs frames 2 ice 128
atoms 128``, because the tool reads the gro twice.

.. figure:: figures/gromacs.svg
   :alt: GROMACS trajectory read: 64 atoms, 0.1588 nm scaled to 1.588 angstrom.

   ``analyzeFrame`` converts the gro from nanometres.

What this tells you
-------------------

The trajectory hook's unit conversion is what makes 0.1588 nm an
oxygen-oxygen contact. With that factor the frame is 64 cubic
molecules.
