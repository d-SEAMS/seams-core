Ice/liquid interface and nucleation
===================================

The interface dump is ``input/traj/interface_f0.lammpstrj``, 864 waters
in a cell 23 by 22 by 53 angstrom. The nucleus is frame 1 of the mW
trajectory released with the paper (``input/traj/nucleation_f1.lammpstrj``,
4096 atoms, timestep 15677000). Both are labelled with published CHILL+
at 3.5 angstrom on the cutoff graph.

Layer cubicity is the cubic share of the basal layers. The period, the
majority class of each layer, and the cubicity follow
``scripts/science_render.py``. The figure itself is the v12 renderer
from ``science/interface-nucleation`` at
``f65657b9b89481e7dd6617b1987bc58b2ec0f01c``, vendored as
``scripts/science_render.py``. Liquid is blended at 0.40, basal
layers are numbered, and the image has no baked caption. Tachyon
draws the raster under Xvfb. An SVG stand-in is logged and fails the
docs job. On the slab
the stacking axis is z, the long edge. On the nucleus it is the box
axis whose crystal atoms have the strongest density wave in that
period window.

Python and C print the same census:

.. code-block:: text

   interface hexagonal 352 water 434 interfacial 49 cubic 1
   interface layers HHHHHHHHC cubicity 0.111111
   nucleus cubic 47 hexagonal 62 water 3589 interfacial 240
   nucleus axis x layers CHHHCHHHHHHHHH cubicity 0.142857

.. literalinclude:: ../../examples/science_tutorials.py
   :language: python
   :start-after: def layer_period
   :end-before: def strongest_axis

.. literalinclude:: ../../examples/science_tutorials.c
   :language: c
   :start-after: static double layer_period
   :end-before: static int cmd_interface

.. figure:: figures/interface.png
   :alt: Interface slab with basal bilayer indices, a shared colour legend, and a label on the periodic image of the interface.

   The view looks along x, so the basal stack is vertical. Each numbered tick is one crystalline basal bilayer: the two puckered halves share an index, and the centres are about 3.67 angstrom apart. A bin with fewer than half the oxygens of a full bilayer is not numbered, so the single cubic oxygen and the interfacial shell have no index. Liquid is drawn like the brine figure, at radius 0.50 and 18 percent opacity. The low-z edge is labelled ``interface (periodic image)`` because that face is the periodic copy of the high-z interface. The legend uses the same order and the same type size as the brine figure.

What this tells you
-------------------

The slab is a hexagonal ice block against liquid. 352 molecules are
hexagonal, 434 are liquid, and 49 are the interfacial shell. One
molecule is cubic. Along z the old floor bins count that molecule as
its own layer and the cubicity of those nine bins is 1/9. The figure
does not number it. The numbered index is the basal stack: six
bilayers of hexagonal ice. The nucleus in the first saved frame is
mixed: 47 cubic and 62 hexagonal molecules in a liquid of 3589.
Molecule cubicity is 47/109. Along x, two of the fourteen layers are
cubic-majority, so the layer cubicity is 2/14. That mixed stacking is
the quantity the layer panel of the paper figure is drawn to show.
