Brine freezing front
====================

The
`d-SEAMS/brine-front-230K <https://huggingface.co/datasets/d-SEAMS/brine-front-230K>`_
dataset provides the initial and final frames for a 20-pair loading.
Download ``T230_m20_r0_t0.xyz`` and ``T230_m20_r0_tend.xyz`` before
running the analysis.

Each file stores a CHILL+ name on every oxygen. The engine is run again
on the oxygen coordinates. All 1496 oxygens match the stored names at
both times. Ice is hexagonal or cubic. Sodium is the ion, water is
type 1, the axis is z, and the bin count is 20. The rejection
coefficient is the bulk-ice molality divided by the bulk-liquid molality.

.. code-block:: text

   brine t0 hexagonal 733 water 619 interfacial 91 cubic 1 agree 1496
   brine t0 thickness 25.5618 k 0.0
   brine tend hexagonal 741 water 621 interfacial 89 cubic 0 agree 1496
   brine tend thickness 25.3740 k 0.0

.. literalinclude:: ../../examples/science_tutorials.py
   :language: python
   :start-after: def brine_frame
   :end-before: def fetch_brine

.. literalinclude:: ../../examples/science_tutorials.c
   :language: c
   :start-after: static void one_brine
   :end-before: static int cmd_brine

.. figure:: figures/brine.png
   :alt: 5 ns brine front with a color legend. Hexagonal ice is orange, the interface is yellow, liquid is a faint gray, clathrate is green, interfacial clathrate is pink, sodium is violet, and chloride is black.

   The 5 ns frame. Hydrogens are omitted. CHILL+ reports 741 hexagonal,
   621 liquid, 89 interfacial, 41 interfacial clathrate, 4 clathrate,
   and no cubic molecules. Sodium and chloride remain in the liquid.

What this tells you
-------------------

Over 5 ns at 230 K the hexagonal count moves from 733 to 741 and the
ice run stays near 25.5 angstrom. The bulk ice holds no sodium, so the
rejection coefficient is 0: the salt stays in the liquid. The stored
labels and a fresh CHILL+ call name every oxygen the same way, which is
the check that the front is the paper's front and not a second
classifier.
