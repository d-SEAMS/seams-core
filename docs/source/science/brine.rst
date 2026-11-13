Brine freezing front
====================

The frames are the 20-pair loading of
`d-SEAMS/brine-front-230K <https://huggingface.co/datasets/d-SEAMS/brine-front-230K/tree/4828e33ecf419c70643b2d466bfa4f283b4138dc>`_
at revision ``4828e33ecf419c70643b2d466bfa4f283b4138dc``. The tutorial
fetches ``T230_m20_r0_t0.xyz`` and ``T230_m20_r0_tend.xyz``, checks
their SHA-256, and caches them under ``docs/.cache``. A later run uses
the cache. If the fetch fails and the cache is empty, the tutorial
prints

.. code-block:: text

   brine tutorial skipped: cannot fetch d-SEAMS/brine-front-230K at 4828e33ecf419c70643b2d466bfa4f283b4138dc

and returns success. It does not treat a missing download as a passed
measurement.

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
   :alt: 5 ns brine front with a colour legend. Hexagonal ice is orange, the interface is yellow, liquid is a faint grey, clathrate is green, interfacial clathrate is pink, sodium is violet, and chloride is black.

   The 5 ns frame. Hydrogens are omitted. The cell origin is the lowest hexagonal oxygen, and the coordinates are then wrapped, so the slab is one piece and the periodic image of the brine is not left under the ice. The frame keeps one sphere radius outside the outline, including above the top edge. Liquid oxygens are small and faint. Sodium and chloride use the oxygen scale: Na is 1.0 Å and Cl is 1.8 Å where oxygen is 1.4 Å, drawn at 0.90 for the ice oxygen. The legend writes Na⁺ and Cl⁻ on two balanced rows. Brackets with end ticks mark the brine and the ice Ih slab, with a gap at the front. A dashed line across the cell is the growth front. The axis reads ``depth from lower ice surface (nm)``. The 1 nm bar sits below the legend. CHILL+ on the whole cell is 741 hexagonal, 621 liquid, 89 interfacial, 41 interfacial clathrate, 4 clathrate and no cubic. The picture keeps those non-hexagonal colours only within 1.5 nm of the growth front: 42 interfacial, 37 interfacial clathrate and all 4 clathrate. Farther from the front they are drawn as liquid. Sodium and chloride are not recoloured. The engine and the frame column agree. Sodium is violet and interfacial clathrate is pink. Those two are 31 CIEDE2000 apart. Chloride is black. ``feat/brine-science`` at ``80201fa2`` imports ``scripts/chill_palette.py``. ``paper/data/brine-frames/render_front.py`` in this tree does too.

Paper Fig 18 and SI Figs 25-26 call the same renderer:

.. code-block:: text

   python docs/examples/render_brine.py --xyz T230_m20_r0_tend.xyz --out brine.png --width-mm 85

``--width-mm`` is the atom panel at 300 ppi. ``--out`` is the png, which
also includes the depth axis and the legend.

What this tells you
-------------------

Over 5 ns at 230 K the hexagonal count moves from 733 to 741 and the
ice run stays near 25.5 angstrom. The bulk ice holds no sodium, so the
rejection coefficient is 0: the salt stays in the liquid. The stored
labels and a fresh CHILL+ call name every oxygen the same way, which is
the check that the front is the paper's front and not a second
classifier.
