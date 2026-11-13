# Ice Ih rings

`seams_rings_find` lists each Franzblau primitive ring once, then labels
it. On ideal ice Ih (oxygen bond 2.75 Å, orthorhombic cell, cutoff
3.5 Å, depth 6) every primitive ring is a hexagon, and there are two
per molecule.

The label is not another ring. A hexagon is an HC basal, an HC
prismatic, or neither. The strict prism-block label is empty on this
crystal: that test wants two rings of one size with no shared vertex
and a perfect bond matching. The cages of ice Ih are the HC labels.

`cycless` `cycles_iter` (maximum size 6, with fractional coordinates so
a cycle that winds the cell is dropped), `rings-statistics-matter`
`ring_statistics` (deduplicated, maximum level 8, `supercell` on when
the cell is narrower than the search), and OVITO `FindRingsModifier`
(maximum size 6) count those same hexagons. The network has girth 6,
so a hexagon has no chord, and this lattice does not contain the spur
that separates a chordless hexagon from a primitive one. That spur is
in `markdown/rings-theory.md`. A chordless search keeps the spur
hexagon. This search does not. The two searches agree on ideal ice Ih
and on ideal ice Ic.

Leaving the cycless coordinate filter off, or leaving
`rings-statistics-matter`'s supercell off on the 128-molecule cell,
also counts squares and extra hexagons. Those cycles wind the periodic
cell. They are not in the primitive census. OVITO FindRings does the
same on that narrow cell (32 squares and 384 hexagons). At 4096
molecules the box is wide enough, and FindRings counts the 8192
primitive hexagons.

The prism pass used to compare every pair of hexagons. A neighbour
index visits only the rings that can pass the basal test, which is the
same set of pairs. Eight threads run that pass, and the stacking-plane
pass, across the rings. The ring list and the labels do not change with
the thread count.
