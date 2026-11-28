<!-- https://docs.dseams.info/science/index.md -->

# Science tutorials {#science-tutorials}

Three calculations from the papers: an ice-water interface, a brine
freezing front, and the sII and sH cages. The printed census is the
check. The figures are those same frames. The gallery collects the
four pictures.

**[Ice/liquid interface and nucleation](interface.md)**

`input/traj/interface_f0.lammpstrj` is 864 waters in a cell 23 by
22 by 53 angstrom. Published CHILL+ at 3.5 angstrom reports 352
hexagonal, 434 liquid, 49 interfacial, 25 interfacial clathrate, 3
clathrate, and 1 cubic. The nucleus is frame 1 of
`input/traj/nucleation_f1.lammpstrj` (4096 atoms): 47 cubic, 62
hexagonal, 3589 liquid, and 240 interfacial.

**[Brine freezing front](brine.md)**

The 230 K frames with 20 ion pairs come from the
`d-SEAMS/brine-front-230K` dataset. All 1496 oxygens match the
stored CHILL+ names. At 5 ns the count is 741 hexagonal, 621 liquid,
89 interfacial, 41 interfacial clathrate, 4 clathrate, and no cubic.
The rejection coefficient is 0.

**[sII and sH hydrate cages](cages.md)**

`input/traj/genice_sII.lammpstrj` and
`input/traj/genice_sH.lammpstrj` are guest-free hosts, 136 waters
each. sII holds 16 of $5^{12}$ and 8 of $5^{12}6^{4}$.
sH holds 12 of $5^{12}$, 8 of $4^{3}5^{6}6^{3}$, and 4
of $5^{12}6^{8}$.

**[Science figures](gallery.md)**

The four molecular views, with the census repeated in the caption.
