<!-- https://docs.dseams.info/user/howto/confined-ice.md -->

# Confined ice {#confined-ice}

Classify ice in a nanotube or a monolayer.

A nanotube is a stack of prism blocks. `ring::prismAnalysis` in
`topo_one_dim.hpp` writes those blocks. A monolayer is a sheet of
primitive rings. `ring::polygonRingAnalysis` in `topo_two_dim.hpp`
writes the sheet. Both headers are on the
[API map](../../reference/api.md).

`seams prisms` is a separate count: perfect and imperfect prisms of
ring size 3 to 6 on the mutual four-nearest graph (`ring::findPrisms`).
Nanotube blocks and monolayer sheets come from the two writers above.

The writers are on the Python Frame in pydseamslib.
With an output directory, `Frame.find_prisms` writes the nanotube
blocks through `prismAnalysis`. `Frame.monolayer_rings` writes the
sheet through `polygonRingAnalysis`. The calls, and a run on
`input/traj/exampleTraj.lammpstrj`, are in
[Python Frame](../../reference/python-frame.md) and
[Confined water](../tutorials/confined-water.md).
Bulk CHILL+ on that dump is
[the mixed TIP4P tutorial](../tutorials/bulk-ice.md).

Lua examples live under `example_lua/iceNanotube/`
(`shapeMatching` and `strictCriterion`) in
[yodaStruct](https://github.com/d-SEAMS/yodaStruct).

The 2020 trajectories are figshare project
[d-SEAMS Datasets](https://figshare.com/projects/d-SEAMS_Datasets/73545):
[Quasi-1D Nanotube LAMMPS Trajectory](https://figshare.com/articles/Quasi-1D_Nanotube_LAMMPS_Trajectory/11448768)
and
[Monolayer LAMMPS Trajectory](https://figshare.com/articles/Monolayer_LAMMPS_Trajectory/11448741).
Download them from those records.
