# Quasi-One-Dimensional Ice Classification

The trajectory file for this example is [here on
figshare](https://figshare.com/articles/Quasi-1D_Nanotube_LAMMPS_Trajectory/11448768).
The trajectory file details a portion of a long run of 750 molecules of
TIP4P/2005 water, inside a smooth featureless (13,0) nanotube. This nanotube
approximates a zigzag carbon nanotube. At 240 K, the temperature of the
simulation, the ice nanotube is primarily composed of tetragonal blocks. In this
example, prism blocks are identified according to a topological network
criterion for confined ice [1]. On running the example, an output top-level
directory named _runOne_ is created.

## Steps to Run the Example

In order to run this example, without making any changes to the example files, please follow the steps below.

- Download the LAMMPS trajectory file from [here on
  figshare](https://figshare.com/articles/Quasi-1D_Nanotube_LAMMPS_Trajectory/11448768).
  Copy the downloaded trajectory file, entitled _dump-240-square.lammpstrj_,
  into the _traj_ folder inside the top-level directory _input_. Alternatively,
  you could change the path to the trajectory file in the _conf.yaml_ file:

```{.lua}
trajectory: "path/to/trajectory/file"
```

- You can obtain the other input files required from _example_lua/iceNanotube_
  folder. Copy the contents of the _iceNanotube_ into the top-level _lua_inputs_
  directory.
- You can change the frames to be analyzed by updating the options in the
  _vars.lua_ file. The starting and ending frames are inclusive, starting from 1
  onwards, irrespective of the timestep number.
- A custom volume slice can also be defined in the _vars.lua_ file.
- The _functions.lua_ file actually contains the Lua functions which interface
  with the C++ backend.

## Analyzing the Output

Inside the output directory, a directory called _topoINT_ is created. Inside the
_topoINT_ directory, a file called _nPrisms.dat_ contains the number of prism
blocks and the corresponding height% [2], for each frame. Inside
_runOne/topoINT/dataFiles_, LAMMPS data files which are numbered according to
the frame number are created. These data files can be visualized in
[OVITO](https://www.ovito.org/) or [VMD](http://www.ks.uiuc.edu/Research/vmd/),
although OVITO is recommended for optimal type visualization.

## References

1. Goswami, A., & Singh, J. K. (2020). A general topological network criterion for exploring the structure of icy nanoribbons and monolayers. Physical Chemistry Chemical Physics. doi:10.1039/c9cp04902a
2. Goswami, R.; Goswami, A.; Singh, J. K. (2019). "d-SEAMS: Deferred Structural Elucidation Analysis for Molecular Simulations". arXiv:1909.09830 [physics.comp-ph].
