# Quasi-One-Dimensional Ice Classification 

The trajectory file for this example is [here on figshare](https://figshare.com/articles/Quasi-1D_Nanotube_LAMMPS_Trajectory/11448768). The trajectory file details a portion of a long run of 750 molecules of TIP4P/2005 water, inside a smooth featureless (13,0) nanotube. This nanotube approximates a zigzag carbon nanotube. At 240 K, the temperature of the simulation, the ice nanotube is primarily composed of tetragonal blocks. In this example, prism blocks are identified according to a topological network criterion for confined ice [1]. On running the example, an output top-level directory named *runOne* is created. 

## Steps to Run the Example

In order to run this example, without making any changes to the example files, please follow the steps below.

- Download the LAMMPS trajectory file from [here on figshare](https://figshare.com/articles/Quasi-1D_Nanotube_LAMMPS_Trajectory/11448768). Copy the downloaded trajectory file, entitled *dump-240-square.lammpstrj*, into the *traj* folder inside the top-level directory *input*. Alternatively, you could change the path to the trajectory file in the *conf.yaml* file:
```{.lua}
trajectory: "path/to/trajectory/file"
```
- You can obtain the other input files required from *example_lua/iceNanotube* folder. Copy the contents of the *iceNanotube* into the top-level *lua_inputs* directory.
- You can change the frames to be analyzed by updating the options in the *vars.lua* file. The starting and ending frames are inclusive, starting from 1 onwards, irrespective of the timestep number.
- A custom volume slice can also be defined in the *vars.lua* file.
- The *functions.lua* file actually contains the Lua functions which interface with the C++ backend.

## Analyzing the Output

Inside the output directory, a directory called *topoINT* is created. Inside the *topoINT* directory, a file called *nPrisms.dat* contains the number of prism blocks and the corresponding height%, for each frame. Inside *runOne/topoINT/dataFiles*, LAMMPS data files which are numbered according to the frame number are created. These data files can be visualized in [OVITO](https://www.ovito.org/) or [VMD](http://www.ks.uiuc.edu/Research/vmd/), although OVITO is recommended for optimal type visualization.

## References 

1. Goswami, A., & Singh, J. K. (2020). A general topological network criterion for exploring the structure of icy nanoribbons and monolayers. Physical Chemistry Chemical Physics. doi:10.1039/c9cp04902a
2. d-SEAMS paper


