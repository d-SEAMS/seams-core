# Bulk Ice Classification Using Topological Network Criteria

The trajectory file for this example is [here on figshare](https://figshare.com/articles/Nucleation_LAMMPS_Trajectory/11448702). The trajectory file details a portion of a long run of 4096 molecules of mW water, which have undergone crystallization. In this example, DDCs (Double-diamond Cages), HCs (Hexagonal Cages), and mixed rings are identified inside the largest ice cluster. On running the example, an output top-level directory named *runOne* is created. 

## Steps to Run the Example

In order to run this example, without making any changes to the example files, please follow the steps below.

- Download the LAMMPS trajectory file from [here on figshare](https://figshare.com/articles/Nucleation_LAMMPS_Trajectory/11448702). Copy the downloaded trajectory file, entitled *nucleation.lammpstrj*, into the *traj* folder inside the top-level directory *input*. Alternatively, you could change the path to the trajectory file in the *conf.yaml* file:
```{.lua}
trajectory: "path/to/trajectory/file"
```
- You can obtain the other input files required from *example_lua/bulkTopologicalCriterion* folder. Copy the contents of the *bulkTopologicalCriterion* into the top-level *lua_inputs* directory.
- You can change the frames to be analyzed by updating the options in the *vars.lua* file. The starting and ending frames are inclusive, starting from 1 onwards, irrespective of the timestep number.
- A custom volume slice can also be defined in the *vars.lua* file.
- The *functions.lua* file actually contains the Lua functions which interface with the C++ backend.

## Analyzing the Output

Inside the output directory, a file called *clusterStats.dat* contains the cluster statistics for each frame. Inside the *bulkTopo* directory, *cageData.dat* contains the number of HCs, DDCs and mixed rings for each frame. Inside *runOne/bulkTopo/dataFiles*, LAMMPS data files which are numbered according to the frame number are created. These data files can be visualized in [OVITO](https://www.ovito.org/) or [VMD](http://www.ks.uiuc.edu/Research/vmd/), although OVITO is recommended for optimal type visualization.

## References 

1. Haji-Akbari, A., & Debenedetti, P. G. (2015). Direct calculation of ice homogeneous nucleation rate for a molecular model of water. Proceedings of the National Academy of Sciences, 112(34), 10582-10588. doi:10.1073/pnas.1509267112


