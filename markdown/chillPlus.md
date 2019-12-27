# Bulk Ice Classification Using the CHILL+ Criterion

The trajectory file for this example is [here on figshare](https://figshare.com/articles/CHILL_LAMMPS_Trajectory/11448720). The trajectory file details a portion of a long run of 4096 molecules of a perfect Ic (cubic ice) lattice. In this example, the all 4096 molecules are identified as Ic, according to the CHILL+ algorithm [1]. On running the example, an output top-level directory named *runOne* is created. 

## Steps to Run the Example

In order to run this example, without making any changes to the example files, please follow the steps below.

- Download the LAMMPS trajectory file from [here on figshare](https://figshare.com/articles/CHILL_LAMMPS_Trajectory/11448720). Copy the downloaded trajectory file, entitled *nucleation.lammpstrj*, into the *traj* folder inside the top-level directory *input*. Alternatively, you could change the path to the trajectory file in the *conf.yaml* file:
```{.lua}
trajectory: "path/to/trajectory/file"
```
- You can obtain the other input files required from *example_lua/chillPlus* folder. Copy the contents of the *chillPlus* into the top-level *lua_inputs* directory.
- You can change the frames to be analyzed by updating the options in the *vars.lua* file. The starting and ending frames are inclusive, starting from 1 onwards, irrespective of the timestep number.
- A custom volume slice can also be defined in the *vars.lua* file.
- The *functions.lua* file actually contains the Lua functions which interface with the C++ backend.

## Analyzing the Output

Inside the output directory, a file called *clusterStats.dat* contains the cluster statistics for each frame. Each frame is also written out to a LAMMPS trajectory file called *waterChillP.lammpstrj*, containing the classified atoms, and the *largestIce.lammpstrj*, which has the particles in the largest ice cluster. These trajectory files can be visualized in [OVITO](https://www.ovito.org/) or [VMD](http://www.ks.uiuc.edu/Research/vmd/). Inside the *bop* directory, *chillPlus.txt* contains the number of particles identified as cubic ice (Ic), hexagonal ice (Ih), interfacial ice (Interfacial), clathrates (Clath), interfacial clathrates (InterClath), water and the total number of molecules (Total). 

## References 

1. Nguyen, A. H., & Molinero, V. (2014). Identification of Clathrate Hydrates, Hexagonal Ice, Cubic Ice, and Liquid Water in Simulations: the CHILL+ Algorithm. The Journal of Physical Chemistry B, 119(29), 9369-9376. doi:10.1021/jp510289t




