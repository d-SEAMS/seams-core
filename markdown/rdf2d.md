# In-Plane Radial Distribution Function 

The trajectory file for this example is [here on
figshare](https://figshare.com/articles/In_plane_2D_RDF_LAMMPS_Trajectory/11448711).
The trajectory file details the fMSI (flat Monolayer Square Ice) formed when quasi-one-dimensional water at 320 K. In this example, the in-plane RDF of the quasi-two-dimensional fMSI is determined.

The radial distribution function \f$g(r)\f$, or pair distribution function. It
can be used to illuminate features of short-range and long-range order. The RDF
is the probability of finding a particle at a distance of \f$r\f$ from a tagged
reference particle, relative to that of an ideal gas. For a system of \f$N\f$
particles, the pair correlation function for \f$N(N-1)\f$ pairs is:

  \f[
  \rho_N^{(2)}(r,r') = ⟨\sum_{i=1}^{N} \sum_{j=1,j \neq i}^{N} \delta
 (r-r_i) \delta (r'-r_j) \langle \f]

The C++ backend code essentially bins distances between pairs of particles, and normalizes the resulting histogram. The normalization is done with respect to an ideal gas.

On running the example, an output top-level directory named _runOne_ is created. 

## Steps to Run the Example

In order to run this example, without making any changes to the example files, please follow the steps below.

- Download the LAMMPS trajectory file from [here on
  figshare](https://figshare.com/articles/In_plane_2D_RDF_LAMMPS_Trajectory/11448711).
  Copy the downloaded trajectory file, entitled _dump-320.lammpstrj_,
  into the _traj_ folder inside the top-level directory _input_. Alternatively,
  you could change the path to the trajectory file in the _conf.yaml_ file:

```{.lua}
trajectory: "path/to/trajectory/file"
```

- You can obtain the other input files required from _example_lua/rdf2D-example_
  folder. Copy the contents of the _rdf2D-example_ into the top-level _lua_inputs_
  directory.
- You can change the frames to be analyzed by updating the options in the
  _vars.lua_ file. The starting and ending frames are inclusive, starting from 1
  onwards, irrespective of the timestep number.
- A custom volume slice can also be defined in the _vars.lua_ file.
- The _functions.lua_ file actually contains the Lua functions which interface
  with the C++ backend.

## Analyzing the Output

Inside the output directory, a directory called _topoMonolayer_ is created. Inside the
_topoMonolayer_ directory, files called _coverageAreaXY.dat_, _coverageAreaXZ.dat_, and _coverageAreaYZ.dat_ contain the number of rings and and the corresponding coverage area% [1], for each frame. Here, the confining sheet is in the XY plane, so the _coverageAreaXY.dat_ contains the coverage area% and quantities of interest. Inside
_runOne/topoMonolayer/dataFiles_, LAMMPS data files which are numbered according to
the frame number are created. These data files can be visualized in
[OVITO](https://www.ovito.org/) or [VMD](http://www.ks.uiuc.edu/Research/vmd/),
although OVITO is recommended for optimal type visualization.

The output file _rdf.dat_ contains the values of $r$ and $g(r)$ in each column, and can be plotted to obtain the in-plane radial distribution function. For a better and less grainy plot, you should use more frames for the RDF.

## References

1. Goswami, A., & Singh, J. K. (2020). A general topological network criterion for exploring the structure of icy nanoribbons and monolayers. Physical Chemistry Chemical Physics. doi:10.1039/c9cp04902a
