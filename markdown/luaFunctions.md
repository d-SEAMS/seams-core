# Lua Function Documentation

In the _functions.lua_ file, lua functions are called, which are registered on the C++ side to inferface with the C++ functions. Here, we document the lua functions currently available to the user. 

## Currently Registered Lua Functions

The workflows for quasi-two-dimensional ice, quasi-one-dimensional ice and bulk systems are separated. The `Lua` functions for each work-flow are registered in different blocks in the `C++` code.

### Common Functions

The following `Lua` functions interface to the same `C++` functions in every work-flow code block.

- **readFrameOnlyOne**: interfaces with \ref sinp::readLammpsTrjreduced <br>The input variables are:
  + *trajectory* - The name of the lammps trajectory file to be read in, obtained from the _config.yml_ file.
  + *frame* - The current frame number whose information will be read in.
  + *resCloud* - The outputted \ref molSys::PointCloud, which has been passed to the `Lua` side. 
  + *atomType* - The type ID of the desired type of atoms. Defined in _vars.lua_ in the examples.
  + *isSlice* - This decides whether a slice will be created (true) or not (false). Defined in _vars.lua_.
  + *sliceLowerLimits* - Contains the lower limits of the slice. If *isSlice* false, this will be ignored.
  + *sliceUpperLimits* - Contains the upper limits of the slice. If *isSlice* false, this will be ignored. <br> This function returns the *resCloud* \ref molSys::PointCloud.

- **neighborList**: interfaces with \ref nneigh::neighListO <br>The input variables are:
  + *cutoffRadius* - The name of the lammps trajectory file to be read in, obtained from the _config.yml_ file.
  + *resCloud* - The input \ref molSys::PointCloud, which has been passed to the `Lua` side. 
  + *atomType* - The type ID of the desired type of atoms. Defined in _vars.lua_ in the examples. <br> This function returns the *nList* vector of vectors for the neighbour list (by atom ID).

- **getHbondNetwork**: interfaces with \ref bond::populateHbonds <br>The input variables are:
  + *trajectory* - The name of the lammps trajectory file to be read in, obtained from the _config.yml_ file.
  + *resCloud* - The input \ref molSys::PointCloud, which has been passed to the `Lua` side. 
  + *nList* - Row-ordered neighbour list by atom ID.
  + *frame* - The current frame number whose information will be read in.
  + *hydrogenAtomType* - The type ID of the hydrogen atoms. Defined in _vars.lua_ in the examples. <br> This function returns the *hbnList* vector of vectors for the hydrogen-bonded neighbour list (by atom ID).
  
- **bondNetworkByIndex**: interfaces with \ref nneigh::neighbourListByIndex <br> The input variables are:
  + *resCloud* - The input \ref molSys::PointCloud, which has been passed to the `Lua` side. 
  + *nList* - Row-ordered neighbour list by atom ID. <br> This function returns the *hbnList* vector of vectors for the row-ordered full neighbour list, by index, NOT atom ID.

- **getPrimitiveRings**: interfaces with \ref primitive::ringNetwork <br>The input variables are:
  + *hbnList* - Row-ordered neighbour list by index (and NOT the atom ID).
  + *maxDepth* - The maximum depth upto which rings will be searched. This means that rings larger than maxDepth in length will not be generated. <br> This function returns the *rings* vector of vectors for the primitive rings; each ring contains the atom indices of the ring members.

### Quasi-Two-Dimensional Ice

- **ringAnalysis**: interfaces with \ref ring::polygonRingAnalysis <br>The input variables are:
  + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
  + *rings* - Row-ordered vector of vectors for rings of a single type. Registered on the `Lua` side. 
  + *nList* - Row-ordered neighbour list by index. 
  + *resCloud* - The input \ref molSys::PointCloud, which has been passed to the `Lua` side.
  + *maxDepth* - The maximum depth upto which rings will be searched. This means that rings larger than maxDepth in length will not be generated.
  + *confiningSheetArea* - Area calculated using the two significant dimensions of the quasi-two-dimensional sheet. Defined in _vars.lua_.

- **calcRDF**: interfaces with \ref rdf2::rdf2Danalysis_AA <br>The input variables are:
  + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
  + *rdf* - Vector containing the RDF values. Must be passed to `Lua`.
  + *resCloud* - The input \ref molSys::PointCloud, which has been passed to the `Lua` side. 
  + *rdfCutoff* - Cutoff for the RDF. This should not be greater than half the box length. Defined in _vars.lua_.
  + *binwidth* - Width of the bin for histogramming. Defined in _vars.lua_.
  + *targetFrame* - The first frame for RDF binning. Defined in _vars.lua_.
  + *finalFrame* - The final frame for RDF binning. Defined in _vars.lua_. <br> This function returns the *hbnList* vector of vectors for the hydrogen-bonded neighbour list (by atom ID).

### Quasi-One-Dimensional Ice
  
- **prismAnalysis**: interfaces with \ref ring::prismAnalysis <br>The input variables are:
  + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
  + *rings* - Row-ordered vector of vectors for rings of a single type. Registered on the `Lua` side. 
  + *nList* - Row-ordered neighbour list by index. 
  + *resCloud* - The input \ref molSys::PointCloud, which has been passed to the `Lua` side.
  + *maxDepth* - The maximum depth upto which rings will be searched. This means that rings larger than maxDepth in length will not be generated.
  + *confiningSheetArea* - Area calculated using the two significant dimensions of the quasi-two-dimensional sheet. Defined in _vars.lua_.

### Bulk Systems

- **readFrame**: interfaces with \ref sinp::readLammpsTrjO <br>The input variables are:
  + *trajectory* - The name of the lammps trajectory file to be read in, obtained from the _config.yml_ file.
  + *frame* - The current frame number whose information will be read in.
  + *resCloud* - The \ref molSys::PointCloud for the current frame, which has been passed to the `Lua` side. 
  + *oxygenAtomType* - The type ID of the Oxygen atoms in the LAMMPS trajectory file. Defined in _vars.lua_ in the examples.
  + *isSlice* - This decides whether a slice will be created (true) or not (false). Defined in _vars.lua_.
  + *sliceLowerLimits* - Contains the lower limits of the slice. If *isSlice* false, this will be ignored.
  + *sliceUpperLimits* - Contains the upper limits of the slice. If *isSlice* false, this will be ignored. <br> This function returns the *resCloud* \ref molSys::PointCloud.

- **chillPlus_cij**: interfaces with \ref chill::getCorrelPlus <br>This function uses the CHILL+ algorithm. <br> The input variables are:
  + *resCloud* - The \ref molSys::PointCloud, which has been passed to the `Lua` side. 
  + *nList* - Row-ordered neighbour list by atom ID. Passed to the `Lua` side.
  + *isSlice* - This decides whether a slice will be created (true) or not (false). Defined in _vars.lua_. <br> This function returns the *resCloud* \ref molSys::PointCloud.

- **chillPlus_iceType**: interfaces with \ref chill::getIceTypePlus <br> This function is for the CHILL+ algorithm. <br>The input variables are:
   + *resCloud* - The \ref molSys::PointCloud, which has been passed to the `Lua` side.
   + *nList* - Row-ordered neighbour list by atom ID. Passed to the `Lua` side.
   + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
   + *isSlice* - This decides whether a slice will be created (true) or not (false). Defined in _vars.lua_. <br> This function returns the *resCloud* \ref molSys::PointCloud.

- **chill_cij**: interfaces with \ref chill::getCorrel <br>This function uses the CHILL algorithm. <br> The input variables are:
  + *resCloud* - The outputted \ref molSys::PointCloud, which has been passed to the `Lua` side. 
  + *nList* - Row-ordered neighbour list by atom ID. Passed to the `Lua` side.
  + *isSlice* - This decides whether a slice will be created (true) or not (false). Defined in _vars.lua_. <br> This function returns the *resCloud* \ref molSys::PointCloud.

- **chill_iceType**: interfaces with \ref chill::getIceType <br> This function is for the CHILL algorithm. <br>The input variables are:
   + *resCloud* - The outputted \ref molSys::PointCloud, which has been passed to the `Lua` side.
   + *nList* - Row-ordered neighbour list by atom ID. Passed to the `Lua` side.
   + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
   + *isSlice* - This decides whether a slice will be created (true) or not (false). Defined in _vars.lua_. <br> This function returns the *resCloud* \ref molSys::PointCloud.

- **percentage_Ice**: interfaces with \ref chill::printIceType <br> This function is for the CHILL algorithm. <br>The input variables are:
   + *resCloud* - The \ref molSys::PointCloud, which has been passed to the `Lua` side.
   + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
   + *isSlice* - This decides whether a slice will be created (true) or not (false). Defined in _vars.lua_. 
   + *outputFileName* - File name of the output file, to which the per-particle ice types will be written out. Defined in _vars.lua_. <br> This function returns the *resCloud* \ref molSys::PointCloud.

- **clusterAnalysis**: interfaces with \ref clump::clusterAnalysis <br> The input variables are:
   + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
   + *iceCloud* - The \ref molSys::PointCloud for the desired largest ice cluster, which has been passed to the `Lua` side.
   + *resCloud* - The \ref molSys::PointCloud for all the particles read in, for the current frame, which has been passed to the `Lua` side. 
   + *nList* - Row-ordered neighbour list by atom ID. Passed to the `Lua` side.
   + *iceNeighbourList* - Row-ordered neighbour list by atom index, not ID, according to the iceCloud atoms. Defined in _vars.lua_. 
   + *cutoff* - Cutoff for the nearest neighbours. Passed to the `Lua` side.
   + *bopAnalysis* - This determines which method to use for determining the ice-like nature of the particles. This can be "q6" or "chill", for using the \f$Q_6\f$ parameter or CHILL algorithm, respectively. Defined in _vars.lua_. <br> This function returns the *iceCloud* \ref molSys::PointCloud for the largest ice cluster.

- **recenterCluster**: interfaces with \ref clump::recenterClusterCloud <br> The input variables are:
   + *iceCloud* - The \ref molSys::PointCloud for particles in the largest ice cluster, which has been passed to the `Lua` side.
   + *iceNeighbourList* - Row-ordered neighbour list by atom index, not ID, according to the iceCloud atoms. Defined in _vars.lua_.

- **bulkTopologicalNetworkCriterion**: interfaces with \ref ring::topoBulkAnalysis <br>The input variables are:
  + *outDir* - The string to the output directory, in which files will be written out. Defined in _vars.lua_.
  + *rings* - Row-ordered vector of vectors for rings of a single type. Registered on the `Lua` side. 
  + *nList* - Row-ordered neighbour list by index. Passed to the `Lua` side. 
  + *resCloud* - The input \ref molSys::PointCloud, which has been passed to the `Lua` side.
  + *printCages* - Flag for printing the information of each cage in the frame (true) or not printing the coordinates/connectivity of each cage (false). Defined in _vars.lua_. 

## Extending d-SEAMS 

d-SEAMS can be extended by writing `C++` functions in pre-existing header files or by creating header files of your own. 

### Registering Lua Functions

In order to call a `C++` function from the `Lua` script, the corresponding `Lua` function must be registered. This can be done simply in the _main.cpp_ file:

```{.cpp}
lua.set_function("luaFuncName", cppFuncName);
```
where _luaFuncName_ is the name of the lua function and _cppFuncName_ is the name of the corresponding C++ function. It should be noted that default arguments must be provided explicitly when the _luaFuncName_ function is called on the `Lua` side, or else the program will break with an error. 

### Passing Variables to Lua

If a `C++` struct or variable is to be defined or passed to the `Lua` side, then the variable must be defined using the following syntax inside _main.cpp_ :

```{.cpp}
lua["luaStruct"] = &cppStruct;
lua["luaVecOfVectors"] = &cppVecOfVectors;
lua["luaString"] = myString;
```
`C++` structs, vectors of vectors and other objects must be passed by reference. This need not be done for strings. 

