print("\n Welcome to the Ice Type Determination Module\n");
cutoffRadius = 3.5; --- This is for H2O
oxygenAtomType = 1; --- This is assigned by LAMMPS
hydrogenAtomType = 1; --- Hydrogen atom type assigned (not used here)
targetFrame=1; --- The first frame
finalFrame=1; --- This is inclusive
frameGap=1; --- The gap between frames
maxDepth = 6; --- The maximum depth upto which rings will be searched.
--- Slice Information
isSlice = false; --- This is true if the analysis is to be done only for a volume slice
sliceLowerLimits = {0,0,0}; --- Lower limit of the slice (for box dim, keep the values the same as 0)
sliceUpperLimits = {0,0,0}; --- Upper limit of the slice   

--- Paths for the output directories and lua scipt
outDir="runOne/"; --- The subdirectory used; Keep the slash at the end 
functionScript="lua_inputs/iceType/functions.lua" --- This is relative to the binary location

-- Variable for the topological network criterion
printCages = false; --- Prints out every cage for every frame if true

-- File names for outputs 
dumpName="wat.lammpstrj"; --- Output file name
chillPlus_mod="supaaChill.txt"; --- This the modified file
chillPlus_noMod="chillPlus.txt"; --- This is the standard file name
chill_noMod="chill.txt"; --- This is the standard file name
largest_ice_cluster_name="largeCluster.txt"; --- This is the standard file name 
dumpChillP= "waterChillP.lammpstrj"; --- Output dump file for CHILL+ classification
dumpSupaaP= "waterSupaaP.lammpstrj"; --- Output dump file for the SUPER CHILL+ classification
largestClusterDump = "largestIce.lammpstrj"; --- Output dump file for the largest ice cluster