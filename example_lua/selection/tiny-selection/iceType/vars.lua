print("\n Welcome to the Ice Type Determination Module\n");
cutoffRadius = 3.5; --- This is for H2O
oxygenAtomType = 1; --- This is assigned by LAMMPS
hydrogenAtomType = 2; --- Hydrogen atom type assigned
targetFrame=1; --- The first frame
finalFrame=1; --- This is inclusive
frameGap=1; --- The gap between frames
maxDepth = 6; --- The maximum depth upto which rings will be searched.
--- Slice Information
isSlice = true; --- This is true if the analysis is to be done only for a volume slice
sliceLowerLimits = {0,0,0}; --- Lower limit of the slice (for box dim, keep the values the same as 0)
sliceUpperLimits = {31.728,0,0}; --- Upper limit of the slice   

--- Paths for the output directories and lua scipt
outDir="runOne/"; --- The subdirectory used; 
functionScript="lua_inputs/iceType/functions.lua" --- This is relative to the binary location
