print("\n Welcome to the Ice Type Determination Module\n");
cutoffRadius = 3.5; --- This is for H2O
oxygenAtomType = 2; --- This is assigned by LAMMPS
hydrogenAtomType = 1; --- Hydrogen atom type assigned
targetFrame=1; --- The first frame
finalFrame=1; --- This is inclusive
frameGap=1; --- The gap between frames
maxDepth = 4; --- The maximum depth upto which rings will be searched.
--- Slice Information
isSlice = true; --- This is true if the analysis is to be done only for a volume slice
sliceLowerLimits = {0,0,0}; --- Lower limit of the slice (for box dim, keep the values the same as 0)
sliceUpperLimits = {50,0,0}; --- Upper limit of the slice   

--- Paths for the output directories and lua scipt
outDir="runOne/"; --- The subdirectory used; 
functionScript="lua_inputs/iceType/functions.lua" --- This is relative to the binary location

--- Variables for the monolayer only:
confiningSheetArea = 50*50; 

--- Variables for the RDF only
rdf = true; --- This should only be set to true if you want to calculate the RDF
rdfCutoff = 12; --- This should be less than half the box length
binwidth = 0.05; --- This is the binwidth or delta_r