print("\n Welcome to the Ice Type Determination Module\n");
cutoffRadius = 3.5; --- This is for H2O
oxygenAtomType = 2; --- This is assigned by LAMMPS
hydrogenAtomType = 1; --- Hydrogen atom type assigned
targetFrame=100; --- The first frame
finalFrame=100; --- This is inclusive
frameGap=1; --- The gap between frames
maxDepth = 7; --- The maximum depth upto which rings will be searched. 

--- DO NOT ENABLE if you have not read the instructions
defineFunctions=true; --- The last test before all hell breaks loose
outDir="../runOne/"; --- The subdirectory used; 
functionScript="../lua_inputs/iceType/functions.lua" --- This is relative to the binary location

--- Variables for creating directories, maybe shift to the config later?
doBOP = false; --- Bond orienational parameter analysis (yes/no)
topoOneDim = true; --- Topological network criterion for one-dimensional prisms
topoTwoDim = false; --- Topological network criterion for two-dimensional ice
topoBulk = false; --- DDC/HC criterion for bulk ice 