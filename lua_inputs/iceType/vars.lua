print("\n Welcome to the Ice Type Determination Module\n");
cutoffRadius = 3.5; --- This is for H2O
oxygenAtomType = 1; --- This is assigned by LAMMPS
targetFrame=1; --- The first frame
finalFrame=1; --- This is inclusive
frameGap=1; --- The gap between frames
dumpName="wat.lammpstrj"; --- Output file name

--- DO NOT ENABLE if you have not read the instructions
defineFunctions=false; --- The last test before all hell breaks loose
functionScript="../lua_inputs/iceType/functions.lua" --- This is relative to the binary location
