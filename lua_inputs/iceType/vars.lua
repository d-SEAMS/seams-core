print("\n Welcome to the Ice Type Determination Module\n");
cutoffRadius = 3.5; --- This is for H2O
oxygenAtomType = 1; --- This is assigned by LAMMPS
targetFrame=1; --- The first frame
finalFrame=10; --- This is inclusive
frameGap=1; --- The gap between frames
dumpName="wat.lammpstrj"; --- Output file name
chillPlus_mod="supaaChill.txt"; --- This the modified file
chillPlus_noMod="chillPlus.txt"; --- This is the standard file name
chill_noMod="chill.txt"; --- This is the standard file name
largest_ice_cluster_name="largeCluster.txt"; --- This is the standard file name

--- DO NOT ENABLE if you have not read the instructions
defineFunctions=true; --- The last test before all hell breaks loose
functionScript="../lua_inputs/iceType/functions.lua" --- This is relative to the binary location
