print("\n Welcome to the manual lua function evaluation environment.\n");

--- Init Modules
local lfs = require"lfs"

--- Strings are immutable, check https://stackoverflow.com/questions/1405583/concatenation-of-strings-in-lua
--- Handle this better
lfs.chdir("../output/");
lfs.mkdir(subdir);

--- Variables again (with subdir)
dumpName= subdir .. "water.lammpstrj";
cpMod= subdir .. chillPlus_mod;
cpNMod= subdir .. chillPlus_noMod;
licn= subdir .. largest_ice_cluster_name;

--- Prep the files (O==one, T==two)
tmpFileO=io.open(cpNMod, "w"); --- Allow overwriting (otherwise use a)
--- sets the default output file as test.lua
io.output(tmpFileO);
--- appends a word test to the last line of the file
io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n")
--- closes the open file
io.close(tmpFileO)
--- Do it again
tmpFileT=io.open(cpMod, "w"); --- Allow overwriting (otherwise use a)
io.output(tmpFileT);
io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n");
io.close(tmpFileT);
--- Once more for the cluster
tmpFileC=io.open(licn, "w"); --- Allow overwriting (otherwise use a)
io.output(tmpFileC);
io.write("Frame number_in_cluster\n");
io.close(tmpFileC);

slice={0,0,0}; --- This is not in use
for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrame(trajectory, frame, resCloud, oxygenAtomType,false,slice,slice) --- Get the frame
   resCloud=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist
   resCloud=chillPlus_cij(resCloud,false); --- Calculate Cij (cloud,slice)
   resCloud=chillPlus_iceType(resCloud,false,chillPlus_noMod); --- Write out data (cloud,slice,name)
   avgQ6=averageQ6(resCloud,false); --- Average Q6 (cloud,slice)
   resCloud=modifyChill(resCloud,avgQ6); --- Modification (cloud,q6)
   percentage_Ice(resCloud,false,chillPlus_mod); --- Post reclassification writeOut
   writeDump(resCloud,dumpName);
   writeHistogram(resCloud,avgQ6);
   --- Do the largest Ice cluster stuff
   clusterCloud=create_cluster(resCloud,clusterCloud);
   largest_ice_cluster=largest_cluster(clusterCloud,cutoffRadius,true,false);
   writeCluster(clusterCloud,largest_ice_cluster_name,false,largest_ice_cluster);
end
print("\nWORM_PWNED\n");
