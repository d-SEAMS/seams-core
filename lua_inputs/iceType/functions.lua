print("\n Welcome to the manual lua function evaluation environment.\n");

--- Variables again
dumpName="water.lammpstrj";
--- Prep the files (O==one, T==two)
tmpFileO=io.open(chillPlus_noMod, "w"); --- Allow overwriting (otherwise use a)
--- sets the default output file as test.lua
io.output(tmpFileO);
--- appends a word test to the last line of the file
io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n")
--- closes the open file
io.close(tmpFileO)
--- Do it again
tmpFileT=io.open(chillPlus_mod, "w"); --- Allow overwriting (otherwise use a)
io.output(tmpFileT);
io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n");
io.close(tmpFileT);

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
end
print("\nWORM_PWNED\n");
