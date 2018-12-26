print("\n Welcome to the manual lua function evaluation environment.\n");

dumpName="water.lammpstrj";
slice={0,0,0}; --- This is not in use
for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrame(trajectory, frame, resCloud, oxygenAtomType,false,x,x) --- Get the frame
   resCloud=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist
   resCloud=chillPlus_cij(resCloud,false); --- Calculate Cij (cloud,slice)
   resCloud=chillPlus_iceType(resCloud,false,"Something"); --- Write out data (cloud,slice,name)
   avgQ6=averageQ6(resCloud,false); --- Average Q6 (cloud,slice)
   resCloud=modifyChill(resCloud,avgQ6); --- Modification (cloud,q6)
   percentage_Ice(resCloud,false,"SomeOther"); --- Post reclassification writeOut
   writeDump(resCloud,dumpName);
   writeHistogram(resCloud,avgQ6);
print("\nWORM\n");
end
print("\nWORM_PWNED\n");
