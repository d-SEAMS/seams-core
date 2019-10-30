print("\n Welcome to the manual lua function evaluation environment.\n");

-- --- Init Modules
-- local lfs = require"lfs"

-- --- Strings are immutable, check https://stackoverflow.com/questions/1405583/concatenation-of-strings-in-lua
-- lfs.mkdir(outDir);
-- lfs.chdir(outDir);
-- os.execute("mkdir " .. outDir)

-- --- Variables again (with subdir)
-- dumpChillP= "waterChillP.lammpstrj";
-- dumpSupaaP= "waterSupaaP.lammpstrj";

-- --- Prep the files (O==one, T==two)
-- tmpFileO=io.open(chillPlus_noMod, "w"); --- Allow overwriting (otherwise use a)
-- --- sets the default output file as test.lua
-- io.output(tmpFileO);
-- --- appends a word test to the last line of the file
-- io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n")
-- --- closes the open file
-- io.close(tmpFileO)
-- --- Do it again
-- tmpFileT=io.open(chillPlus_mod, "w"); --- Allow overwriting (otherwise use a)
-- io.output(tmpFileT);
-- io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n");
-- io.close(tmpFileT);
-- --- Once more for the cluster
-- tmpFileC=io.open(largest_ice_cluster_name, "w"); --- Allow overwriting (otherwise use a)
-- io.output(tmpFileC);
-- io.write("Frame number_in_cluster\n");
-- io.close(tmpFileC);

slice={0,0,0}; --- This is not in use
for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrameOnlyOne(trajectory,frame,resCloud,oxygenAtomType,false,slice,slice) --- Get the frame
   nList=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist
   hbnList=getHbondNetwork(trajectory,resCloud,nList,frame,hydrogenAtomType) --- Get the hydrogen-bonded network for the current frame
   -- graph=countEveryRing(resCloud, nList,7); --- Gets every ring (non-primitives included)
end
print("\nFinito\n");
