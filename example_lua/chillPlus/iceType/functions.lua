print("\n Welcome to the manual lua function evaluation environment.\n");

--- Init Modules
local lfs = require"lfs"

--- Functions for creating directories
function make_output_dirs( doBOP, topoOneDim, topoTwoDim, topoBulk )
  -- Make the main output folder:
  lfs.mkdir(outDir);
  -- Bond orientational parameter folder
  if not not doBOP then
    bopDir = outDir .. "bop";
    lfs.mkdir(bopDir);
    --- Bulk BOP files:
    --- Prep the files (O==one, T==two)
    tmpFileO=io.open(outDir .. "bop/"..chillPlus_noMod, "w"); --- Allow overwriting (otherwise use a)
    --- sets the default output file as test.lua
    io.output(tmpFileO);
    --- appends a word test to the last line of the file
    io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n")
    --- closes the open file
    io.close(tmpFileO)
    --- Do it again
    tmpFileT=io.open(outDir .. "bop/"..chillPlus_mod, "w"); --- Allow overwriting (otherwise use a)
    io.output(tmpFileT);
    io.write("Frame Ic Ih Interfacial Clath InterClath Water Total\n");
    io.close(tmpFileT);
    ---
  end --- end of bop dir
  -- Topological network criterion for quasi-one-dimensional INTs
  if not not topoOneDim then
    topoOneDimDir = outDir .. "topoINT";
    lfs.mkdir(topoOneDimDir);
    topoOneDimDir = outDir .. "topoINT/prisms";
    lfs.mkdir(topoOneDimDir);
    topoOneDimData = outDir .. "topoINT/dataFiles";
    lfs.mkdir(topoOneDimData);
    -- Create file for nPrisms (no. of prisms)
    prismFileName = outDir .. "topoINT/nPrisms.dat";
    prismFile=io.open(prismFileName, "w"); --- Allow overwriting (otherwise use a)
    io.output(prismFile);
    --- appends a word test to the last line of the file
    io.write("Frame RingSize Num_of_prisms Height% RingSize ... Height%\n");
    --- closes the open file
    io.close(prismFile);
  end --- end of topo one dimensional dir
  -- Topological network criterion for quasi-two-dimensional INTs
  if not not topoTwoDim then
    topoTwoDimDir = outDir .. "topoMonolayer";
    lfs.mkdir(topoTwoDimDir);
    topoTwoDimData = outDir .. "topoMonolayer/dataFiles";
    lfs.mkdir(topoTwoDimData);
    -- Create file for coverageAreaXY.dat 
    areaFileName = outDir .. "topoMonolayer/coverageAreaXY.dat";
    areaFile=io.open(areaFileName, "w"); --- Allow overwriting (otherwise use a)
    io.output(areaFile);
    --- appends a word test to the last line of the file
    io.write("Frame RingSize Num_of_rings CoverageAreaXY% RingSize ... CoverageAreaXY%\n");
    --- closes the open file
    io.close(areaFile);
    -- Create file for coverageAreaXZ.dat 
    areaFileNameXZ = outDir .. "topoMonolayer/coverageAreaXZ.dat";
    areaFileXZ=io.open(areaFileNameXZ, "w"); --- Allow overwriting (otherwise use a)
    io.output(areaFileXZ);
    --- appends a word test to the last line of the file
    io.write("Frame RingSize Num_of_rings CoverageAreaXZ% RingSize ... CoverageAreaXZ%\n");
    --- closes the open file
    io.close(areaFileXZ);
    -- Create file for coverageAreaYZ.dat 
    areaFileNameYZ = outDir .. "topoMonolayer/coverageAreaYZ.dat";
    areaFileYZ=io.open(areaFileNameYZ, "w"); --- Allow overwriting (otherwise use a)
    io.output(areaFileYZ);
    --- appends a word test to the last line of the file
    io.write("Frame RingSize Num_of_rings CoverageAreaYZ% RingSize ... CoverageAreaYZ%\n");
    --- closes the open file
    io.close(areaFileYZ);
  end --- end of topo two dimensional dir
  -- Bulk topological network criterion 
  if not not topoBulk then
    topoBulkDir = outDir .. "bulkTopo";
    lfs.mkdir(topoBulkDir);
    topoBulkDataDir = outDir .. "bulkTopo/dataFiles";
    lfs.mkdir(topoBulkDataDir);
    -- Create file for cageData (no. of cages, rings etc.)
    topoFileName = outDir .. "bulkTopo/cageData.dat";
    topoFile=io.open(prismFileName, "w"); --- Allow overwriting (otherwise use a)
    io.output(topoFile);
    --- appends a word test to the last line of the file
    io.write("Frame HCnumber DDCnumber MixedRingNumber PrismaticRings basalRings\n");
    --- closes the open file
    io.close(topoFile);
  end --- end of topo bulk dir creation
end
---
function clusterStatsFile()
  -- Create the file
  clusterFileName = outDir .. "clusterStats.dat";
  clusterFile=io.open(clusterFileName, "w"); --- Allow overwriting (otherwise use a)
    io.output(clusterFile);
    --- appends a word test to the last line of the file
    io.write("Frame largestCluster numOfClusters smallestCluster avgClusterSize\n");
    --- closes the open file
    io.close(clusterFile);
end

--- Make the directories
make_output_dirs( doBOP, topoOneDim, topoTwoDim, topoBulk );
clusterStatsFile();

for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrame(trajectory,frame,resCloud,oxygenAtomType,isSlice,sliceLowerLimits,sliceUpperLimits) --- Get the frame
   nList=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist by ID
   ---
   resCloud=chillPlus_cij(resCloud,nList,isSlice); --- Calculate Cij (cloud,slice)
   resCloud=chillPlus_iceType(resCloud,nList,outDir,isSlice,chillPlus_noMod); --- Write out data (cloud,slice,name)
   writeDump(resCloud,outDir,dumpChillP); --- Dump the rescloud which currently has CHILL Plus classifications
   --- Modified CHILL+
   avgQ6=averageQ6(resCloud,nList,isSlice); --- Average Q6 (cloud,slice)
   resCloud=modifyChill(resCloud,avgQ6); --- Modification (cloud,q6)
   percentage_Ice(resCloud,outDir,isSlice,chillPlus_mod); --- Post reclassification writeOut
   writeDump(resCloud,outDir,dumpSupaaP); --- Dump the rescloud which now has the supaa CHILL Plus Trajectory
   --- Get the largest ice cluster. Here, iceNeighbourList is the neighbour list by index.
   clusterAnalysis(outDir, clusterCloud, resCloud, nList, iceNeighbourList, cutoffRadius, "q6");
   --- Recenter the cluster such that the centroid is at the center of the simulation box 
   recenterCluster(clusterCloud, iceNeighbourList);
   writeDump(resCloud,outDir,largestClusterDump); --- Dump the recentered largest ice cluster
end
