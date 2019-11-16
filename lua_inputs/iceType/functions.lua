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
  end --- end of topo two dimensional dir
end
---

--- Make the directories
make_output_dirs( doBOP, topoOneDim, topoTwoDim, topoBulk );

for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrameOnlyOne(trajectory,frame,resCloud,oxygenAtomType,isSlice,sliceLowerLimits,sliceUpperLimits) --- Get the frame
   nList=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist by ID
   hbnList=getHbondNetwork(trajectory,resCloud,nList,frame,hydrogenAtomType) --- Get the hydrogen-bonded network for the current frame
   hbnList=bondNetworkByIndex(resCloud,hbnList) --- Hydrogen-bonded network using indices not IDs
   rings=getPrimitiveRings(hbnList,maxDepth); --- Gets every ring (non-primitives included)
   ringAnalysis(outDir, rings, hbnList, resCloud, maxDepth, confiningSheetArea); --- Does the ring analysis for quasi-two-dimensional ice
end
print("\nFinito\n");
