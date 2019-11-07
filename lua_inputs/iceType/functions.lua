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
    topoTwoDimDir = outDir .. "topoMonolayer/";
    lfs.mkdir(topoTwoDimDir);
    topoTwoDimData = outDir .. "topoMonolayer/dataFiles";
    lfs.mkdir(topoTwoDimData);
    --- Printing out cages per frame
    if not not printCages then
      for frame=targetFrame,finalFrame, frameGap
       do
        lfs.mkdir(topoTwoDimDir .. frame);
      end --- end of for loop through frames
    end --- End of making folders for every frame
    ---
  end --- end of topo two dimensional dir
  -- Bulk topological network criterion 
  if not not topoBulk then
    topoBulkDir = outDir .. "bulkTopo";
    lfs.mkdir(topoBulkDir);
    topoBulkDataDir = outDir .. "bulkTopo/dataFiles";
    lfs.mkdir(topoBulkDataDir);
    -- Create file for cageData (no. of cages, rings etc.)
    prismFileName = outDir .. "bulkTopo/cageData.dat";
    prismFile=io.open(prismFileName, "w"); --- Allow overwriting (otherwise use a)
    io.output(prismFile);
    --- appends a word test to the last line of the file
    io.write("Frame HCnumber DDCnumber MixedRingNumber PrismaticRings basalRings\n");
    --- closes the open file
    io.close(prismFile);
  end --- end of topo two dimensional dir
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

--- Make the directories and files
make_output_dirs( doBOP, topoOneDim, topoTwoDim, topoBulk );
clusterStatsFile(); --- For cluster stats file

slice={0,0,0}; --- This is not in use
for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrameOnlyOne(trajectory,frame,resCloud,oxygenAtomType,false,slice,slice) --- Get the frame
   nList=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist by ID
   ---
   --- Since the bulk topological network criteria are slow for dense systems,
   --- it is recommended to apply it on the largest ice cluster.
   --- If you want to apply it on all the particles in the box, uncomment the following lines
   --- and use resCloud as the pointCloud for the ring analyses. 
   clusterAnalysis(outDir, clusterCloud, resCloud, nList, iceNeighbourList, cutoffRadius, "q6");
   --- Recenter the cluster such that the centroid is at the center of the simulation box 
   recenterCluster(clusterCloud, iceNeighbourList);
   --- End of getting the largest ice cluster
   ---
   --- Start of analysis using rings (by index from here onwards.)
   rings=getPrimitiveRings(iceNeighbourList,maxDepth); --- Gets every ring (non-primitives included)
   bulkTopologicalNetworkCriterion(outDir, rings, iceNeighbourList, clusterCloud, printCages); --- Finds DDCs and HCs
   --- Clean up for ice cluster and neighbour
end
print("\nFinito\n");
