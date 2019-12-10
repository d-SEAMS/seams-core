print("\n Welcome to the manual lua function evaluation environment.\n");

--- Init Modules
local lfs = require"lfs"

--- Call functions defined in script file
package.path = './../lua_inputs/luaModules/?.lua;' .. package.path
local luaFunctions = require("scripts");

--- Make the directories
luaFunctions.make_output_dirs( doBOP, topoOneDim, topoTwoDim, topoBulk );
luaFunctions.clusterStatsFile(); --- For cluster stats file

for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrameOnlyOne(trajectory,frame,resCloud,oxygenAtomType,isSlice,sliceLowerLimits,sliceUpperLimits) --- Get the frame
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
end
