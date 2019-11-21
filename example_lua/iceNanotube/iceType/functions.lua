print("\n Welcome to the manual lua function evaluation environment.\n");

--- Init Modules
local lfs = require"lfs"

--- Call functions defined in script file
package.path = './../lua_inputs/luaFunctions/?.lua;' .. package.path
local luaFunctions = require("scripts");

--- Make the directories
make_output_dirs( doBOP, topoOneDim, topoTwoDim, topoBulk );

for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrameOnlyOne(trajectory,frame,resCloud,oxygenAtomType,isSlice,sliceLowerLimits,sliceUpperLimits) --- Get the frame
   nList=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist by ID
   hbnList=getHbondNetwork(trajectory,resCloud,nList,frame,hydrogenAtomType) --- Get the hydrogen-bonded network for the current frame
   hbnList=bondNetworkByIndex(resCloud,hbnList) --- Hydrogen-bonded network using indices not IDs
   rings=getPrimitiveRings(hbnList,maxDepth); --- Gets every ring (non-primitives included)
   prismAnalysis(outDir, rings, hbnList, resCloud, maxDepth); --- Does the prism analysis for quasi-one-dimensional ice
end
