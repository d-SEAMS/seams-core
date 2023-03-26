print("\n Welcome to the manual lua function evaluation environment.\n");

for frame=targetFrame,finalFrame,frameGap do
   resCloud=readFrameOnlyOne(trajectory,frame,resCloud,oxygenAtomType,isSlice,sliceLowerLimits,sliceUpperLimits) --- Get the frame
   nList=neighborList(cutoffRadius, resCloud, oxygenAtomType); --- Calculate the neighborlist by ID
   hbnList=getHbondNetwork(trajectory,resCloud,nList,frame,hydrogenAtomType) --- Get the hydrogen-bonded network for the current frame
   hbnList=bondNetworkByIndex(resCloud,hbnList) --- Hydrogen-bonded network using indices not IDs
   rings=getPrimitiveRings(hbnList,maxDepth); --- Gets every ring (non-primitives included)
   bulkRingNumberAnalysis(outDir, rings, hbnList, resCloud, maxDepth, targetFrame); --- Writes out primitive rings for a bulk system 
end
