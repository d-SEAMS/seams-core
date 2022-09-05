print("\n Welcome to the manual lua function evaluation environment.\n");

for frame=targetFrame,finalFrame,frameGap do
    -- Read in all atom types without applying the slice 
    resCloud=readFrameOnlyOneAllAtoms(trajectory, frame,resCloud,false, {0,0,0}, {0,0,0});
    -- Get just the O atoms from the PointCloud with all the atoms 
    oCloud=getPointCloudAtomsOfOneAtomType(resCloud,oCloud,oxygenAtomType,isSlice,sliceLowerLimits,sliceUpperLimits);
    -- Get the H atoms (including those not in the slice) for the entire system  
    hCloud=getPointCloudAtomsOfOneAtomType(resCloud,hCloud,hydrogenAtomType,false,{0,0,0}, {0,0,0});
    nList=neighborList(cutoffRadius, oCloud, oxygenAtomType); --- Calculate the neighborlist by ID (check this)
    hbnList=getHbondNetworkFromClouds(oCloud,hCloud,nList); --- Get the hydrogen-bonded network for the current frame
    hbnList=bondNetworkByIndex(oCloud,hbnList); --- Hydrogen-bonded network using indices not IDs
    rings=getPrimitiveRings(hbnList,maxDepth); --- Gets every ring (non-primitives included)
    bulkRingNumberAnalysis(outDir, rings, hbnList, oCloud, maxDepth, targetFrame); --- Writes out primitive rings for a bulk system 
end
