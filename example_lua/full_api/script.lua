print("full_api: exercising the registered Lua surface")

local function histogram(t)
  local counts = {}
  for i = 1, #t do
    counts[t[i]] = (counts[t[i]] or 0) + 1
  end
  local parts = {}
  for name, n in pairs(counts) do
    parts[#parts + 1] = string.format("%s=%d", name, n)
  end
  table.sort(parts)
  return table.concat(parts, " ")
end

local scratch = (os.getenv("TMPDIR") or "/tmp")
  .. "/dseams_full_api_" .. tostring(os.time()) .. "/"

-- Readers -------------------------------------------------------------------
local water = readLammpsTrjO(trajectory, targetFrame, oxygenAtomType)
local waterReduced = readLammpsTrjreduced(trajectory, targetFrame,
  oxygenAtomType, false, {0, 0, 0}, {0, 0, 0})
local cube = readXYZ(xyzFile)
assert(water.nop == waterReduced.nop)
assert(cube.nop == 8)
local box = water:box()
print(string.format("read water_nop=%d xyz_nop=%d frame=%d box={%.2f, %.2f, %.2f}",
  water.nop, cube.nop, water.currentFrame, box[1], box[2], box[3]))

-- Error propagation: a bad path raises a Lua error, caught by pcall ---------
local ok, err = pcall(readXYZ, "no/such/file.xyz")
assert(not ok and err ~= nil)
print("pcall_missing_file caught=true")

-- Neighbour lists -----------------------------------------------------------
local nlWater = neighListO(cutoffRadius, water, oxygenAtomType)
local nlWaterIdx = neighbourListByIndex(water, nlWater)
local nlWaterIdx2 = getNewNeighbourListByIndex(water, cutoffRadius)
assert(#nlWaterIdx == water.nop and #nlWaterIdx2 == water.nop)
print(string.format("neighbours rows=%d row1_len=%d row1_len_direct=%d",
  #nlWater, #nlWaterIdx[1], #nlWaterIdx2[1]))

-- Rings ---------------------------------------------------------------------
local rings = ringNetwork(nlWaterIdx, maxDepth)
local updater = RingUpdater.new(maxDepth)
local ringsFirst = updater:update(nlWaterIdx)
local sourcesFirst = updater:lastRecomputedSources()
local ringsSecond = updater:update(nlWaterIdx)
assert(#ringsFirst == #rings and #ringsSecond == #rings)
assert(updater:lastRecomputedSources() == 0)
print(string.format("rings n=%d updater_first_sources=%d repeat_sources=%d repeat_balls=%d",
  #rings, sourcesFirst, updater:lastRecomputedSources(),
  updater:lastBallsRefreshed()))

-- CHILL+ on the water frame -------------------------------------------------
getCorrelPlus(water, nlWater, false)
local plusTypes = getIceTypePlus(water, nlWater, scratch, targetFrame, false,
  "chillPlus.txt")
assert(#plusTypes == water.nop)
print("chillPlus " .. histogram(plusTypes))

-- Plain CHILL on the mW frame -----------------------------------------------
local mw = readLammpsTrjO(mwTrajectory, targetFrame, mwAtomType)
local nlMw = neighListO(cutoffRadius, mw, mwAtomType)
getCorrel(mw, nlMw, false)
local chillTypes = getIceType(mw, nlMw, scratch, targetFrame, false,
  "chill.txt")
assert(#chillTypes == mw.nop)
assert(#mw:iceTypes() == mw.nop)
print("chill " .. histogram(chillTypes))

-- Steinhardt, plain and Voronoi-weighted ------------------------------------
-- The candidate cutoff bounds the bisector-plane search and must exceed the
-- largest facet neighbour distance; twice the first-shell distance is ample.
local voroCut = 2 * cutoffRadius
for _, l in ipairs({3, 4, 6, 8}) do
  local q = steinhardtQl(mw, nlMw, l)
  local qv = steinhardtQlVoronoi(mw, voroCut, l)
  assert(#q.ql == mw.nop and #q.qlBar == mw.nop and #qv.ql == mw.nop)
  print(string.format("steinhardt l=%d ql1=%.4f qlBar1=%.4f voronoi_ql1=%.4f",
    l, q.ql[1], q.qlBar[1], qv.ql[1]))
end

-- Voronoi facet weights -----------------------------------------------------
local cells = voronoiFacetWeights(mw, voroCut)
local weightSum = 0.0
for i = 1, #cells[1].weights do
  weightSum = weightSum + cells[1].weights[i]
end
print(string.format("voronoi cells=%d facets1=%d weight_sum1=%.3f",
  #cells, #cells[1].neighbours, weightSum))

-- Hydrogen bonds with explicit geometric thresholds --------------------------
-- The defaults are the water criterion (2.42 A, 30 deg); a cutoff below any
-- O-H separation must produce an edge-free network
local hbWater = getHbondNetwork(trajectory, water, nlWater, targetFrame, 1)
local hbNone = getHbondNetwork(trajectory, water, nlWater, targetFrame, 1,
  0.1, 30.0)
local edges, noneEdges = 0, 0
for i = 1, #hbWater do
  edges = edges + #hbWater[i] - 1
  noneEdges = noneEdges + #hbNone[i] - 1
end
assert(edges > 0 and noneEdges == 0)
print(string.format("hbonds default_edges=%d tight_cutoff_edges=%d",
  edges, noneEdges))

-- Bond-classification registry ----------------------------------------------
-- CHILL and CHILL+ are registered rule sets; scripts add their own
local ruleNames = bondClassifierNames()
local haveChill, haveChillPlus = false, false
for i = 1, #ruleNames do
  if ruleNames[i] == "CHILL" then haveChill = true end
  if ruleNames[i] == "CHILL+" then haveChillPlus = true end
end
assert(haveChill and haveChillPlus)

local mwByName = readLammpsTrjO(mwTrajectory, targetFrame, mwAtomType)
classifyBonds(mwByName, nlMw, "CHILL")
local namedTypes = getIceType(mwByName, nlMw, scratch, targetFrame, false,
  "chillNamed.txt")
assert(#namedTypes == mw.nop)

registerBondClassifier("everything-staggered",
  { staggeredMax = 1.0, eclipsedMin = 2.0, eclipsedMax = 3.0,
    coordinationNumber = 4 })
local mwCustom = readLammpsTrjO(mwTrajectory, targetFrame, mwAtomType)
classifyBonds(mwCustom, nlMw, "everything-staggered")
print(string.format("bond_rules n=%d chill_by_name=%s custom_registered=true",
  #ruleNames, histogram(namedTypes)))

-- The .con reader, when the readcon-core backend is compiled in --------------
if readCon ~= nil then
  local con = readCon(conFile, 2)
  assert(con.nop == 4 and con.currentFrame == 2)
  print(string.format("readcon nop=%d frame=%d", con.nop, con.currentFrame))
else
  print("readcon backend absent; readCon not registered")
end

-- Structure descriptors -----------------------------------------------------
local hits = classifyTemplates(mw, nlMw, 12)
local names = {}
for i = 1, #hits do
  names[i] = hits[i].name
end
print("templates " .. histogram(names))

local spec = soapSpectrum(mw, 0, nlMw, 3, 6, cutoffRadius)
local specAll = soapSpectrumAll(water, nlWater, 2, 2, cutoffRadius)
local vFeat = voronoiFeatures(water, cutoffRadius)
assert(#specAll == water.nop and #vFeat == water.nop)
print(string.format("soap len=%d all_rows=%d voronoi_feat_len=%d",
  #spec, #specAll, #vFeat[1]))

-- Clean up the scratch output -----------------------------------------------
os.remove(scratch .. "bop/chillPlus.txt")
os.remove(scratch .. "bop/chill.txt")
os.remove(scratch .. "bop/chillNamed.txt")
os.remove(scratch .. "bop")
os.remove(scratch)

print("full_api: done")
