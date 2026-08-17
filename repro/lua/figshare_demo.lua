-- Figshare v1 deposits through require("dseams").
-- Env: FIGSHARE_TRAJ, FIGSHARE_OUTDIR, FIGSHARE_FRAME (optional).
-- Arg 1: demo key (chillPlus, bulkTopologicalCriterion, iceNanotube, monolayer, rdf2D).

local dseams = require("dseams")
local core = dseams.core

local key = arg[1]
local traj = os.getenv("FIGSHARE_TRAJ")
local outdir = os.getenv("FIGSHARE_OUTDIR") or "."
local frame = tonumber(os.getenv("FIGSHARE_FRAME") or "1")
assert(key and traj, "usage: FIGSHARE_TRAJ=... lua run_demo.lua <key>")

local function ntrue(t)
  local n = 0
  if t == nil then
    return 0
  end
  for i = 1, #t do
    if t[i] then
      n = n + 1
    end
  end
  return n
end

local function count_types(types)
  local c = {}
  for i = 1, #types do
    local name = types[i]
    c[name] = (c[name] or 0) + 1
  end
  return c
end

local function emit(fields)
  local parts = {string.format("demo=%s", key)}
  local keys = {}
  for k in pairs(fields) do
    keys[#keys + 1] = k
  end
  table.sort(keys)
  for i = 1, #keys do
    parts[#parts + 1] = string.format("%s=%s", keys[i], tostring(fields[keys[i]]))
  end
  print(table.concat(parts, " "))
end

if key == "chillPlus" then
  local cloud = dseams.read(traj, {type = 1, frame = frame})
  assert(cloud.nop > 0, "empty cloud")
  local types = dseams.chill_plus(cloud, {cutoff = 3.5, type = 1})
  local counts = count_types(types)
  local cages = dseams.cages(cloud, {type = 1, cutoff = 5.0})
  local n_ddc = ntrue(cages.ddc)
  local n_hc = ntrue(cages.hc)
  assert((counts.cubic or 0) == cloud.nop, "Ic lattice is not all cubic")
  assert(n_ddc > 0 and n_hc == 0, "Ic cages are not all DDC")
  emit({
    nop = cloud.nop,
    cubic = counts.cubic or 0,
    n_ddc = n_ddc,
    n_hc = n_hc,
  })
elseif key == "bulkTopologicalCriterion" then
  local cloud = dseams.read(traj, {type = 1, frame = frame})
  assert(cloud.nop > 0, "empty cloud")
  local cages = dseams.cages(cloud, {type = 1, cutoff = 5.0})
  local n_ddc = ntrue(cages.ddc)
  local n_hc = ntrue(cages.hc)
  assert(n_ddc + n_hc > 0, "crystallized frame has no cage labels")
  emit({
    nop = cloud.nop,
    frame = frame,
    n_ddc = n_ddc,
    n_hc = n_hc,
  })
elseif key == "iceNanotube" then
  local cloud = dseams.read(traj, {type = 2, frame = frame})
  assert(cloud.nop > 0, "empty cloud")
  local nList = dseams.neighbors(cloud, {cutoff = 3.5, type = 2})
  local hbn = core.getHbondNetwork(traj, cloud, nList, frame, 1)
  hbn = core.bondNetworkByIndex(cloud, hbn)
  local rings = core.getPrimitiveRings(hbn, 6)
  core.prismAnalysis(outdir .. "/", rings, hbn, cloud, 6, 1, frame, frame, false)
  emit({nop = cloud.nop, n_rings = #rings})
elseif key == "monolayer" then
  -- x in [0, 50]; equal y/z bounds leave those axes unconstrained
  -- (same region as the Python monolayer notebook).
  local cloud = core.readLammpsTrjO(
    traj, frame, 2, true, {0.0, 0.0, 0.0}, {50.0, 0.0, 0.0}
  )
  assert(
    cloud.nop > 250 and cloud.nop < 400,
    string.format("slice did not select the sheet (nop=%d)", cloud.nop)
  )
  local nList = dseams.neighbors(cloud, {cutoff = 3.5, type = 2})
  local hbn = core.getHbondNetwork(traj, cloud, nList, frame, 1)
  hbn = core.bondNetworkByIndex(cloud, hbn)
  local rings = core.getPrimitiveRings(hbn, 4)
  core.ringAnalysis(outdir .. "/", rings, hbn, cloud, 4, 50.0 * 50.0, frame)
  emit({nop = cloud.nop, n_rings = #rings})
elseif key == "rdf2D" then
  local cloud = core.readLammpsTrjO(
    traj, frame, 2, true, {0.0, 0.0, 0.0}, {50.0, 0.0, 0.0}
  )
  assert(cloud.nop > 0, "empty cloud")
  local rdf = core.calcRDF3D(cloud, 2, 2, 12.0, 80)
  local gmax, rpeak = 0.0, 0.0
  for i = 1, #rdf.g do
    if rdf.r[i] < 4.0 and rdf.g[i] > gmax then
      gmax = rdf.g[i]
      rpeak = rdf.r[i]
    end
  end
  assert(gmax > 1.5, "O-O first peak is missing")
  emit({nop = cloud.nop, gmax = string.format("%.4f", gmax), rpeak = string.format("%.4f", rpeak)})
else
  error("unknown demo " .. tostring(key))
end
