--- d-SEAMS Lua front.
-- Every call goes through the C API. Ice Ih
-- (`input/traj/ice1h.lammpstrj`, every atom) is 128 hexagonal waters
-- and 256 primitive six-rings. `docs/examples/core_counts.lua` asserts
-- both counts.
-- @module dseams

local core = require("dseams_core")

local dseams = {
  core = core,
  PUBLISHED = core.PUBLISHED,
  OVITO = core.OVITO,
  CUTOFF = core.CUTOFF,
  MUTUAL = core.MUTUAL,
}

local cloud_mt = {}
cloud_mt.__index = cloud_mt

function cloud_mt:box()
  return self._box
end

function cloud_mt:close()
  if self._sys ~= nil then
    self._sys:close()
    self._sys = nil
  end
end

function cloud_mt:__close()
  self:close()
end

local function wrap(sys, xyz, box, types, ids)
  return setmetatable({
    nop = sys:count(),
    _sys = sys,
    _xyz = xyz,
    _box = box,
    _types = types,
    _ids = ids,
  }, cloud_mt)
end

local function opts(value)
  if value == nil then
    return {}
  end
  return value
end

local function graph_of(value)
  if value == nil or value == "" or value == "cutoff" or value == core.CUTOFF then
    return core.CUTOFF
  end
  if value == "mutual" or value == "knn" or value == core.MUTUAL then
    return core.MUTUAL
  end
  if value == "seeded" or value == core.SEEDED then
    return core.SEEDED
  end
  error("graph is cutoff, mutual, or seeded", 3)
end

local function rule_of(value)
  if type(value) == "table" then
    value = value.rule
  end
  if value == nil or value == "" or value == "chill+" or value == "published" or value == core.PUBLISHED then
    return core.PUBLISHED
  end
  if value == "ovito-chill+" or value == "ovito" or value == core.OVITO then
    return core.OVITO
  end
  error("rule is published or ovito", 3)
end

local label_keys = {
  "cubic",
  "hexagonal",
  "water",
  "interfacial",
  "clathrate",
  "interClathrate",
  "unclassified",
  "reCubic",
  "reHex",
}

local function named_counts(codes)
  local named = {}
  for i, key in ipairs(label_keys) do
    named[key] = codes[i] or 0
  end
  return named
end

--- dseams.system.
function dseams.system(xyz, box, types, ids)
  local sys = core.system()
  if xyz == nil then
    return wrap(sys, nil, nil, nil, nil)
  end
  sys:set_frame(xyz, box, types, ids)
  return wrap(sys, xyz, box, types, ids)
end

--- dseams.set_num_threads.
function dseams.set_num_threads(n)
  return core.set_num_threads(n)
end

--- dseams.num_threads.
function dseams.num_threads()
  return core.num_threads()
end

--- dseams.version.
function dseams.version()
  return core.version()
end

--- dseams.label_name.
function dseams.label_name(code)
  return core.label_name(code)
end

--- Read one frame. `type = -1` keeps every atom.
function dseams.read(path, o)
  o = opts(o)
  local filter = -1
  if o.all == true then
    filter = -1
  elseif o.type ~= nil then
    filter = o.type
  else
    filter = 0
  end
  local sys = core.system()
  sys:read(path, o.frame or 1, filter)
  return wrap(sys)
end

--- Number of frames in a trajectory.
function dseams.nframes(path)
  return core.nframes(path)
end

--- Alias dseams.readLammpsTrj.
dseams.readLammpsTrj = function(path, o)
  o = opts(o)
  o.all = true
  return dseams.read(path, o)
end

--- Alias dseams.readLammpsTrjO.
dseams.readLammpsTrjO = dseams.read

local function rows_of(nb, cloud)
  local offsets, cols = nb:csr()
  nb:close()
  local rows = {}
  local ids = cloud._ids
  for i = 1, cloud.nop do
    local row = {ids and ids[i] or i}
    local start_at = offsets[i] + 1
    local stop_at = offsets[i + 1]
    for k = start_at, stop_at do
      local index = cols[k] + 1
      row[#row + 1] = ids and ids[index] or index
    end
    rows[i] = row
  end
  return rows
end

--- dseams.neighbors.
function dseams.neighbors(cloud, o)
  o = opts(o)
  return rows_of(cloud._sys:neighbours(o.cutoff or 3.5, o.type or 0, graph_of(o.graph)), cloud)
end

--- Alias dseams.neighbours.
dseams.neighbours = dseams.neighbors
--- Alias dseams.neighListO.
dseams.neighListO = dseams.neighbors

--- dseams.knn.
function dseams.knn(cloud, o)
  o = opts(o)
  local mutual = o.mutual == false and 0 or 1
  return rows_of(cloud._sys:knn(o.k or 4, o.cutoff or 3.5, o.type or 0, mutual), cloud)
end

--- Alias dseams.kNearestNeighbourList.
dseams.kNearestNeighbourList = dseams.knn

local function ice_result(labels, counts)
  local names = {}
  for i, code in ipairs(labels) do
    names[i] = core.label_name(code)
  end
  names.codes = labels
  names.counts = named_counts(counts)
  return names
end

--- CHILL+ labels. On Ice Ih `counts.hexagonal` is 128.
function dseams.chill_plus(cloud, o)
  o = opts(o)
  local labels, counts = cloud._sys:chill_plus(
    o.cutoff or 3.5, o.type or 0, graph_of(o.graph), rule_of(o.rule))
  return ice_result(labels, counts)
end

--- Alias dseams.getIceTypePlus.
dseams.getIceTypePlus = dseams.chill_plus
--- Alias dseams.getCorrelPlus.
dseams.getCorrelPlus = dseams.chill_plus
--- Alias dseams.getIceTypePlusNoPrint.
dseams.getIceTypePlusNoPrint = dseams.chill_plus

--- dseams.chill.
function dseams.chill(cloud, o)
  o = opts(o)
  local labels, counts = cloud._sys:chill(
    o.cutoff or 3.5, o.type or 0, graph_of(o.graph), rule_of(o.rule))
  return ice_result(labels, counts)
end

--- Alias dseams.getIceType.
dseams.getIceType = dseams.chill
--- Alias dseams.getCorrel.
dseams.getCorrel = dseams.chill

--- dseams.q4q6.
function dseams.q4q6(cloud, o)
  o = opts(o)
  local q4, q6 = cloud._sys:q4q6(o.cutoff or 3.5, o.type or 0)
  return {q4 = q4, q6 = q6}
end

--- dseams.steinhardtQl.
function dseams.steinhardtQl(cloud, o)
  o = opts(o)
  local order = o.l or o.order or 6
  local average = (o.average or o.qlbar) and 1 or 0
  local ql, bar = cloud._sys:steinhardt(o.cutoff or 3.5, o.type or 0, order, average)
  return {ql = ql, qlBar = bar}
end

--- Primitive rings. On Ice Ih `by_size[6]` is 256.
function dseams.rings(cloud, o)
  o = opts(o)
  local depth = o.max_depth or o.max_ring or o.depth or 6
  local rings = cloud._sys:rings(o.cutoff or 3.5, o.type or 0, graph_of(o.graph), depth)
  local counts = rings:histogram(depth + 1)
  local by_size = {}
  for index, count in ipairs(counts) do
    by_size[index - 1] = count
  end
  local sizes, symmetry = rings:labels()
  local fraction = rings:fraction(0)
  local vertices = {}
  for index = 0, rings:count() - 1 do
    vertices[#vertices + 1] = rings:vertices(index)
  end
  rings:close()
  return {
    counts = counts,
    by_size = by_size,
    sizes = sizes,
    symmetry = symmetry,
    fraction = fraction,
    vertices = vertices,
    n = #sizes,
  }
end

--- dseams.ringNetwork.
function dseams.ringNetwork(cloud, depth)
  if type(depth) == "number" then
    return dseams.rings(cloud, {max_depth = depth})
  end
  return dseams.rings(cloud, depth)
end

--- Alias dseams.bulkPrimitiveRings.
dseams.bulkPrimitiveRings = dseams.rings
--- Alias dseams.getPrimitiveRings.
dseams.getPrimitiveRings = function(cloud, depth)
  return dseams.ringNetwork(cloud, depth).vertices
end

--- Hexagonal and cubic cage membership.
function dseams.cages(cloud, o)
  o = opts(o)
  local graph = o.graph or "seeded"
  local handle = cloud._sys:cages(
    o.cutoff or 3.5, o.type or 0, graph_of(graph), o.k or 4, o.complete and 1 or 0)
  local n_hex, n_cubic, n_water, n_both = handle:counts()
  local labels = handle:labels(cloud.nop)
  handle:close()
  return {
    n_hex = n_hex,
    n_cubic = n_cubic,
    n_water = n_water,
    n_both = n_both,
    labels = labels,
  }
end

--- Alias dseams.seededCageAffiliation.
dseams.seededCageAffiliation = dseams.cages

--- dseams.cages_by_signature.
function dseams.cages_by_signature(cloud, spec, o)
  o = opts(o)
  local n_cages, n_atoms = cloud._sys:signature(
    o.cutoff or 3.5, o.type or 0, graph_of(o.graph), spec)
  return {n_cages = n_cages, n_atoms = n_atoms}
end

--- Alias dseams.findBySignature.
dseams.findBySignature = dseams.cages_by_signature

--- dseams.rdf.
function dseams.rdf(cloud, o)
  o = opts(o)
  local r, g, volume = cloud._sys:rdf(
    o.cutoff or o.rmax or 12.0, o.bins or 0, o.type_i or 0, o.type_j or 0)
  return {r = r, g = g, volume = volume}
end

--- Alias dseams.calcRDF3D.
dseams.calcRDF3D = dseams.rdf

--- dseams.running_cn.
function dseams.running_cn(cloud, o)
  o = opts(o)
  local r, cn = cloud._sys:running_cn(
    o.cutoff or o.rmax or 12.0, o.bins or 0, o.type_i or 0, o.type_j or 0)
  return {r = r, cn = cn}
end

--- Alias dseams.calcRunningCN.
dseams.calcRunningCN = dseams.running_cn

--- dseams.cn.
function dseams.cn(cloud, o)
  o = opts(o)
  return cloud._sys:coordination(
    o.cutoff or o.rmax or 4.5, o.bins or 0, o.type_i or 0, o.type_j or 0)
end

--- Alias dseams.calcCN.
dseams.calcCN = dseams.cn

--- dseams.hbonds.
function dseams.hbonds(cloud, o)
  o = opts(o)
  return cloud._sys:hbonds(
    o.cutoff or 3.5, o.oxygen_type or o.type or 0, o.h_type or o.hydrogen_type or 1,
    o.dist or 2.42, o.angle or 30.0)
end

--- Alias dseams.getHbondNetwork.
dseams.getHbondNetwork = dseams.hbonds

--- dseams.density.
function dseams.density(cloud, o)
  o = opts(o)
  local axis = o.axis or 2
  if axis == "x" then axis = 0 elseif axis == "y" then axis = 1 elseif axis == "z" then axis = 2 end
  local z, rho = cloud._sys:density(o.type or 0, o.bins or 0, axis)
  return {centres = z, rho = rho, axis = axis}
end

--- Alias dseams.densityByType.
dseams.densityByType = dseams.density

--- dseams.fingerprint.
function dseams.fingerprint(cloud, o)
  o = opts(o)
  local key, rings = cloud._sys:fingerprint(
    o.cutoff or 3.5, o.type or 0, o.k or 4, graph_of(o.graph), o.hops or 2,
    o.colours and 1 or 0)
  return {key = key, rings = rings}
end

--- Alias dseams.topologyFingerprint.
dseams.topologyFingerprint = dseams.fingerprint

--- dseams.site_table.
function dseams.site_table(cloud, spec)
  if type(cloud) == "string" then
    error("site_table needs a frame and a spec", 2)
  end
  return cloud._sys:site_table(spec)
end

--- dseams.pairs.
function dseams.pairs(cloud, o)
  o = opts(o)
  local n_pairs, n_cation, n_anion = cloud._sys:contact_pairs(o.spec or o.table)
  return {count = n_pairs, n_cation = n_cation, n_anion = n_anion}
end

--- Alias dseams.contactPairs.
dseams.contactPairs = dseams.pairs

--- dseams.domain.
function dseams.domain(cloud, o)
  o = opts(o)
  local subset = o.subset
  if subset == nil then
    subset = (o.kind == "polar" or o.kind == 1) and 1 or 0
  end
  local n_subset, largest, percolation = cloud._sys:domains(
    o.cutoff or 3.5, o.spec or o.table, subset)
  return {n = n_subset, largest = largest, percolation = percolation}
end

--- Alias dseams.domainStats.
dseams.domainStats = dseams.domain

--- dseams.prisms.
function dseams.prisms(cloud, o)
  o = opts(o)
  local perfect, imperfect = cloud._sys:prisms(o.cutoff or 3.5, o.type or 0)
  return {perfect = perfect, imperfect = imperfect}
end

--- dseams.monolayer.
function dseams.monolayer(cloud, o)
  o = opts(o)
  local depth = o.max_depth or 6
  local counts = cloud._sys:monolayer(
    o.cutoff or 3.5, o.type or 0, depth, o.sheet_area or o.area or 1.0)
  local by_size = {}
  for index, count in ipairs(counts) do
    by_size[index - 1] = count
  end
  return {counts = counts, by_size = by_size}
end

--- Alias dseams.shapeMatching.
dseams.shapeMatching = dseams.monolayer

--- Alias dseams.SEEDED.
dseams.SEEDED = core.SEEDED
--- Alias dseams.WATER_O.
dseams.WATER_O = core.WATER_O
dseams["chill-plus"] = dseams.chill_plus
dseams["ring-network"] = dseams.ringNetwork
dseams["cages-by-signature"] = dseams.cages_by_signature
dseams["running-cn"] = dseams.running_cn
dseams["site-table"] = dseams.site_table
dseams["set-num-threads"] = dseams.set_num_threads

return dseams
