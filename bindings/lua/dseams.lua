-- Copyright (c) d-SEAMS developers
-- SPDX-License-Identifier: MIT
--
-- 1.x names on top of dseams_core. The core module calls seams.h.
-- A missing C entry raises a Lua error that quotes the signature
-- recorded in docs/orgmode/explanation/c-api-gaps.org.
-- Orthorhombic LAMMPS dumps and "# dseams-frame" files are parsed
-- here and passed to seams_system_set_frame.

local core = require("dseams_core")

local gaps = {
  read = "seams_system_read(seams_system *sys, const char *path, int frame, int type_filter)",
  chill = "seams_chill(seams_system *sys, double cutoff, int type, int graph, int *labels, int counts9[9])",
  steinhardt = "seams_steinhardt(seams_system *sys, double cutoff, int type, int order_l, double *ql, double *qlbar)",
  cages = "seams_cages(seams_system *sys, double cutoff, int type, int k, int graph, int complete, int *per_atom, int *n_hex, int *n_cubic, int *n_water, int *n_both)",
  signature = "seams_signature_cages(seams_system *sys, double cutoff, int type, int graph, const char *spec, const char *edge, int *n_cages, int *n_atoms)",
  knn = "seams_neighbours_knn(seams_system *sys, int k, double cutoff, int type, int mutual, seams_neighbours **out)",
  rdf = "seams_rdf(const seams_system *sys, double rmax, int nbins, int type_i, int type_j, double *r, double *g, int *count)",
  running_cn = "seams_running_cn(const seams_system *sys, double rmax, int nbins, int type_i, int type_j, double *r, double *cn)",
  cn = "seams_coordination(const seams_system *sys, int type_i, int type_j, double cutoff, double *cn)",
  hbonds = "seams_hbonds(seams_system *sys, double cutoff, int oxygen_type, int hydrogen_type, int *degree)",
  density = "seams_density_axis(const seams_system *sys, int type, int axis, int nbins, double *centres, double *rho)",
  fingerprint = "seams_fingerprint(seams_system *sys, double cutoff, int type, int hops, int max_ring, int *n_classes)",
  pairs = "seams_contact_pairs(seams_system *sys, const int *kind_of_type, int n_types, int *n_pairs)",
  domain = "seams_domains(seams_system *sys, const int *kind_of_type, int n_types, int kind, double cutoff, int *n_sites, int *largest, double *percolation)",
  prisms = "seams_prisms(seams_system *sys, double cutoff, int type, int max_depth, int *sizes, int *perfect, int *imperfect, int *n_rows)",
  monolayer = "seams_monolayer(seams_system *sys, double cutoff, int type, int max_depth, double sheet_area, int *counts, int counts_len)",
  vertices = "seams_rings_vertices(const seams_rings *rings, int index, int *vertices, int cap, int *n_vertices)",
  site = "seams_site_table(const char *spec, int *types, int *kinds, int *n)",
}

local function missing(key)
  local sig = gaps[key] or key
  error(sig .. ": not in the public C API", 2)
end

local function stub(key)
  return function()
    missing(key)
  end
end

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
  error("graph is cutoff or mutual", 3)
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

function dseams.system(xyz, box, types, ids)
  local sys = core.system()
  if xyz == nil then
    return wrap(sys, nil, nil, nil, nil)
  end
  sys:set_frame(xyz, box, types, ids)
  return wrap(sys, xyz, box, types, ids)
end

function dseams.set_num_threads(n)
  return core.set_num_threads(n)
end

function dseams.num_threads()
  return core.num_threads()
end

function dseams.version()
  return core.version()
end

function dseams.label_name(code)
  return core.label_name(code)
end

local function fields_of(line)
  local out = {}
  for word in line:gmatch("%S+") do
    out[#out + 1] = word
  end
  return out
end

local function read_frame_text(text)
  local lines = {}
  for line in (text .. "\n"):gmatch("(.-)\n") do
    if line ~= "" and line:sub(1, 1) ~= "#" then
      lines[#lines + 1] = line
    end
  end
  local n = tonumber(lines[1])
  local box = {}
  for word in lines[2]:gmatch("%S+") do
    box[#box + 1] = tonumber(word)
  end
  if n == nil or #box ~= 9 or #lines ~= n + 2 then
    error("dseams-frame needs a count, nine box numbers, and that many x y z rows", 3)
  end
  local xyz = {}
  for i = 1, n do
    local x, y, z = lines[i + 2]:match("^%s*(%S+)%s+(%S+)%s+(%S+)")
    xyz[3 * i - 2] = tonumber(x)
    xyz[3 * i - 1] = tonumber(y)
    xyz[3 * i] = tonumber(z)
  end
  return xyz, box, n
end

local function header_columns(line)
  local words = fields_of(line)
  local cols = {}
  for i = 3, #words do
    cols[words[i]] = i - 2
  end
  return cols
end

local function read_lammps(text, frame_index, type_filter)
  local frames = {}
  local current = nil
  local mode = nil
  local columns = nil
  for line in (text .. "\n"):gmatch("(.-)\n") do
    if line:find("^ITEM:%s*TIMESTEP") then
      current = {atoms = {}, bounds = {}}
      frames[#frames + 1] = current
      mode = "skip"
      columns = nil
    elseif current == nil then
      -- preamble before the first frame
    elseif line:find("^ITEM:%s*NUMBER") then
      mode = "number"
    elseif line:find("^ITEM:%s*BOX") then
      mode = "box"
      if line:find("xy") or line:find("xz") or line:find("yz") then
        missing("read")
      end
    elseif line:find("^ITEM:%s*ATOMS") then
      mode = "atoms"
      columns = header_columns(line)
      if columns.x == nil or columns.y == nil or columns.z == nil then
        error("ATOMS header needs x y z columns", 3)
      end
    elseif mode == "number" then
      current.n = tonumber(line)
      mode = nil
    elseif mode == "box" then
      local words = fields_of(line)
      if #words >= 3 then
        missing("read")
      end
      if #words >= 2 then
        current.bounds[#current.bounds + 1] = {tonumber(words[1]), tonumber(words[2])}
      end
    elseif mode == "atoms" and line:find("%S") then
      local words = fields_of(line)
      local function num(name, fallback)
        local index = columns[name] or fallback
        return tonumber(words[index])
      end
      local atom = {
        id = num("id", 1),
        type = num("type", 2),
        x = num("x"),
        y = num("y"),
        z = num("z"),
      }
      if atom.x ~= nil and atom.y ~= nil and atom.z ~= nil then
        current.atoms[#current.atoms + 1] = atom
      end
    end
  end
  local frame = frames[frame_index]
  if frame == nil then
    error(string.format("frame %d is past the end of the dump", frame_index), 3)
  end
  if #frame.bounds ~= 3 then
    error("BOX BOUNDS needs three lines", 3)
  end
  local xyz, types, ids = {}, {}, {}
  for _, atom in ipairs(frame.atoms) do
    if type_filter == nil or atom.type == type_filter then
      xyz[#xyz + 1] = atom.x
      xyz[#xyz + 1] = atom.y
      xyz[#xyz + 1] = atom.z
      types[#types + 1] = atom.type or 1
      ids[#ids + 1] = atom.id or #types
    end
  end
  local box = {
    frame.bounds[1][1],
    frame.bounds[2][1],
    frame.bounds[3][1],
    frame.bounds[1][2] - frame.bounds[1][1],
    frame.bounds[2][2] - frame.bounds[2][1],
    frame.bounds[3][2] - frame.bounds[3][1],
    0.0,
    0.0,
    0.0,
  }
  return xyz, box, types, ids
end

local function nonempty(path, cloud, detail)
  if cloud.nop == 0 then
    error(string.format("%s has no atoms%s", path, detail or ""), 3)
  end
  return cloud
end

function dseams.read(path, o)
  o = opts(o)
  local handle = assert(io.open(path, "r"))
  local text = handle:read("*a")
  handle:close()
  if text:sub(1, 14) == "# dseams-frame" then
    local xyz, box = read_frame_text(text)
    return nonempty(path, dseams.system(xyz, box))
  end
  local suffix = (path:match("%.([^%./]+)$") or ""):lower()
  if suffix ~= "lammpstrj" and suffix ~= "dump" and not text:find("ITEM:%s*TIMESTEP") then
    missing("read")
  end
  local frame = o.frame or 1
  local function load(type_filter)
    local xyz, box, types, ids = read_lammps(text, frame, type_filter)
    if #types == 0 then
      return {nop = 0, close = function() end, _box = box}
    end
    return dseams.system(xyz, box, types, ids)
  end
  if o.all == true then
    return nonempty(path, load(nil), " for all atoms")
  end
  if o.type ~= nil then
    return nonempty(path, load(o.type), " for type " .. tostring(o.type))
  end
  local guessed = load(2)
  if guessed.nop > 0 then
    return guessed
  end
  guessed:close()
  return nonempty(path, load(1), " for types 2 or 1")
end

dseams.readLammpsTrj = function(path, o)
  o = opts(o)
  o.all = true
  return dseams.read(path, o)
end

dseams.readLammpsTrjO = function(path, o)
  return dseams.read(path, o)
end

local function rows_from(cloud, o)
  o = opts(o)
  local nb = cloud._sys:neighbours(o.cutoff or 3.5, o.type or 0, graph_of(o.graph))
  local offsets, cols = nb:csr()
  nb:close()
  local rows = {}
  local ids = cloud._ids
  for i = 1, cloud.nop do
    local self_id = ids and ids[i] or i
    local row = {self_id}
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

function dseams.neighbors(cloud, o)
  return rows_from(cloud, o)
end

dseams.neighbours = dseams.neighbors
dseams.neighListO = dseams.neighbors

function dseams.knn(cloud, o)
  o = opts(o)
  local k = o.k or 4
  local mutual = o.mutual ~= false
  if k ~= 4 or not mutual then
    missing("knn")
  end
  return rows_from(cloud, {
    cutoff = (o.cutoff or 5.0) - 1.5,
    type = o.type or 0,
    graph = "mutual",
  })
end

dseams.kNearestNeighbourList = dseams.knn

function dseams.chill_plus(cloud, o)
  o = opts(o)
  local labels, counts = cloud._sys:chill_plus(
    o.cutoff or 3.5,
    o.type or 0,
    graph_of(o.graph),
    rule_of(o.rule)
  )
  local names = {}
  for i, code in ipairs(labels) do
    names[i] = core.label_name(code)
  end
  names.codes = labels
  names.counts = named_counts(counts)
  return names
end

dseams.getIceTypePlus = dseams.chill_plus
dseams.getCorrelPlus = dseams.chill_plus
dseams.getIceTypePlusNoPrint = dseams.chill_plus

function dseams.q4q6(cloud, o)
  o = opts(o)
  local q4, q6 = cloud._sys:q4q6(o.cutoff or 3.5, o.type or 0)
  return {q4 = q4, q6 = q6}
end

function dseams.steinhardtQl(cloud, o)
  o = opts(o)
  local order = o.l or o.order or 6
  if o.average or o.qlbar or (order ~= 4 and order ~= 6) then
    missing("steinhardt")
  end
  local pair = dseams.q4q6(cloud, o)
  if order == 4 then
    return {ql = pair.q4}
  end
  return {ql = pair.q6}
end

function dseams.rings(cloud, o)
  o = opts(o)
  local depth = o.max_depth or o.max_ring or o.depth or 6
  local rings = cloud._sys:rings(
    o.cutoff or 3.5,
    o.type or 0,
    graph_of(o.graph),
    depth
  )
  local counts = rings:histogram(depth + 1)
  local by_size = {}
  for index, count in ipairs(counts) do
    by_size[index - 1] = count
  end
  local sizes, symmetry = rings:labels()
  local fraction = rings:fraction(0)
  local sym_names = {}
  for i, code in ipairs(symmetry) do
    sym_names[i] = core.ring_symmetry_name(code)
  end
  rings:close()
  return {
    counts = counts,
    by_size = by_size,
    sizes = sizes,
    symmetry = symmetry,
    symmetry_names = sym_names,
    fraction = fraction,
    n = #sizes,
  }
end

function dseams.ringNetwork(cloud, depth)
  if type(cloud) ~= "table" or cloud._sys == nil then
    missing("vertices")
  end
  if type(depth) == "number" then
    return dseams.rings(cloud, {max_depth = depth})
  end
  return dseams.rings(cloud, depth)
end

dseams.bulkPrimitiveRings = dseams.rings
dseams.getPrimitiveRings = function(rows, depth)
  if type(rows) ~= "table" or rows._sys == nil then
    missing("vertices")
  end
  return dseams.ringNetwork(rows, depth)
end

function dseams.chill()
  missing("chill")
end

dseams.getIceType = dseams.chill
dseams.getCorrel = dseams.chill

function dseams.cages()
  missing("cages")
end

dseams.seededCageAffiliation = dseams.cages

function dseams.cages_by_signature()
  missing("signature")
end

dseams.findBySignature = dseams.cages_by_signature

function dseams.rdf()
  missing("rdf")
end

dseams.calcRDF3D = dseams.rdf

function dseams.running_cn()
  missing("running_cn")
end

dseams.calcRunningCN = dseams.running_cn

function dseams.cn()
  missing("cn")
end

dseams.calcCN = dseams.cn

function dseams.hbonds()
  missing("hbonds")
end

dseams.getHbondNetwork = dseams.hbonds

function dseams.density()
  missing("density")
end

dseams.densityByType = dseams.density

function dseams.fingerprint()
  missing("fingerprint")
end

dseams.topologyFingerprint = dseams.fingerprint

function dseams.site_table()
  missing("site")
end

function dseams.pairs()
  missing("pairs")
end

dseams.contactPairs = dseams.pairs

function dseams.domain()
  missing("domain")
end

dseams.domainStats = dseams.domain

function dseams.prisms()
  missing("prisms")
end

function dseams.monolayer()
  missing("monolayer")
end

dseams.shapeMatching = dseams.monolayer

-- Fennel calls the same table. Kebab-case is the Fennel spelling.
dseams["chill-plus"] = dseams.chill_plus
dseams["ring-network"] = dseams.ringNetwork
dseams["cages-by-signature"] = dseams.cages_by_signature
dseams["running-cn"] = dseams.running_cn
dseams["site-table"] = dseams.site_table
dseams["set-num-threads"] = dseams.set_num_threads

local core_stubs = {
  readXYZ = "read",
  readCon = "read",
  readChemfiles = "read",
  readLammpsTrjreduced = "read",
  calcRDF3D = "rdf",
  calcRunningCN = "running_cn",
  calcCN = "cn",
  densityByType = "density",
  densityByKind = "density",
  contactPairs = "pairs",
  domainStats = "domain",
  mutualNearestUnlike = "pairs",
  getHbondNetwork = "hbonds",
  topologyFingerprint = "fingerprint",
  findBySignature = "signature",
  seededCageAffiliation = "cages",
  cageAffiliation = "cages",
  prismAnalysis = "prisms",
  ringAnalysis = "monolayer",
  bulkTopologicalNetworkCriterion = "cages",
  steinhardtQl = "steinhardt",
  getCorrel = "chill",
  getIceType = "chill",
}

for name, key in pairs(core_stubs) do
  if core[name] == nil then
    core[name] = stub(key)
  end
end

return dseams
