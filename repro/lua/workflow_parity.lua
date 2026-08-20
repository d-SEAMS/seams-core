-- Normalize the stable d-SEAMS workflows for the parity manifest.
-- Arguments: water trajectory, ice trajectory, ionic-site trajectory.

local dseams = require("dseams")

local water_path, ice_path, ions_path = arg[1], arg[2], arg[3]
assert(water_path and ice_path and ions_path, "three trajectory paths are required")

local labels = {
  "cubic", "hexagonal", "water", "interfacial", "clathrate",
  "interClathrate", "reCubic", "reHex", "unclassified",
}

local function emit(key, value)
  print(string.format("%s\t%s", key, tostring(value)))
end

local function emit_array(key, values)
  local out = {}
  for index = 1, #values do
    out[index] = string.format("%.17g", values[index])
  end
  emit(key, table.concat(out, ","))
end

local function count_labels(values, prefix)
  local counts = {}
  for _, label in ipairs(labels) do
    counts[label] = 0
  end
  for index = 1, #values do
    counts[values[index]] = (counts[values[index]] or 0) + 1
  end
  emit(prefix .. ".nop", #values)
  for _, label in ipairs(labels) do
    emit(prefix .. "." .. label, counts[label])
  end
end

local function cage_counts(score)
  local ih, ic, water = 0, 0, 0
  for index = 1, #score.hc do
    if score.hc[index] then
      ih = ih + 1
    elseif score.ddc[index] then
      ic = ic + 1
    else
      water = water + 1
    end
  end
  emit("cages.nop", #score.hc)
  emit("cages.ih", ih)
  emit("cages.ic", ic)
  emit("cages.water", water)
end

local function hbond_count(network)
  local directed = 0
  for index = 1, #network do
    directed = directed + math.max(0, #network[index] - 1)
  end
  return math.floor(directed / 2)
end

local ice = dseams.read(ice_path, {type = 1})
local water = dseams.read(water_path, {type = 2})
local ions = dseams.read(ions_path, {all = true})
emit("read.water_nop", water.nop)
emit("read.ice_nop", ice.nop)

count_labels(dseams.chill(ice, {type = 1, cutoff = 3.5}), "chill")
count_labels(dseams.chill_plus(ice, {type = 1, cutoff = 3.5}), "chill_plus")
cage_counts(dseams.cages(ice, {type = 1, cutoff = 5.0, k = 4}))

local rdf = dseams.rdf(water, {type_i = 2, type_j = 2, cutoff = 6.0, bins = 60})
emit_array("rdf.r", rdf.r)
emit_array("rdf.g", rdf.g)
emit("cn.value", dseams.cn(water, {type_i = 2, type_j = 2, cutoff = 6.0, bins = 60}))

local hbonds = dseams.hbonds(water, {
  path = water_path,
  frame = 1,
  type = 2,
  h_type = 1,
  cutoff = 3.5,
})
emit("hbonds.nop", water.nop)
emit("hbonds.count", hbond_count(hbonds))

local density = dseams.density(water, {type = 2, axis = "z", bins = 20})
emit_array("density.centres", density.centres)
emit_array("density.rho", density.rho)

local pairs = dseams.pairs(ions, {
  table = dseams.site_table("1=cationHead,2=anion"),
})
emit("pairs.count", pairs.count)
emit("pairs.n_cation", pairs.n_cation)
emit("pairs.n_anion", pairs.n_anion)

local domain = dseams.domain(ions, {
  table = dseams.site_table("1=polar,2=apolar"),
  kind = dseams.core.Kind.polar,
  cutoff = 1.5,
})
emit("domain.n", domain.n)
emit("domain.largest", domain.largest)
emit("domain.percolation", domain.percolation)
