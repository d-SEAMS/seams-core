-- Ice Ih and sII census through require("dseams").
-- Ice Ih: 128 hexagonal, 256 six-rings.
-- sII: 136 clathrate, 144 five-rings, 16 six-rings.
local dseams = require("dseams")

dseams.set_num_threads(1)

local function census(path)
  local cloud = dseams.read(path, {all = true})
  local chill = dseams.chill_plus(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED, graph = "cutoff"})
  local rings = dseams.rings(cloud, {cutoff = 3.5, graph = "cutoff", max_depth = 6})
  cloud:close()
  return chill.counts, rings.by_size
end

local ice, ice_rings = census("input/traj/ice1h.lammpstrj")
assert(ice.hexagonal == 128, "ice Ih hexagonal")
assert(ice_rings[6] == 256, "ice Ih six-rings")
local sii, sii_rings = census("input/traj/genice_sII.lammpstrj")
assert(sii.clathrate == 136, "sII clathrate")
assert(sii_rings[5] == 144, "sII five-rings")
assert(sii_rings[6] == 16, "sII six-rings")
print("lua hexagonal 128 rings6 256 clathrate 136 rings5 144")
