-- Neighbours, rings, and CHILL+ from Lua.
local dseams = require("dseams")

local path = os.getenv("DSEAMS_ROOT") .. "/input/traj/ice1h.lammpstrj"
local cloud = dseams.read(path, {type = -1})
local rows = dseams.neighbours(cloud, {cutoff = 3.5, graph = "cutoff"})
local pairs = 0
for _, row in ipairs(rows) do
  pairs = pairs + #row - 1
end
local rings = dseams.rings(cloud, {cutoff = 3.5, max_depth = 6, graph = "cutoff"})
local ice = dseams.chill_plus(cloud, {cutoff = 3.5, graph = "cutoff"})
local keys = {
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
local parts = {"pairs", pairs, "rings"}
for size = 0, 6 do
  parts[#parts + 1] = rings.by_size[size] or 0
end
parts[#parts + 1] = "chill"
for _, key in ipairs(keys) do
  parts[#parts + 1] = ice.counts[key] or 0
end
print(table.concat(parts, " "))
cloud:close()
