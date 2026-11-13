-- Copyright (c) d-SEAMS developers
-- Primitive-ring census on the oxygen network of the mixed TIP4P frame.
local dseams = require("dseams")

local root = arg[1] or "."
local cloud = dseams.read(root .. "/input/traj/exampleTraj.lammpstrj", {type = 2})
local rings = dseams.ringNetwork(cloud, 8)
local parts = {}
for size = 3, 8 do
  parts[#parts + 1] = string.format("%d %d", size, rings.by_size[size] or 0)
end
print(string.format("rings %s fraction %.6f", table.concat(parts, " "), rings.fraction))
cloud:close()
