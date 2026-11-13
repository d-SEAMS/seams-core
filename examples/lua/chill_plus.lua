-- Copyright (c) d-SEAMS developers
-- CHILL+ on the sI cage. The options table is examples/chill-plus.lua.
local dseams = require("dseams")

local root = arg[1] or "."
local rule = dofile(root .. "/examples/chill-plus.lua")
local cloud = dseams.read(root .. "/input/traj/genice_sI.lammpstrj", {type = 1})
local labels = dseams.chill_plus(cloud, {
  cutoff = rule.cutoff,
  type = 1,
  graph = "cutoff",
  rule = rule,
})
assert(#labels == cloud.nop)
local line = {}
for _, name in ipairs({
  "cubic", "hexagonal", "water", "interfacial", "clathrate",
  "interClathrate", "unclassified", "reCubic", "reHex",
}) do
  line[#line + 1] = string.format("%s %d", name, labels.counts[name])
end
print(table.concat(line, " "))
cloud:close()
