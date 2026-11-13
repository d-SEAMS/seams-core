-- Copyright (c) d-SEAMS developers
-- 1.x read of the mixed TIP4P frame, oxygen type 2.
local dseams = require("dseams")

local root = arg[1] or "."
local cloud = dseams.read(root .. "/input/traj/exampleTraj.lammpstrj", {type = 2})
assert(cloud.nop == 250, cloud.nop)
print(string.format("dseams.read nop=%d", cloud.nop))
cloud:close()
