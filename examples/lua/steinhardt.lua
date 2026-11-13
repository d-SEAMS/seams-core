-- Copyright (c) d-SEAMS developers
-- Local q6 on the cutoff shell. Degrees other than 4 and 6 are a C API gap.
local dseams = require("dseams")

local root = arg[1] or "."
local cloud = dseams.read(root .. "/input/traj/genice_sI.lammpstrj", {type = 1})
local q6 = dseams.steinhardtQl(cloud, {cutoff = 3.5, type = 1, l = 6})
local total = 0
for _, value in ipairs(q6.ql) do
  total = total + value
end
print(string.format("q6 %.12g", total / cloud.nop))
cloud:close()
