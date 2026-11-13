-- Copyright (c) d-SEAMS developers
-- Cutoff neighbour rows. Each row leads with the atom, matching neighListO.
local dseams = require("dseams")

local root = arg[1] or "."
local cloud = dseams.read(root .. "/input/traj/genice_sI.lammpstrj", {type = 1})
local rows = dseams.neighbors(cloud, {cutoff = 3.5, type = 1, graph = "cutoff"})
assert(#rows == cloud.nop)
local nnz = 0
for i, row in ipairs(rows) do
  assert(row[1] == cloud._ids[i])
  nnz = nnz + #row - 1
end
print(string.format("neighbors rows=%d nnz=%d", #rows, nnz))
cloud:close()
