-- Copyright (c) d-SEAMS developers
-- A rejected C call is still a Lua error carrying the C message.
local dseams = require("dseams")

local ok, err = pcall(function()
  dseams.read("input/traj/genice_sI.lammpstrj", {frame = 0, type = 1})
end)
assert(not ok)
assert(tostring(err):find("frame index starts at 1", 1, true), tostring(err))

local cloud = dseams.read("input/traj/genice_sI.lammpstrj", {type = 1})
local bad, bad_err = pcall(function()
  dseams.chill(cloud, {graph = "seeded"})
end)
assert(not bad)
assert(tostring(bad_err):find("graph is cutoff or mutual", 1, true), tostring(bad_err))
cloud:close()
print("errors ok")
