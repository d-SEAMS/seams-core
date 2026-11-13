-- Copyright (c) d-SEAMS developers
-- Same census as tests/lua_oracle.cpp, through require("dseams").
local dseams = require("dseams")

local frame = arg[1]
local root = arg[2]
dseams.set_num_threads(1)

local function report(name, cloud, depth)
  local labels = dseams.chill_plus(cloud, {
    cutoff = 3.5,
    type = 0,
    graph = "cutoff",
    rule = dseams.PUBLISHED,
  })
  local rings = dseams.rings(cloud, {
    cutoff = 3.5,
    type = 0,
    graph = "cutoff",
    max_depth = depth,
  })
  local q6 = dseams.steinhardtQl(cloud, {cutoff = 3.5, type = 0, l = 6})
  local rows = dseams.neighbors(cloud, {cutoff = 3.5, type = 0, graph = "cutoff"})
  local total = 0
  local nnz = 0
  for _, value in ipairs(q6.ql) do
    total = total + value
  end
  for _, row in ipairs(rows) do
    nnz = nnz + #row - 1
  end
  local order = {
    "cubic", "hexagonal", "water", "interfacial", "clathrate",
    "interClathrate", "unclassified", "reCubic", "reHex",
  }
  io.write(string.format("%s n %d\n%s chill", name, cloud.nop, name))
  for _, key in ipairs(order) do
    io.write(" ", labels.counts[key])
  end
  io.write("\n", name, " rings")
  for size = 0, depth do
    io.write(" ", rings.by_size[size] or 0)
  end
  io.write(string.format("\n%s q6 %.12g\n%s nnz %d\n", name, total / cloud.nop, name, nnz))
  cloud:close()
end

report("ih", dseams.read(frame), 6)
report("si", dseams.read(root .. "/input/traj/genice_sI.lammpstrj", {type = 1}), 6)
report("interface", dseams.read(root .. "/input/traj/interface_f0.lammpstrj", {type = 1}), 6)
report("tip", dseams.read(root .. "/input/traj/exampleTraj.lammpstrj", {type = 2}), 8)
