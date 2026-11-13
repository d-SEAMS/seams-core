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

local function nop_line(tag, path, options)
  local cloud = dseams.read(path, options)
  io.write(string.format("read %s %d\n", tag, cloud.nop))
  cloud:close()
end

local function mean(values)
  local total = 0
  for _, value in ipairs(values) do
    total = total + value
  end
  return total / #values
end

local function extra(name, cloud)
  local chill = dseams.chill(cloud, {cutoff = 3.5, type = 0, graph = "cutoff", rule = dseams.PUBLISHED})
  io.write(name, " chill0")
  for _, key in ipairs({
    "cubic", "hexagonal", "water", "interfacial", "clathrate",
    "interClathrate", "unclassified", "reCubic", "reHex",
  }) do
    io.write(" ", chill.counts[key])
  end
  io.write("\n")
  local q3 = dseams.steinhardtQl(cloud, {cutoff = 3.5, l = 3})
  local q6 = dseams.steinhardtQl(cloud, {cutoff = 3.5, l = 6, average = true})
  local q4 = dseams.steinhardtQl(cloud, {cutoff = 3.5, l = 4})
  local q8 = dseams.steinhardtQl(cloud, {cutoff = 3.5, l = 8})
  local q12 = dseams.steinhardtQl(cloud, {cutoff = 3.5, l = 12})
  io.write(string.format(
    "%s q3 %.12g\n%s q4 %.12g\n%s q6bar %.12g\n%s q8 %.12g\n%s q12 %.12g\n",
    name, mean(q3.ql), name, mean(q4.ql), name, mean(q6.qlBar), name, mean(q8.ql),
    name, mean(q12.ql)))
  local rings = dseams.rings(cloud, {cutoff = 3.5, max_depth = 6, graph = "cutoff"})
  local first = rings.vertices[1] or {}
  io.write(string.format("%s vert %d %d\n", name, #first, first[1] or -1))
  local cages = dseams.cages(cloud, {cutoff = 3.5, graph = "seeded", k = 4})
  io.write(string.format("%s cages %d %d %d %d\n", name, cages.n_hex, cages.n_cubic, cages.n_water, cages.n_both))
  local knn = dseams.knn(cloud, {k = 4, cutoff = 3.5, mutual = true})
  local nnz = 0
  for _, row in ipairs(knn) do
    nnz = nnz + #row - 1
  end
  io.write(string.format("%s knn %d\n", name, nnz))
  local rdf = dseams.rdf(cloud, {rmax = 6.0, bins = 8})
  io.write(string.format("%s rdf %d %.12g %.12g\n", name, #rdf.r, rdf.volume, rdf.g[1] or 0))
  local running = dseams.running_cn(cloud, {rmax = 6.0, bins = 8})
  io.write(string.format("%s rcn %.12g\n", name, running.cn[#running.cn] or 0))
  io.write(string.format("%s cn %.12g\n", name, dseams.cn(cloud, {cutoff = 3.5, bins = 8})))
  io.write(string.format("%s hbond %d\n", name, dseams.hbonds(cloud, {oxygen_type = 1, h_type = 2})))
  local density = dseams.density(cloud, {type = 0, bins = 4, axis = 2})
  io.write(string.format("%s rho %d %.12g\n", name, #density.rho, density.rho[1] or 0))
  local fp = dseams.fingerprint(cloud, {cutoff = 3.5, hops = 2})
  local ring_sum = 0
  for _, count in ipairs(fp.rings) do
    ring_sum = ring_sum + count
  end
  io.write(string.format("%s fp %d %d\n", name, #fp.key, ring_sum))
  local kinds = dseams.site_table(cloud, "1=waterO")
  io.write(string.format("%s site %d\n", name, kinds[1]))
  local pairs = dseams.pairs(cloud, {spec = "1=cationHead"})
  io.write(string.format("%s pairs %d %d %d\n", name, pairs.count, pairs.n_cation, pairs.n_anion))
  local domain = dseams.domain(cloud, {cutoff = 3.5, spec = "1=waterO", subset = 0})
  io.write(string.format("%s domain %d %d %.12g\n", name, domain.n, domain.largest, domain.percolation))
  local prisms = dseams.prisms(cloud, {cutoff = 3.5})
  io.write(string.format("%s prism %d\n", name, prisms.perfect[7] or 0))
  local sheet = dseams.monolayer(cloud, {cutoff = 3.5, max_depth = 6, sheet_area = 100})
  io.write(string.format("%s sheet %d\n", name, sheet.by_size[6] or 0))
  cloud:close()
end

report("ih", dseams.read(frame, {all = true}), 6)
report("si", dseams.read(root .. "/input/traj/genice_sI.lammpstrj", {type = 1}), 6)
report("interface", dseams.read(root .. "/input/traj/interface_f0.lammpstrj", {type = 1}), 6)
report("tip", dseams.read(root .. "/input/traj/exampleTraj.lammpstrj", {type = 2}), 8)

local xyz = (frame:match("^(.*)/") or ".") .. "/sample.xyz"
nop_line("si", root .. "/input/traj/genice_sI.lammpstrj", {type = 1})
nop_line("tip", root .. "/input/traj/exampleTraj.lammpstrj", {type = 2})
nop_line("xyz", xyz, {all = true})
nop_line("con", root .. "/input/con/tiny_multi_cuh2.con", {all = true})
nop_line("tri", frame, {all = true})

local guest = dseams.read(root .. "/input/traj/genice_sI.lammpstrj", {type = 1})
local signed = dseams.cages_by_signature(guest, "512", {cutoff = 3.5, graph = "cutoff"})
guest:close()
io.write(string.format("si sig %d %d\n", signed.n_cages, signed.n_atoms))

extra("ih", dseams.read(frame, {all = true}))
