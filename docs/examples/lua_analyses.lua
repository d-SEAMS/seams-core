-- Copyright (c) d-SEAMS developers
-- Every analysis require("dseams") forwards to seams.h, on the ice Ih frame.
local dseams = require("dseams")

dseams.set_num_threads(1)

local function mean(values)
  local total = 0
  for _, value in ipairs(values) do
    total = total + value
  end
  return total / #values
end

local cloud = dseams.read("input/traj/ice1h.lammpstrj", {all = true})
local chill = dseams.chill(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED})
print(string.format("hexagonal %d", chill.counts.hexagonal))
local q6 = dseams.steinhardtQl(cloud, {cutoff = 3.5, l = 6, average = true})
print(string.format("q6bar %.12g", mean(q6.qlBar)))
local cages = dseams.cages(cloud, {cutoff = 3.5, graph = "seeded", k = 4})
print(string.format("cages %d %d %d %d", cages.n_hex, cages.n_cubic, cages.n_water, cages.n_both))
local knn = dseams.knn(cloud, {k = 4, cutoff = 3.5, mutual = true})
local nnz = 0
for _, row in ipairs(knn) do
  nnz = nnz + #row - 1
end
print(string.format("knn %d", nnz))
local rdf = dseams.rdf(cloud, {rmax = 6.0, bins = 8})
print(string.format("rdf %d %.12g %.12g", #rdf.r, rdf.volume, rdf.g[1] or 0))
local running = dseams.running_cn(cloud, {rmax = 6.0, bins = 8})
print(string.format("rcn %.12g", running.cn[#running.cn] or 0))
print(string.format("cn %.12g", dseams.cn(cloud, {cutoff = 3.5, bins = 8})))
print(string.format("hbond %d", dseams.hbonds(cloud, {oxygen_type = 1, h_type = 2})))
local density = dseams.density(cloud, {type = 0, bins = 4, axis = 2})
print(string.format("rho %d %.12g", #density.rho, density.rho[1] or 0))
local fp = dseams.fingerprint(cloud, {cutoff = 3.5, hops = 2})
local ring_sum = 0
for _, count in ipairs(fp.rings) do
  ring_sum = ring_sum + count
end
print(string.format("fp %d %d", #fp.key, ring_sum))
local kinds = dseams.site_table(cloud, "1=waterO")
print(string.format("site %d", kinds[1]))
local pairs = dseams.pairs(cloud, {spec = "1=cationHead"})
print(string.format("pairs %d %d %d", pairs.count, pairs.n_cation, pairs.n_anion))
local domain = dseams.domain(cloud, {cutoff = 3.5, spec = "1=waterO", subset = 0})
print(string.format("domain %d %d %.12g", domain.n, domain.largest, domain.percolation))
local prisms = dseams.prisms(cloud, {cutoff = 3.5})
print(string.format("prism %d", prisms.perfect[7] or 0))
local sheet = dseams.monolayer(cloud, {cutoff = 3.5, max_depth = 6, sheet_area = 100})
print(string.format("sheet %d", sheet.by_size[6] or 0))
local rings = dseams.rings(cloud, {cutoff = 3.5, max_depth = 6})
local first = rings.vertices[1] or {}
print(string.format("vert %d %d", #first, first[1] or -1))
cloud:close()

local guest = dseams.read("input/traj/genice_sI.lammpstrj", {type = 1})
local signed = dseams.cages_by_signature(guest, "512", {cutoff = 3.5, graph = "cutoff"})
print(string.format("sig %d %d", signed.n_cages, signed.n_atoms))
guest:close()
