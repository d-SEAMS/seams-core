-- Copyright (c) d-SEAMS developers
-- Eight-atom cubic cell. The C program docs/examples/ice.c labels every atom cubic.
local dseams = require("dseams")

local a = 2.75 * 4.0 / math.sqrt(3.0)
local fcc = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}
local basis = {{0, 0, 0}, {0.25, 0.25, 0.25}}
local xyz = {}
for _, site in ipairs(fcc) do
  for _, shift in ipairs(basis) do
    xyz[#xyz + 1] = (site[1] + shift[1]) * a
    xyz[#xyz + 1] = (site[2] + shift[2]) * a
    xyz[#xyz + 1] = (site[3] + shift[3]) * a
  end
end
local box = {0, 0, 0, a, a, a, 0, 0, 0}
local cloud = dseams.system(xyz, box)
local labels = dseams.chill_plus(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED})
assert(labels.counts.cubic == 8, labels.counts.cubic)
print(string.format("%s cubic %d", dseams.version(), labels.counts.cubic))
cloud:close()
