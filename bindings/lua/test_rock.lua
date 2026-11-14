-- Copyright (c) d-SEAMS developers
-- SPDX-License-Identifier: MIT
-- `luarocks test` for the dseams rock: modules installed by the rockspec,
-- loaded the way a user loads them. arg[1] is the repository root; the
-- ice Ih check is skipped when its input/ is not there.
local dseams = require("dseams")
local root = arg and arg[1] or "."

local a = 2.75 * 4.0 / math.sqrt(3.0)
local fcc = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}
local basis = {{0, 0, 0}, {0.25, 0.25, 0.25}}
local xyz = {}
for i = 0, 1 do
  for j = 0, 1 do
    for k = 0, 1 do
      for _, site in ipairs(fcc) do
        for _, shift in ipairs(basis) do
          xyz[#xyz + 1] = (i + site[1] + shift[1]) * a
          xyz[#xyz + 1] = (j + site[2] + shift[2]) * a
          xyz[#xyz + 1] = (k + site[3] + shift[3]) * a
        end
      end
    end
  end
end
local edge = 2 * a
local cloud = dseams.system(xyz, {0, 0, 0, edge, edge, edge, 0, 0, 0})
assert(cloud.nop == 64, cloud.nop)
local labels = dseams.chill_plus(cloud, {cutoff = 3.5, rule = dseams.PUBLISHED})
assert(labels.counts.cubic == 64, labels.counts.cubic)
local q = dseams.q4q6(cloud, {cutoff = 3.5})
assert(#q.q6 == 64)
local rings = dseams.ringNetwork(cloud, 6)
assert((rings.by_size[6] or 0) > 0)
cloud:close()
print(string.format("dseams %s %s cubic %d six-rings %d",
  dseams.version(), jit and jit.version or _VERSION, labels.counts.cubic, rings.by_size[6]))

assert(require("yoda") == dseams)

local frame = root .. "/input/traj/ice1h.lammpstrj"
local probe = io.open(frame)
if probe then
  probe:close()
  local ice = dseams.read(frame)
  local ih = dseams.chill_plus(ice, {cutoff = 3.5, rule = dseams.PUBLISHED})
  assert(ice.nop == 128 and ih.counts.hexagonal == 128, ih.counts.hexagonal)
  print("dseams ice Ih hexagonal " .. ih.counts.hexagonal)
  ice:close()
else
  print("dseams " .. frame .. " not found; ice Ih check skipped")
end

local has_fennel, fennel = pcall(require, "fennel")
if has_fennel then
  fennel.path = package.path:gsub("%.lua", ".fnl") .. ";" .. fennel.path
  local fnl = fennel.dofile(assert(package.searchpath("dseams_fnl", fennel.path)))
  assert(fnl["chill-plus"] == dseams.chill_plus)
  print("dseams_fnl " .. fennel.version)
else
  print("fennel not installed; dseams_fnl check skipped")
end
print("dseams rock ok")
