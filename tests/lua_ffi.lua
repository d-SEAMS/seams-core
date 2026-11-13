-- Copyright (c) d-SEAMS developers
-- LuaJIT FFI module generated from include/seams.h.
local ffi = require("ffi")
local C = dofile(arg[1])

local a = 2.75 * 4.0 / math.sqrt(3.0)
local fcc = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}
local basis = {{0, 0, 0}, {0.25, 0.25, 0.25}}
local xyz = ffi.new("double[24]")
local n = 0
for _, site in ipairs(fcc) do
  for _, shift in ipairs(basis) do
    xyz[3 * n] = (site[1] + shift[1]) * a
    xyz[3 * n + 1] = (site[2] + shift[2]) * a
    xyz[3 * n + 2] = (site[3] + shift[3]) * a
    n = n + 1
  end
end

local sys = C.seams_system_create()
assert(sys ~= nil)
local status = C.seams_system_set_ortho(sys, xyz, 8, 0, 0, 0, a, a, a, nil, nil)
assert(status == 0, ffi.string(C.seams_last_error()))
local labels = ffi.new("int[8]")
local counts = ffi.new("int[9]")
status = C.seams_chill_plus(sys, 3.5, 0, 0, 0, labels, counts)
assert(status == 0, ffi.string(C.seams_last_error()))
assert(counts[0] == 8)
assert(C.seams_abi() == 1)
C.seams_system_destroy(sys)
print("ffi ok")
