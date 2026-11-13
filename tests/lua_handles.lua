-- Copyright (c) d-SEAMS developers
-- Userdata lifetime and the C error string on a rejected box.
local dseams = require("dseams")
local core = dseams.core

local sys = core.system()
local ok, err = pcall(function()
  sys:set_ortho({0, 0, 0}, 0, 0, 0, 0, 1, 1)
end)
assert(not ok)
assert(tostring(err):find("box edge lengths must be positive", 1, true), tostring(err))

sys:close()
local closed, closed_err = pcall(function()
  sys:count()
end)
assert(not closed)
assert(tostring(closed_err):find("system is closed", 1, true), tostring(closed_err))

local dropped = core.system()
dropped:set_ortho({0, 0, 0, 1, 0, 0}, 0, 0, 0, 4, 4, 4)
dropped = nil
collectgarbage()
collectgarbage()

local yoda = require("yoda")
assert(yoda.version() == dseams.version())
assert(yoda.read == dseams.read)
print("handles ok")
