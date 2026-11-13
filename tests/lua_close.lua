-- Copyright (c) d-SEAMS developers
-- Lua 5.4 to-be-closed variable. __close and __gc both release the handle.
local core = require("dseams_core")

local function scoped()
  local held <close> = core.system()
  held:set_ortho({0, 0, 0}, 0, 0, 0, 4, 4, 4)
  assert(held:count() == 1)
  held:close()
end

scoped()
print("close ok")
