-- Copyright (c) d-SEAMS developers
-- Times the Lua module over the C API. Commands are "threads op path".

local core = require("dseams_core")

local frames = {}
local sink = 0

local function retune(threads)
  -- seams_set_num_threads updates the OpenMP team for the next analysis.
  core.set_num_threads(threads)
end

local function system_for(path)
  local held = frames[path]
  if held == nil then
    held = {sys = core.system(), loaded = false}
    frames[path] = held
  end
  return held
end

local function ensure(path)
  local held = system_for(path)
  if not held.loaded then
    held.sys:read(path, 1, -1)
    held.loaded = true
  end
  return held
end

while true do
  local line = io.read()
  if line == nil or line:sub(1, 4) == "quit" then
    break
  end
  if line == "" then
    goto next_line
  end
  local threads, op, path = line:match("^(%d+)%s+(%S+)%s+(%S+)$")
  if threads == nil then
    io.stderr:write("error command\n")
    os.exit(1)
  end
  retune(tonumber(threads))
  local t0 = core.monotonic()
  if op == "read" then
    local held = system_for(path)
    local count = held.sys:read(path, 1, -1)
    held.loaded = true
    sink = sink + count
  elseif op == "knn" then
    local held = ensure(path)
    local nb = held.sys:knn(4, 3.5, 0, 1)
    sink = sink + nb:count()
    nb:close()
  elseif op == "steinhardt" then
    local held = ensure(path)
    local ql, bar = held.sys:steinhardt(3.5, 0, 6, 1)
    sink = sink + ql[1] + bar[1]
  elseif op == "chill" then
    local held = ensure(path)
    local labels, counts = held.sys:chill(3.5, 0, 0, 0)
    sink = sink + counts[2] + labels[1]
  elseif op == "rings" then
    local held = ensure(path)
    local rings = held.sys:rings(3.5, 0, 0, 6)
    sink = sink + rings:count()
    rings:close()
  else
    io.stderr:write("error " .. op .. "\n")
    os.exit(1)
  end
  local ms = (core.monotonic() - t0) * 1e3
  io.write(string.format("%.9f\n", ms))
  io.flush()
  ::next_line::
end
io.stderr:write(string.format("sink %.6f\n", sink))
