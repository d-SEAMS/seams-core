-- Copyright (c) d-SEAMS developers
-- Compile one Fennel file with the compiler named by arg[1].
local fennel = dofile(arg[1])
local source = assert(io.open(arg[2]))
local text = source:read("*a")
source:close()
local compiled = fennel.compileString(text, {filename = arg[2]})
io.write(compiled)
