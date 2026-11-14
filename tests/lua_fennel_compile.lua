-- Copyright (c) d-SEAMS developers
-- Compile one Fennel file with the compiler named by arg[1].
local compiler, source_path = arg[1], arg[2]
-- fennel.lua is a CLI after its preload table. A script arg makes that CLI
-- compile the compiler. Clear arg, load the preload, then require the module.
local saved = arg
-- "--compile" with no files installs the module and leaves the CLI.
arg = {"--compile"}
local returned = dofile(compiler)
arg = saved
-- share/lua/5.4/fennel.lua, which meson also accepts, is the library
-- form: it returns the module instead of preloading it.
local found, fennel = pcall(require, "fennel")
if not found then
  assert(type(returned) == "table", fennel)
  fennel = returned
end
local source = assert(io.open(source_path))
local text = source:read("*a")
source:close()
local compiled = fennel.compileString(text, {filename = source_path})
io.write(compiled)
