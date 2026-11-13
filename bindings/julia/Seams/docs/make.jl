# Build the Seams manual and run every export's doctest.
using Pkg

docsdir = @__DIR__
Pkg.activate(docsdir)
Pkg.develop(PackageSpec(path = joinpath(docsdir, "..")))
Pkg.instantiate()

root = normpath(joinpath(docsdir, "..", "..", "..", ".."))
lib = joinpath(root, "bbdir", "src", "libseams.so")
if isfile(lib)
    ENV["SEAMS_LIB"] = lib
end
engine = joinpath(root, "bbdir", "src")
if isdir(engine)
    ENV["LD_LIBRARY_PATH"] = engine * ":" * get(ENV, "LD_LIBRARY_PATH", "")
end

using Documenter
using Seams

makedocs(;
    sitename = "Seams",
    modules = [Seams],
    checkdocs = :exports,
    doctest = true,
    pages = ["index.md"],
)
