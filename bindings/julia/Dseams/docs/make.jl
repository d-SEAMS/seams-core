# Build the Dseams manual and run every export's doctest.
using Pkg

docsdir = @__DIR__
Pkg.activate(docsdir)
Pkg.develop(PackageSpec(path = joinpath(docsdir, "..")))
Pkg.instantiate()

root = normpath(joinpath(docsdir, "..", "..", "..", ".."))
lang = joinpath(root, "bbdir", "bindings", "libseams_lang.so")
if isfile(lang)
    ENV["SEAMS_LANG_LIB"] = lang
end
engine = joinpath(root, "bbdir", "src")
if isdir(engine)
    ENV["LD_LIBRARY_PATH"] = engine * ":" * joinpath(root, "bbdir", "bindings") *
                             ":" * get(ENV, "LD_LIBRARY_PATH", "")
end

using Documenter
using Dseams

makedocs(;
    sitename = "Dseams",
    modules = [Dseams],
    checkdocs = :exports,
    doctest = true,
    pages = ["index.md"],
)
