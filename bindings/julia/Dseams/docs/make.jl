# Build the Dseams manual and run every doctest.
#
#   julia --project=docs docs/make.jl
#
# Before libseams_jll is registered, first run dev/develop_jll.jl in the
# docs project (see its header).
using Pkg

docsdir = @__DIR__
Pkg.activate(docsdir)
Pkg.develop(PackageSpec(path = joinpath(docsdir, "..")))
Pkg.instantiate()

using Documenter
using Dseams

DocMeta.setdocmeta!(Dseams, :DocTestSetup, :(using Dseams); recursive = true)

makedocs(;
    sitename = "Dseams",
    modules = [Dseams],
    checkdocs = :exports,
    doctest = true,
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "") == "true"),
    pages = ["Home" => "index.md", "API" => "api.md"],
)
