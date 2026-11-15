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

root = normpath(joinpath(docsdir, "..", "..", "..", ".."))
lang = joinpath(root, "bbdir", "bindings", "libseams_lang.so")
if isfile(lang)
    ENV["SEAMS_LANG_LIB"] = lang
    ENV["LD_LIBRARY_PATH"] = joinpath(root, "bbdir", "bindings") * ":" * joinpath(root, "bbdir", "src")
end

DocMeta.setdocmeta!(Dseams, :DocTestSetup, :(using Dseams); recursive = true)

# Same remote as the Seams manual: the file on main, with no line fragment.
struct PublishedRemote <: Documenter.Remotes.Remote end
Documenter.Remotes.repourl(::PublishedRemote) = "https://github.com/d-SEAMS/seams-core"
function Documenter.Remotes.fileurl(::PublishedRemote, ::AbstractString, filename::AbstractString, linerange)
    filename = lstrip(replace(filename, '\\' => '/'), '/')
    return "https://github.com/d-SEAMS/seams-core/blob/main/$(filename)"
end

makedocs(;
    sitename = "Dseams",
    modules = [Dseams],
    checkdocs = :exports,
    doctest = true,
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "") == "true"),
    repo = PublishedRemote(),
    pages = ["Home" => "index.md", "API" => "api.md"],
)
