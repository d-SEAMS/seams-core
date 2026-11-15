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

# GitHub's static page has no id for a #Lnn fragment, so a line link fails
# the HTML checker. A commit that has not been pushed 404s as well. The
# published file is the one on main.
struct PublishedRemote <: Documenter.Remotes.Remote end
Documenter.Remotes.repourl(::PublishedRemote) = "https://github.com/HaoZeke/seams-core"
function Documenter.Remotes.fileurl(::PublishedRemote, ::AbstractString, filename::AbstractString, linerange)
    filename = lstrip(replace(filename, '\\' => '/'), '/')
    return "https://github.com/HaoZeke/seams-core/blob/main/$(filename)"
end

makedocs(;
    sitename = "Seams",
    modules = [Seams],
    checkdocs = :exports,
    doctest = true,
    repo = PublishedRemote(),
    pages = ["index.md"],
)
