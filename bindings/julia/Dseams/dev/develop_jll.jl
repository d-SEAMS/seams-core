# Make libseams_jll resolvable in the active project before it is in the
# General registry, then instantiate. A project that depends on Dseams
# (docs/, for one) gets this checkout of Dseams too.
#
#   julia --project=bindings/julia/Dseams develop_jll.jl
#       the stand-in in dev/libseams_jll, which has no binaries
#   julia --project=ENV develop_jll.jl PATH
#       a JLL from `julia build_tarballs.jl --deploy=local`
#       (~/.julia/dev/libseams_jll)
using Pkg

path = isempty(ARGS) ? joinpath(@__DIR__, "libseams_jll") : abspath(expanduser(ARGS[1]))
paths = [path]
project = Pkg.project()
if project.name != "Dseams" && haskey(project.dependencies, "Dseams")
    push!(paths, dirname(@__DIR__))
end

# Julia 1.11 and later write developed paths under [sources] in the project
# file. The projects in this checkout (Dseams, docs) are tracked, so they get
# relative paths, and the [sources] table is taken out again afterwards; the
# gitignored Manifest keeps the developed copies.
tracked = startswith(project.path, joinpath(dirname(@__DIR__), ""))
projdir = dirname(project.path)
cd(projdir) do
    specs = [PackageSpec(; path = tracked ? relpath(p, projdir) : p) for p in paths]
    Pkg.develop(specs)
    Pkg.instantiate()
end
if tracked
    text = read(project.path, String)
    stripped = replace(text, r"^\[sources\]\n(?:[^\[\n].*\n|\n)*"m => "")
    stripped == text || write(project.path, stripped)
end
