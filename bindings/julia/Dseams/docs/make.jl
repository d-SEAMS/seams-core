# Copyright (c) d-SEAMS developers
# Native reference. The Sphinx page is the how-to.
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter
using Dseams

makedocs(
    sitename = "Dseams",
    modules = [Dseams],
    format = Documenter.HTML(prettyurls = false, ansicolor = true),
    pages = ["Home" => "index.md"],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
    doctest = false,
)
