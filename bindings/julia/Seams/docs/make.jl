# Copyright (c) d-SEAMS developers
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Documenter
using Seams

makedocs(
    sitename = "Seams",
    modules = [Seams],
    format = Documenter.HTML(prettyurls = false, ansicolor = true),
    pages = ["Home" => "index.md"],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
    doctest = false,
)
