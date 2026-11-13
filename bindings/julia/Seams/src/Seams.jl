# Copyright (c) d-SEAMS developers
"""
    Seams

Julia front for the d-SEAMS C API (`include/seams.h`).

Set `SEAMS_LIB` to `libseams.so` when it is not on the loader path.
Handles are mutable structs with finalizers. A failing call throws
[`SeamsError`](@ref). The message is `seams_handle_error`.

`System`, `Neighbours`, `Rings`, and `CageSet` from the included C layer
own their pointers. Result structs iterate in field order, so
`labels, counts = chill_plus(sys)` still unpacks.
"""
module Seams

export System, version, set_num_threads, num_threads, set_ortho!, chill_plus, q4q6, six_rings
export seams_system_count, seams_system_read, seams_chill, seams_steinhardt, seams_cages, seams_signature_cages
export seams_neighbours_knn, seams_rdf, seams_running_cn, seams_coordination, seams_hbonds
export seams_density_axis, seams_fingerprint, seams_contact_pairs, seams_domains, seams_prisms
export seams_monolayer, seams_rings_vertices, seams_site_table
export SeamsError, ChillLabels, Q4Q6Pair, SixRings

include(joinpath(@__DIR__, "..", "..", "capi.jl"))

"""Published CHILL+ rule. Matches `SEAMS_CHILL_PUBLISHED`."""
const PUBLISHED = 0

"""OVITO CHILL+ rule. Matches `SEAMS_CHILL_OVITO`."""
const OVITO = 1

"""Cutoff neighbour graph. Matches `SEAMS_GRAPH_CUTOFF`."""
const CUTOFF = 0

"""Mutual four-nearest neighbour graph. Matches `SEAMS_GRAPH_MUTUAL`."""
const MUTUAL = 1

"""
    version()

Version text from `seams_version_string`. The caller does not free it.
"""
version() = unsafe_string(ccall(_seams_sym(:seams_version_string), Ptr{Cchar}, ()))

"""
    set_num_threads(n)

OpenMP width for later analysis. `n = 0` restores the saved default.
A negative width throws [`SeamsError`](@ref). The width is process-global.
"""
function set_num_threads(n::Integer)
    _seams_check(ccall(_seams_sym(:seams_set_num_threads), Cint, (Cint,), Cint(n)))
    return nothing
end

"""
    num_threads()

Width the next analysis call will request.
"""
num_threads() = Int(ccall(_seams_sym(:seams_get_num_threads), Cint, ()))

"""
    set_ortho!(sys, xyz, lengths)

Store a `3 × n` coordinate matrix in an orthorhombic cell at the origin.
`lengths` is `(lx, ly, lz)` in angstrom. The library copies the arrays.
"""
function set_ortho!(sys::System, xyz::AbstractMatrix{<:Real}, lengths::NTuple{3,Real})
    size(xyz, 1) == 3 || throw(SeamsError("xyz is 3 by n"))
    flat = collect(vec(Float64.(xyz)))
    _seams_check(ccall(_seams_sym(:seams_system_set_ortho), Cint,
                       (Ptr{Cvoid}, Ptr{Float64}, Cint, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cint}, Ptr{Cint}),
                       sys.ptr, flat, Cint(size(xyz, 2)),
                       0.0, 0.0, 0.0, Float64(lengths[1]), Float64(lengths[2]), Float64(lengths[3]),
                       C_NULL, C_NULL), sys.ptr)
    return sys
end

"""CHILL+ labels and the nine population counts. Iterates as `(labels, counts)`."""
struct ChillLabels
    labels::Vector{Int32}
    counts::Vector{Int32}
end

Base.length(result::ChillLabels) = nfields(result)
Base.iterate(result::ChillLabels, i::Int=1) = i > nfields(result) ? nothing : (getfield(result, i), i + 1)

function Base.show(io::IO, result::ChillLabels)
    print(io, "ChillLabels(", length(result.labels), " atoms, cubic ", Int(result.counts[1]), ")")
end

"""
    chill_plus(sys; cutoff=3.5, graph=CUTOFF, rule=PUBLISHED)
    chill_plus(xyz, lengths; kwargs...)

CHILL+ on a [`System`](@ref), or on a `3 × n` matrix in an orthorhombic cell.
`cutoff` is angstrom.
"""
function chill_plus(sys::System; cutoff=3.5, graph=CUTOFF, rule=PUBLISHED)
    n = Int(ccall(_seams_sym(:seams_system_count), Cint, (Ptr{Cvoid},), sys.ptr))
    labels = Vector{Int32}(undef, n)
    counts = Vector{Int32}(undef, 9)
    _seams_check(ccall(_seams_sym(:seams_chill_plus), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Int32}, Ptr{Int32}),
                       sys.ptr, Float64(cutoff), Cint(0), Cint(graph), Cint(rule), labels, counts), sys.ptr)
    return ChillLabels(labels, counts)
end

function chill_plus(xyz::AbstractMatrix{<:Real}, lengths::NTuple{3,Real}; kwargs...)
    sys = System()
    set_ortho!(sys, xyz, lengths)
    return chill_plus(sys; kwargs...)
end

"""Local q4 and q6. Iterates as `(q4, q6)`."""
struct Q4Q6Pair
    q4::Vector{Float64}
    q6::Vector{Float64}
end

Base.length(result::Q4Q6Pair) = nfields(result)
Base.iterate(result::Q4Q6Pair, i::Int=1) = i > nfields(result) ? nothing : (getfield(result, i), i + 1)

function Base.show(io::IO, result::Q4Q6Pair)
    print(io, "Q4Q6Pair(n=", length(result.q6), ")")
end

"""
    q4q6(sys; cutoff=3.5)

Local q4 and q6 on the cutoff shell. `cutoff` is angstrom.
"""
function q4q6(sys::System; cutoff=3.5)
    n = Int(ccall(_seams_sym(:seams_system_count), Cint, (Ptr{Cvoid},), sys.ptr))
    q4 = Vector{Float64}(undef, n)
    q6 = Vector{Float64}(undef, n)
    _seams_check(ccall(_seams_sym(:seams_q4q6), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Ptr{Float64}, Ptr{Float64}),
                       sys.ptr, Float64(cutoff), Cint(0), q4, q6), sys.ptr)
    return Q4Q6Pair(q4, q6)
end

"""Six-ring count and the fraction of atoms on a six-ring. Iterates as `(count, fraction)`."""
struct SixRings
    count::Int
    fraction::Float64
end

Base.length(result::SixRings) = nfields(result)
Base.iterate(result::SixRings, i::Int=1) = i > nfields(result) ? nothing : (getfield(result, i), i + 1)

function Base.show(io::IO, result::SixRings)
    print(io, "SixRings(", result.count, ", ", result.fraction, ")")
end

"""
    six_rings(sys; cutoff=3.5)

Primitive six-ring count and the fraction of atoms on a six-ring.
The ring handle is freed before this function returns.
"""
function six_rings(sys::System; cutoff=3.5)
    rings = Ref{Ptr{Cvoid}}(C_NULL)
    _seams_check(ccall(_seams_sym(:seams_rings_find), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Ptr{Cvoid}}),
                       sys.ptr, Float64(cutoff), Cint(0), Cint(CUTOFF), Cint(6), rings), sys.ptr)
    hist = zeros(Int32, 7)
    fraction = Ref{Float64}(0.0)
    try
        _seams_check(ccall(_seams_sym(:seams_rings_histogram), Cint,
                           (Ptr{Cvoid}, Ptr{Int32}, Cint), rings[], hist, Cint(7)), rings[])
        _seams_check(ccall(_seams_sym(:seams_rings_fraction), Cint,
                           (Ptr{Cvoid}, Cint, Ptr{Float64}), rings[], Cint(6), fraction), rings[])
    finally
        ccall(_seams_sym(:seams_rings_destroy), Cvoid, (Ptr{Cvoid},), rings[])
    end
    return SixRings(Int(hist[7]), fraction[])
end

end
