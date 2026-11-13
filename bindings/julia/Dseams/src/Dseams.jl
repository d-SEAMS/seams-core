# Copyright (c) d-SEAMS developers
"""
    Dseams

Julia front for the d-SEAMS engine. Calls `libseams_lang` with `ccall`.
CxxWrap is not required. Set `SEAMS_LANG_LIB` to the shared library if it
is not on the loader path.

Ice Ih (`input/traj/ice1h.lammpstrj`, every atom) is 128 hexagonal waters
and 256 primitive six-rings.

The C API names re-exported here are System, seams_system_count,
seams_system_read, seams_chill, seams_steinhardt, seams_cages,
seams_signature_cages, seams_neighbours_knn, seams_rdf, seams_running_cn,
seams_coordination, seams_hbonds, seams_density_axis, seams_fingerprint,
seams_contact_pairs, seams_domains, seams_prisms, seams_monolayer,
seams_rings_vertices, seams_site_table, SeamsError, ChillPlus,
RingHistogram, CageScore, SteinhardtOrder, SignatureCages, and Q4Q6.
"""
module Dseams

using Libdl

"""
    ice_ih_path()

Path of the committed Ice Ih dump, `input/traj/ice1h.lammpstrj`.
"""
ice_ih_path() = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "input", "traj", "ice1h.lammpstrj"))

export ice_ih_path
export Frame, nop, nframes, read_frame, chill_plus, cages, primitive_rings, steinhardt, signature_cages, symmetry, jump_tau90, q4_q6
export System, seams_system_count, seams_system_read, seams_chill, seams_steinhardt, seams_cages, seams_signature_cages
export seams_neighbours_knn, seams_rdf, seams_running_cn, seams_coordination, seams_hbonds
export seams_density_axis, seams_fingerprint, seams_contact_pairs, seams_domains, seams_prisms
export seams_monolayer, seams_rings_vertices, seams_site_table
export SeamsError, ChillPlus, RingHistogram, CageScore, SteinhardtOrder, SignatureCages, Q4Q6

const ICE_NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "unclassified",
    "reCubic",
    "reHex",
)

function _lib()
    get(ENV, "SEAMS_LANG_LIB", "libseams_lang")
end

# Julia 1.13 rejects a runtime library name inside ccall. The handle is
# opened once and each call goes through the symbol pointer.
const _dl = Ref{Ptr{Cvoid}}(C_NULL)

function _so()
    if _dl[] == C_NULL
        _dl[] = dlopen(_lib())
    end
    _dl[]
end

function _check(code::Integer)
    code == 0 && return
    err = unsafe_string(ccall(dlsym(_so(), :seams_last_error), Ptr{Cchar}, ()))
    isempty(err) && (err = "dseams call failed")
    throw(SeamsError(err))
end

"""
    Frame

One frame from `libseams_lang`. The finalizer calls `seams_frame_free`.
Do not use `ptr` after the frame is finalized.
"""
mutable struct Frame
    ptr::Ptr{Cvoid}
    function Frame(ptr::Ptr{Cvoid})
        ptr == C_NULL && error("null frame")
        frame = new(ptr)
        finalizer(frame) do f
            f.ptr == C_NULL && return
            ccall(dlsym(_so(), :seams_frame_free), Cvoid, (Ptr{Cvoid},), f.ptr)
            f.ptr = C_NULL
        end
        frame
    end
end

"""One frame. `ptr` is owned by the finalizer, which calls `seams_frame_free`."""
function Base.show(io::IO, frame::Frame)
    n = frame.ptr == C_NULL ? 0 : nop(frame)
    print(io, "Frame(", n, " atoms)")
end

"""CHILL+ labels, per-atom names, and the nine-label histogram."""
struct ChillPlus
    labels::Vector{Int32}
    names::Vector{String}
    counts::Dict{String,Int}
end

function Base.show(io::IO, ice::ChillPlus)
    water = get(ice.counts, "water", 0)
    print(io, "ChillPlus(", length(ice.labels), " atoms, water ", water, ")")
end

"""Hexagonal and double-diamond membership. The four counts are disjoint."""
struct CageScore
    cage::Vector{Int32}
    n_ih::Int
    n_ic::Int
    n_water::Int
    n_both::Int
end

function Base.show(io::IO, found::CageScore)
    print(io, "CageScore(hex ", found.n_ih, ", cubic ", found.n_ic, ", both ", found.n_both, ", water ", found.n_water, ")")
end

"""Primitive-ring histogram. `counts[size]` is the number of rings of that size."""
struct RingHistogram
    counts::Dict{Int,Int}
end

Base.getindex(hist::RingHistogram, size::Integer) = hist.counts[size]

function Base.show(io::IO, hist::RingHistogram)
    print(io, "RingHistogram(")
    join(io, ["$(k)=>$(v)" for (k, v) in sort!(collect(hist.counts))], ", ")
    print(io, ")")
end

"""Local Steinhardt `ql` and the Lechner-Dellago `qlbar`."""
struct SteinhardtOrder
    ql::Vector{Float64}
    qlbar::Vector{Float64}
end

function Base.show(io::IO, order::SteinhardtOrder)
    print(io, "SteinhardtOrder(n=", length(order.ql), ")")
end

"""Closed cages of one ring census."""
struct SignatureCages
    cages::Int
    atoms::Int
end

function Base.show(io::IO, sig::SignatureCages)
    print(io, "SignatureCages(", sig.cages, " cages, ", sig.atoms, " atoms)")
end

"""Local and averaged q4 and q6."""
struct Q4Q6
    q4::Vector{Float64}
    q6::Vector{Float64}
    q4bar::Vector{Float64}
    q6bar::Vector{Float64}
end

function Base.show(io::IO, order::Q4Q6)
    print(io, "Q4Q6(n=", length(order.q6), ")")
end

"""
    nframes(path)

Number of frames in a trajectory. An XYZ file counts as one.
"""
function nframes(path::AbstractString)
    n = ccall(dlsym(_so(), :seams_path_nframes), Cint, (Cstring,), path)
    n < 0 && _check(1)
    Int(n)
end

"""
    read_frame(path; frame=1, type=0)

Read one frame. `frame` is 1-based. `type = 0` guesses oxygen. A positive
type keeps that species. A negative type keeps every atom.
"""
function read_frame(path::AbstractString; frame::Integer = 1, type::Integer = 0)
    ptr = ccall(
        dlsym(_so(), :seams_frame_read),
        Ptr{Cvoid},
        (Cstring, Cint, Cint),
        path,
        frame,
        type,
    )
    ptr == C_NULL && _check(1)
    Frame(ptr)
end

"""
    nop(frame)

Number of atoms in `frame`.
"""
nop(frame::Frame) = Int(ccall(dlsym(_so(), :seams_frame_nop), Cint, (Ptr{Cvoid},), frame.ptr))

"""
    chill_plus(frame; cutoff = 3.5, type = 0, graph = "cutoff")

CHILL+ labels. `graph = "cutoff"` (the default) bonds every pair closer
than `cutoff`. `graph = "mutual"` keeps the mutual four nearest
neighbours inside `cutoff + 1.5`.
"""
function chill_plus(frame::Frame; cutoff::Real = 3.5, type::Integer = 0, graph::AbstractString = "cutoff")
    n = nop(frame)
    labels = Vector{Cint}(undef, n)
    counts = Vector{Cint}(undef, 9)
    code = ccall(
        dlsym(_so(), :seams_frame_chill_plus),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Cstring, Ptr{Cint}, Ptr{Cint}),
        frame.ptr,
        cutoff,
        type,
        graph,
        labels,
        counts,
    )
    _check(code)
    names = [ICE_NAMES[Int(label) + 1] for label in labels]
    histogram = Dict(ICE_NAMES[i] => Int(counts[i]) for i in 1:9)
    ChillPlus(labels, names, histogram)
end

"""
    chill_plus(path; frame=1, type=0, cutoff=3.5, graph="cutoff")

Read `path` and label it. `type` is the trajectory filter, not the
CHILL+ species filter.
"""
function chill_plus(path::AbstractString; frame::Integer=1, type::Integer=0, cutoff::Real=3.5, graph::AbstractString="cutoff")
    chill_plus(read_frame(path; frame=frame, type=type); cutoff=cutoff, graph=graph)
end

"""
    cages(frame; cutoff=3.5, type=0, k=4, graph="seeded")

Hexagonal and double-diamond membership. `cage` is 0 water, 1 hexagonal
only, 2 double-diamond only, 3 both. The four counts are disjoint and
match `seams cages`. `cutoff` is angstrom.
"""
function cages(frame::Frame; cutoff::Real = 3.5, type::Integer = 0, k::Integer = 4, graph::AbstractString = "seeded")
    n = nop(frame)
    per = Vector{Cint}(undef, n)
    n_hex = Ref{Cint}(0)
    n_cubic = Ref{Cint}(0)
    n_water = Ref{Cint}(0)
    n_both = Ref{Cint}(0)
    code = ccall(
        dlsym(_so(), :seams_frame_cages),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Cint, Cstring, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        frame.ptr,
        cutoff,
        type,
        k,
        graph,
        0,
        per,
        n_hex,
        n_cubic,
        n_water,
        n_both,
    )
    _check(code)
    # Disjoint counters, the same four numbers `seams cages` prints.
    # cage is 0 water, 1 hexagonal only, 2 double-diamond only, 3 both.
    CageScore(per, Int(n_hex[]), Int(n_cubic[]), Int(n_water[]), Int(n_both[]))
end

"""
    primitive_rings(frame; cutoff=3.5, type=0, graph="mutual", max_depth=6)

Histogram of primitive rings. Index with the ring size: `rings[6]`.
"""
function primitive_rings(frame::Frame; cutoff::Real = 3.5, type::Integer = 0, graph::AbstractString = "mutual", max_depth::Integer = 6)
    counts = zeros(Cint, max_depth + 1)
    code = ccall(
        dlsym(_so(), :seams_frame_rings),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Cstring, Cint, Ptr{Cint}, Cint),
        frame.ptr,
        cutoff,
        type,
        graph,
        max_depth,
        counts,
        length(counts),
    )
    _check(code)
    RingHistogram(Dict(size => Int(counts[size + 1]) for size in 3:max_depth if counts[size + 1] != 0))
end

"""
    steinhardt(frame, order_l; cutoff=3.5, type=0, graph=nothing)

Local Steinhardt parameter. `order_l` is 3, 4, 6, 8, or 12. Omitting
`graph` uses the cutoff shell. `"mutual"` is the four nearest neighbours
inside `cutoff + 1.5`.
"""
function steinhardt(frame::Frame, order_l::Integer; cutoff::Real = 3.5, type::Integer = 0, graph::Union{Nothing,AbstractString} = nothing)
    n = nop(frame)
    ql = Vector{Cdouble}(undef, n)
    qlbar = Vector{Cdouble}(undef, n)
    # No graph argument is the cutoff shell. A name selects that graph.
    code = if graph === nothing
        ccall(
            dlsym(_so(), :seams_frame_steinhardt),
            Cint,
            (Ptr{Cvoid}, Cdouble, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
            frame.ptr,
            cutoff,
            type,
            order_l,
            ql,
            qlbar,
        )
    else
        ccall(
            dlsym(_so(), :seams_frame_steinhardt_graph),
            Cint,
            (Ptr{Cvoid}, Cdouble, Cint, Cint, Cstring, Ptr{Cdouble}, Ptr{Cdouble}),
            frame.ptr,
            cutoff,
            type,
            order_l,
            graph,
            ql,
            qlbar,
        )
    end
    _check(code)
    SteinhardtOrder(ql, qlbar)
end

"""
    signature_cages(frame, spec; cutoff=3.5, type=0, graph="cutoff")

Closed cages whose ring census matches `spec`, for example `"512"`.
"""
function signature_cages(frame::Frame, spec::AbstractString; cutoff::Real = 3.5, type::Integer = 0, graph::AbstractString = "cutoff")
    n_cages = Ref{Cint}(0)
    n_atoms = Ref{Cint}(0)
    code = ccall(
        dlsym(_so(), :seams_frame_signature),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Cstring, Cstring, Ptr{Cint}, Ptr{Cint}),
        frame.ptr,
        cutoff,
        type,
        graph,
        spec,
        n_cages,
        n_atoms,
    )
    _check(code)
    SignatureCages(Int(n_cages[]), Int(n_atoms[]))
end

"""
    jump_tau90(first, second, dt=1.0; oxygen_type=1, hydrogen_type=2)

Time for an H-H vector to rotate 90 degrees between two frames.
Returns `dt` when any water rotates that far, and `NaN` otherwise.
"""
function jump_tau90(first::Frame, second::Frame, dt::Real = 1.0; oxygen_type::Integer = 1, hydrogen_type::Integer = 2)
    tau = Ref{Cdouble}(0.0)
    code = ccall(
        dlsym(_so(), :seams_frames_jump_tau90),
        Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble, Cint, Cint, Ptr{Cdouble}),
        first.ptr,
        second.ptr,
        dt,
        oxygen_type,
        hydrogen_type,
        tau,
    )
    _check(code)
    tau[]
end

"""
    symmetry(frame; shell_rmsd=0.12)

Per-atom shell label. 0 is ice Ih, 1 is ice Ic, 2 is stacking disorder,
and -1 is anything else.
"""
function symmetry(frame::Frame; shell_rmsd::Real = 0.12)
    n = nop(frame)
    labels = Vector{Cint}(undef, n)
    code = ccall(
        dlsym(_so(), :seams_frame_symmetry),
        Cint,
        (Ptr{Cvoid}, Cdouble, Ptr{Cint}),
        frame.ptr,
        shell_rmsd,
        labels,
    )
    _check(code)
    labels
end

"""
    q4_q6(frame; cutoff=3.5, type=0)

Local and Lechner-Dellago q4 and q6 on the cutoff shell. `cutoff` is angstrom.
"""
function q4_q6(frame::Frame; cutoff::Real = 3.5, type::Integer = 0)
    n = nop(frame)
    q4 = Vector{Cdouble}(undef, n)
    q6 = Vector{Cdouble}(undef, n)
    q4bar = Vector{Cdouble}(undef, n)
    q6bar = Vector{Cdouble}(undef, n)
    code = ccall(
        dlsym(_so(), :seams_frame_q4_q6),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        frame.ptr,
        cutoff,
        type,
        q4,
        q6,
        q4bar,
        q6bar,
    )
    _check(code)
    Q4Q6(q4, q6, q4bar, q6bar)
end

include(joinpath(@__DIR__, "..", "..", "capi.jl"))

end
