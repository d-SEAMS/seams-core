# Copyright (c) d-SEAMS developers
"""
    Dseams

Julia front for the d-SEAMS engine. Calls `libseams_lang` with `ccall`.
CxxWrap is not required. Set `SEAMS_LANG_LIB` to the shared library if it
is not on the loader path.
"""
module Dseams

using Libdl

export Frame, nop, nframes, read_frame, chill_plus, cages, primitive_rings, steinhardt, signature_cages, symmetry, jump_tau90, q4_q6

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
    error(err)
end

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

function nframes(path::AbstractString)
    n = ccall(dlsym(_so(), :seams_path_nframes), Cint, (Cstring,), path)
    n < 0 && _check(1)
    Int(n)
end

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
    (labels = labels, names = names, counts = histogram)
end

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
    (
        cage = per,
        # Disjoint counters, the same four numbers `seams cages` prints.
        # cage is 0 water, 1 hexagonal only, 2 double-diamond only, 3 both.
        n_ih = Int(n_hex[]),
        n_ic = Int(n_cubic[]),
        n_water = Int(n_water[]),
        n_both = Int(n_both[]),
    )
end

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
    Dict(size => Int(counts[size + 1]) for size in 3:max_depth if counts[size + 1] != 0)
end

function steinhardt(frame::Frame, order_l::Integer; cutoff::Real = 3.5, type::Integer = 0, graph::AbstractString = "mutual")
    n = nop(frame)
    ql = Vector{Cdouble}(undef, n)
    qlbar = Vector{Cdouble}(undef, n)
    code = ccall(
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
    _check(code)
    (ql = ql, qlbar = qlbar)
end

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
    (cages = Int(n_cages[]), atoms = Int(n_atoms[]))
end

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
    (q4 = q4, q6 = q6, q4bar = q4bar, q6bar = q6bar)
end

end
