"""
    Dseams

Julia front for the d-SEAMS engine. Calls `libseams_lang` with `ccall`.
CxxWrap is not required. Set `SEAMS_LANG_LIB` to the shared library if it
is not on the loader path.
"""
module Dseams

export Frame, nop, nframes, read_frame, chill_plus, cages, primitive_rings, steinhardt

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

function _check(code::Integer)
    code == 0 && return
    err = unsafe_string(ccall((:seams_last_error, _lib()), Ptr{Cchar}, ()))
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
            ccall((:seams_frame_free, _lib()), Cvoid, (Ptr{Cvoid},), f.ptr)
            f.ptr = C_NULL
        end
        frame
    end
end

function nframes(path::AbstractString)
    n = ccall((:seams_path_nframes, _lib()), Cint, (Cstring,), path)
    n < 0 && _check(1)
    Int(n)
end

function read_frame(path::AbstractString; frame::Integer = 1, type::Integer = 0)
    ptr = ccall(
        (:seams_frame_read, _lib()),
        Ptr{Cvoid},
        (Cstring, Cint, Cint),
        path,
        frame,
        type,
    )
    ptr == C_NULL && _check(1)
    Frame(ptr)
end

nop(frame::Frame) = Int(ccall((:seams_frame_nop, _lib()), Cint, (Ptr{Cvoid},), frame.ptr))

"""
    chill_plus(frame; cutoff = 3.5, type = 0, graph = "cutoff", rule = "chill+")

CHILL+ labels. `rule = "chill+"` is the published shell gate, four-nearest
q3, and the Table 1 neighbour test. `rule = "ovito-chill+"` averages q3
over the cutoff shell and skips the neighbour test.
`graph = "cutoff"` (the default) bonds every pair closer
than `cutoff`. `graph = "mutual"` keeps the mutual four nearest
neighbours inside `cutoff + 1.5`.
"""
function chill_plus(frame::Frame; cutoff::Real = 3.5, type::Integer = 0, graph::AbstractString = "cutoff", rule::AbstractString = "chill+")
    n = nop(frame)
    labels = Vector{Cint}(undef, n)
    counts = Vector{Cint}(undef, 9)
    code = ccall(
        (:seams_frame_chill_plus_rule, _lib()),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Cstring, Cstring, Ptr{Cint}, Ptr{Cint}),
        frame.ptr,
        cutoff,
        type,
        graph,
        rule,
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
        (:seams_frame_cages, _lib()),
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
        (:seams_frame_rings, _lib()),
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

function steinhardt(frame::Frame, order_l::Integer; cutoff::Real = 3.5, type::Integer = 0)
    n = nop(frame)
    ql = Vector{Cdouble}(undef, n)
    qlbar = Vector{Cdouble}(undef, n)
    code = ccall(
        (:seams_frame_steinhardt, _lib()),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
        frame.ptr,
        cutoff,
        type,
        order_l,
        ql,
        qlbar,
    )
    _check(code)
    (ql = ql, qlbar = qlbar)
end

end
