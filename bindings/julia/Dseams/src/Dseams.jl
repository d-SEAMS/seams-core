# Copyright (c) d-SEAMS developers
"""
    Dseams

Julia front for the d-SEAMS engine. Calls `libseams_lang` with `ccall`.
CxxWrap is not required. Set `SEAMS_LANG_LIB` to the shared library if it
is not on the loader path.
"""
module Dseams

using Libdl

"""
    ice_ih_path()

Path of the committed Ice Ih dump, `input/traj/ice1h.lammpstrj`.
"""
ice_ih_path() = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "input", "traj", "ice1h.lammpstrj"))


export Frame, nop, nframes, read_frame, neighbour_pairs, chill_plus, cages, primitive_rings, steinhardt, signature_cages, symmetry, jump_tau90, q4_q6

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

"""
    Frame(ptr)

One trajectory frame. Freed when the object is collected.

```jldoctest
julia> using Dseams

julia> read_frame(Dseams.ice_ih_path(); type=-1) isa Frame
true
```
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

"""
    nframes(path)

Number of frames in a trajectory.

```jldoctest
julia> using Dseams

julia> nframes(Dseams.ice_ih_path())
1
```
"""
function nframes(path::AbstractString)
    n = ccall(dlsym(_so(), :seams_path_nframes), Cint, (Cstring,), path)
    n < 0 && _check(1)
    Int(n)
end

"""
    read_frame(path; frame=1, type=0)

Read one frame. `type = -1` keeps every atom.

```jldoctest
julia> using Dseams

julia> nop(read_frame(Dseams.ice_ih_path(); type=-1))
128
```
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

Number of atoms in a frame.

```jldoctest
julia> using Dseams

julia> nop(read_frame(Dseams.ice_ih_path(); type=-1))
128
```
"""
nop(frame::Frame) = Int(ccall(dlsym(_so(), :seams_frame_nop), Cint, (Ptr{Cvoid},), frame.ptr))

"""
    chill_plus(frame; cutoff = 3.5, type = 0, graph = "cutoff")

CHILL+ labels. `graph = "cutoff"` (the default) bonds every pair closer
than `cutoff`. `graph = "mutual"` keeps the mutual four nearest
neighbours inside `cutoff + 1.5`.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> chill_plus(f; graph="cutoff").counts["hexagonal"]
128
```
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

"""
    cages(frame; cutoff=3.5, type=0, k=4, graph="seeded")

Hexagonal and cubic cage membership.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> cages(f; graph="cutoff").n_ih
128
```
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

"""
    neighbour_pairs(frame; cutoff=3.5, type=0, graph="cutoff")

Directed pairs in the neighbour list. The self entry is omitted.
"""
function neighbour_pairs(frame::Frame; cutoff::Real = 3.5, type::Integer = 0, graph::AbstractString = "cutoff")
    code = ccall(
        dlsym(_so(), :seams_frame_neighbours),
        Cint,
        (Ptr{Cvoid}, Cdouble, Cint, Cstring),
        frame.ptr,
        Cdouble(cutoff),
        Cint(type),
        graph,
    )
    _check(code)
    Int(ccall(dlsym(_so(), :seams_frame_neighbour_pairs), Cint, (Ptr{Cvoid},), frame.ptr))
end

"""
    primitive_rings(frame; cutoff=3.5, type=0, graph="mutual", max_depth=6)

Primitive-ring histogram. The key is the ring size.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> primitive_rings(f; graph="cutoff", max_depth=6)[6]
256
```
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
    Dict(size => Int(counts[size + 1]) for size in 3:max_depth if counts[size + 1] != 0)
end

"""
    steinhardt(frame, order_l; cutoff=3.5, type=0, graph=nothing)

Local Steinhardt parameter and the Lechner-Dellago average.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> length(steinhardt(f, 6).ql)
128
```
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
    (ql = ql, qlbar = qlbar)
end

"""
    signature_cages(frame, spec; cutoff=3.5, type=0, graph="cutoff")

Closed cages whose faces match `spec`.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> signature_cages(f, "512").cages isa Int
true
```
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
    (cages = Int(n_cages[]), atoms = Int(n_atoms[]))
end

"""
    jump_tau90(first, second, dt=1.0; oxygen_type=1, hydrogen_type=2)

Jump-rotor tau_90 between two frames.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> jump_tau90(f, f) isa Float64
true
```
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

Per-atom shell label. 0 is ice Ih.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> length(symmetry(f))
128
```
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

Local q4 and q6 and the Lechner-Dellago averages.

```jldoctest
julia> using Dseams

julia> f = read_frame(Dseams.ice_ih_path(); type=-1);

julia> length(q4_q6(f).q6)
128
```
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
    (q4 = q4, q6 = q6, q4bar = q4bar, q6bar = q6bar)
end

end
