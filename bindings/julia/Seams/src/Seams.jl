# Copyright (c) d-SEAMS developers
"""
Julia front for the d-SEAMS C API. Set `SEAMS_LIB` to `libseams.so`
when it is not on the loader path.

Ice Ih (`input/traj/ice1h.lammpstrj`, every atom) is 128 hexagonal
waters and 256 primitive six-rings. The frame front asserts both:

```julia
f = read_frame(Dseams.ice_ih_path(); type=-1)
chill_plus(f; graph="cutoff").counts["hexagonal"] == 128
primitive_rings(f; graph="cutoff", max_depth=6)[6] == 256
```
"""
module Seams

export System, version, set_num_threads, num_threads, set_ortho!, chill_plus, q4q6, six_rings

const PUBLISHED = 0
const OVITO = 1
const CUTOFF = 0
const MUTUAL = 1

function _lib()
    get(ENV, "SEAMS_LIB", "libseams")
end

function _check(code::Integer)
    code == 0 && return
    err = unsafe_string(ccall((:seams_last_error, _lib()), Ptr{Cchar}, ()))
    isempty(err) && (err = unsafe_string(ccall((:seams_status_string, _lib()), Ptr{Cchar}, (Cint,), Cint(code))))
    error(err)
end

"""
    version()

Engine version string.

```jldoctest
julia> using Seams

julia> version() isa String
true
```
"""
version() = unsafe_string(ccall((:seams_version_string, _lib()), Ptr{Cchar}, ()))

"""
    set_num_threads(n)

Request `n` worker threads.

```jldoctest
julia> using Seams

julia> set_num_threads(1)

julia> num_threads() >= 1
true
```
"""
set_num_threads(n::Integer) = _check(ccall((:seams_set_num_threads, _lib()), Cint, (Cint,), Cint(n)))

"""
    num_threads()

Worker threads the engine is using.

```jldoctest
julia> using Seams

julia> num_threads() isa Int
true
```
"""
num_threads() = Int(ccall((:seams_get_num_threads, _lib()), Cint, ()))

"""
    System()

An empty system. Positions are added with [`set_ortho!`](@ref).

```jldoctest
julia> using Seams

julia> System() isa System
true
```
"""
mutable struct System
    ptr::Ptr{Cvoid}
    function System()
        ptr = ccall((:seams_system_create, _lib()), Ptr{Cvoid}, ())
        ptr == C_NULL && error("null system")
        sys = new(ptr)
        finalizer(sys) do s
            s.ptr == C_NULL && return
            ccall((:seams_system_destroy, _lib()), Cvoid, (Ptr{Cvoid},), s.ptr)
            s.ptr = C_NULL
        end
        return sys
    end
end

"""
    set_ortho!(sys, xyz, lengths)

Set an orthogonal cell. `xyz` is 3 by n.

```jldoctest
julia> using Seams

julia> sys = System();

julia> set_ortho!(sys, [0.0 1.0; 0.0 0.0; 0.0 0.0], (5.0, 5.0, 5.0)) isa System
true
```
"""
function set_ortho!(sys::System, xyz::AbstractMatrix{<:Real}, lengths::NTuple{3,Real})
    size(xyz, 1) == 3 || error("xyz is 3 by n")
    flat = collect(vec(Float64.(xyz)))
    _check(ccall((:seams_system_set_ortho, _lib()), Cint,
                 (Ptr{Cvoid}, Ptr{Float64}, Cint, Float64, Float64, Float64, Float64, Float64, Float64, Ptr{Cint}, Ptr{Cint}),
                 sys.ptr, flat, Cint(size(xyz, 2)),
                 0.0, 0.0, 0.0, Float64(lengths[1]), Float64(lengths[2]), Float64(lengths[3]),
                 C_NULL, C_NULL))
    return sys
end

"""
    chill_plus(sys; cutoff=3.5, graph=CUTOFF, rule=PUBLISHED)

CHILL+ labels and the nine-bin census. On Ice Ih the hexagonal bin
(index 2) is 128. That count is asserted by `Dseams.chill_plus`.

```jldoctest
julia> using Seams

julia> methods(chill_plus).ms[1].name == :chill_plus
true
```
"""
function chill_plus(sys::System; cutoff=3.5, graph=CUTOFF, rule=PUBLISHED)
    n = Int(ccall((:seams_system_count, _lib()), Cint, (Ptr{Cvoid},), sys.ptr))
    labels = Vector{Int32}(undef, n)
    counts = Vector{Int32}(undef, 9)
    _check(ccall((:seams_chill_plus, _lib()), Cint,
                 (Ptr{Cvoid}, Float64, Cint, Cint, Cint, Ptr{Int32}, Ptr{Int32}),
                 sys.ptr, Float64(cutoff), Cint(0), Cint(graph), Cint(rule), labels, counts))
    return labels, counts
end

"""
    q4q6(sys; cutoff=3.5)

Local q4 and q6.

```jldoctest
julia> using Seams

julia> methods(q4q6).ms[1].name == :q4q6
true
```
"""
function q4q6(sys::System; cutoff=3.5)
    n = Int(ccall((:seams_system_count, _lib()), Cint, (Ptr{Cvoid},), sys.ptr))
    q4 = Vector{Float64}(undef, n)
    q6 = Vector{Float64}(undef, n)
    _check(ccall((:seams_q4q6, _lib()), Cint,
                 (Ptr{Cvoid}, Float64, Cint, Ptr{Float64}, Ptr{Float64}),
                 sys.ptr, Float64(cutoff), Cint(0), q4, q6))
    return q4, q6
end

"""
    six_rings(sys; cutoff=3.5)

Number of primitive six-rings and their fraction. Ice Ih has 256
six-rings; `Dseams.primitive_rings` asserts that count.

```jldoctest
julia> using Seams

julia> methods(six_rings).ms[1].name == :six_rings
true
```
"""
function six_rings(sys::System; cutoff=3.5)
    rings = Ref{Ptr{Cvoid}}(C_NULL)
    _check(ccall((:seams_rings_find, _lib()), Cint,
                 (Ptr{Cvoid}, Float64, Cint, Cint, Cint, Ptr{Ptr{Cvoid}}),
                 sys.ptr, Float64(cutoff), Cint(0), Cint(CUTOFF), Cint(6), rings))
    hist = zeros(Int32, 7)
    fraction = Ref{Float64}(0.0)
    try
        _check(ccall((:seams_rings_histogram, _lib()), Cint,
                     (Ptr{Cvoid}, Ptr{Int32}, Cint), rings[], hist, Cint(7)))
        _check(ccall((:seams_rings_fraction, _lib()), Cint,
                     (Ptr{Cvoid}, Cint, Ptr{Float64}), rings[], Cint(6), fraction))
    finally
        ccall((:seams_rings_destroy, _lib()), Cvoid, (Ptr{Cvoid},), rings[])
    end
    return Int(hist[7]), fraction[]
end

end
