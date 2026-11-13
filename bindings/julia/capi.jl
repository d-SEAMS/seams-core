# Copyright (c) d-SEAMS developers
# ccall layer over include/seams.h. Included inside Seams and Dseams.

using Libdl

"""A failed `libseams` call. `message` is the text from `seams_handle_error`."""
struct SeamsError <: Exception
    message::String
end

Base.showerror(io::IO, err::SeamsError) = print(io, err.message)

const _seams_dl = Ref{Ptr{Cvoid}}(C_NULL)

function _seams_so()
    if _seams_dl[] == C_NULL
        _seams_dl[] = dlopen(get(ENV, "SEAMS_LIB", "libseams"))
    end
    _seams_dl[]
end

_seams_sym(name::Symbol) = dlsym(_seams_so(), name)

function _seams_check(code::Integer, handle::Ptr{Cvoid}=C_NULL)
    code == 0 && return
    err = unsafe_string(ccall(_seams_sym(:seams_handle_error), Ptr{Cchar}, (Ptr{Cvoid},), handle))
    if isempty(err)
        err = unsafe_string(ccall(_seams_sym(:seams_status_string), Ptr{Cchar}, (Cint,), Cint(code)))
    end
    throw(SeamsError(isempty(err) ? "dseams call failed" : err))
end

"""
    System()

One frame from `seams_system_create`. The finalizer calls `seams_system_destroy`.
Neighbour, ring, and cage handles copied out of it stay valid after `System` is freed.
"""
mutable struct System
    ptr::Ptr{Cvoid}
    function System()
        ptr = ccall(_seams_sym(:seams_system_create), Ptr{Cvoid}, ())
        ptr == C_NULL && throw(SeamsError("null system"))
        sys = new(ptr)
        finalizer(sys) do s
            s.ptr == C_NULL && return
            ccall(_seams_sym(:seams_system_destroy), Cvoid, (Ptr{Cvoid},), s.ptr)
            s.ptr = C_NULL
        end
        return sys
    end
end

function Base.show(io::IO, sys::System)
    n = sys.ptr == C_NULL ? 0 : seams_system_count(sys)
    print(io, "System(", n, " atoms)")
end

"""
    seams_system_count(sys)

Atom count stored on `sys`, or 0 when the handle is empty.
"""
function seams_system_count(sys::System)
    Int(ccall(_seams_sym(:seams_system_count), Cint, (Ptr{Cvoid},), sys.ptr))
end

function _count(sys::System)
    seams_system_count(sys)
end

"""
    seams_system_read(sys, path, frame=1, type_filter=-1)

Read one frame into `sys`. `frame` is 1-based. A negative filter keeps every atom.
The previous frame is kept when the call throws [`SeamsError`](@ref).
"""
function seams_system_read(sys::System, path::AbstractString, frame::Integer=1, type_filter::Integer=-1)
    _seams_check(ccall(_seams_sym(:seams_system_read), Cint,
                       (Ptr{Cvoid}, Cstring, Cint, Cint),
                       sys.ptr, path, frame, type_filter), sys.ptr)
    return sys
end

"""
    seams_chill(sys; cutoff=3.5, type=0, graph=0, rule=0)

Original CHILL labels and nine population counts. `cutoff` is angstrom.
`graph` 0 is cutoff, 1 is mutual. `rule` 0 is published, 1 is OVITO.
"""
function seams_chill(sys::System; cutoff::Real=3.5, type::Integer=0, graph::Integer=0, rule::Integer=0)
    n = _count(sys)
    labels = Vector{Int32}(undef, n)
    counts = Vector{Int32}(undef, 9)
    _seams_check(ccall(_seams_sym(:seams_chill), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Int32}, Csize_t, Ptr{Int32}),
                       sys.ptr, cutoff, type, graph, rule, labels, n, counts), sys.ptr)
    return labels, counts
end

"""
    seams_steinhardt(sys, order; cutoff=3.5, average=false)

Local Steinhardt `q_l` on the cutoff shell. `order` is 3, 4, 6, 8, or 12.
`average=true` also returns the Lechner-Dellago mean.
"""
function seams_steinhardt(sys::System, order::Integer; cutoff::Real=3.5, average::Bool=false)
    n = _count(sys)
    ql = Vector{Float64}(undef, n)
    bar = average ? Vector{Float64}(undef, n) : Float64[]
    _seams_check(ccall(_seams_sym(:seams_steinhardt), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Ptr{Float64}, Csize_t, Ptr{Float64}, Csize_t),
                       sys.ptr, cutoff, 0, order, ql, n,
                       average ? pointer(bar) : C_NULL, average ? n : 0), sys.ptr)
    return ql, bar
end

"""
    seams_cages(sys; cutoff=3.5, graph=2, k=4, complete=false)

Disjoint cage counts `(n_hex, n_cubic, n_water, n_both)`. The cage handle is freed before return.
"""
function seams_cages(sys::System; cutoff::Real=3.5, graph::Integer=2, k::Integer=4, complete::Bool=false)
    out = Ref{Ptr{Cvoid}}(C_NULL)
    _seams_check(ccall(_seams_sym(:seams_cages), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Cint, Ptr{Ptr{Cvoid}}),
                       sys.ptr, cutoff, 0, graph, k, complete ? 1 : 0, out), sys.ptr)
    n_hex = Ref{Cint}(0)
    n_cubic = Ref{Cint}(0)
    n_water = Ref{Cint}(0)
    n_both = Ref{Cint}(0)
    try
        _seams_check(ccall(_seams_sym(:seams_cages_counts), Cint,
                           (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                           out[], n_hex, n_cubic, n_water, n_both), out[])
    finally
        ccall(_seams_sym(:seams_cages_destroy), Cvoid, (Ptr{Cvoid},), out[])
    end
    return Int(n_hex[]), Int(n_cubic[]), Int(n_water[]), Int(n_both[])
end

"""
    seams_signature_cages(sys, spec; cutoff=3.5, graph=0)

Closed cages of one ring census. Returns `(n_cages, n_atoms)`.
"""
function seams_signature_cages(sys::System, spec::AbstractString; cutoff::Real=3.5, graph::Integer=0)
    n_cages = Ref{Csize_t}(0)
    n_atoms = Ref{Csize_t}(0)
    _seams_check(ccall(_seams_sym(:seams_signature_cages), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cstring, Ptr{Csize_t}, Ptr{Csize_t}),
                       sys.ptr, cutoff, 0, graph, spec, n_cages, n_atoms), sys.ptr)
    return Int(n_cages[]), Int(n_atoms[])
end

"""
    seams_neighbours_knn(sys; k=4, cutoff=3.5, mutual=true)

Number of neighbour-list columns. The handle is freed before return.
"""
function seams_neighbours_knn(sys::System; k::Integer=4, cutoff::Real=3.5, mutual::Bool=true)
    out = Ref{Ptr{Cvoid}}(C_NULL)
    _seams_check(ccall(_seams_sym(:seams_neighbours_knn), Cint,
                       (Ptr{Cvoid}, Cint, Cdouble, Cint, Cint, Ptr{Ptr{Cvoid}}),
                       sys.ptr, k, cutoff, 0, mutual ? 1 : 0, out), sys.ptr)
    n_offsets = Ref{Csize_t}(0)
    nnz = Ref{Csize_t}(0)
    try
        _seams_check(ccall(_seams_sym(:seams_neighbours_copy), Cint,
                           (Ptr{Cvoid}, Ptr{Cint}, Csize_t, Ptr{Csize_t}, Ptr{Cint}, Csize_t, Ptr{Csize_t}),
                           out[], C_NULL, 0, n_offsets, C_NULL, 0, nnz), out[])
    finally
        ccall(_seams_sym(:seams_neighbours_destroy), Cvoid, (Ptr{Cvoid},), out[])
    end
    return Int(nnz[])
end

"""
    seams_rdf(sys; rmax=6.0, bins=8)

Partial radial distribution. Returns `(r, g, volume)`. `r` is angstrom.
"""
function seams_rdf(sys::System; rmax::Real=6.0, bins::Integer=8)
    n_out = Ref{Csize_t}(0)
    volume = Ref{Float64}(0.0)
    _seams_check(ccall(_seams_sym(:seams_rdf), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Csize_t, Ptr{Csize_t}, Ptr{Float64}),
                       sys.ptr, rmax, bins, 0, 0, C_NULL, C_NULL, C_NULL, 0, n_out, volume), sys.ptr)
    r = Vector{Float64}(undef, n_out[])
    g = Vector{Float64}(undef, n_out[])
    _seams_check(ccall(_seams_sym(:seams_rdf), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Csize_t, Ptr{Csize_t}, Ptr{Float64}),
                       sys.ptr, rmax, bins, 0, 0, r, g, C_NULL, n_out[], n_out, C_NULL), sys.ptr)
    return r, g, volume[]
end

"""
    seams_running_cn(sys; rmax=6.0, bins=8)

Running coordination number along the partial RDF.
"""
function seams_running_cn(sys::System; rmax::Real=6.0, bins::Integer=8)
    n_out = Ref{Csize_t}(0)
    _seams_check(ccall(_seams_sym(:seams_running_cn), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Csize_t, Ptr{Csize_t}),
                       sys.ptr, rmax, bins, 0, 0, C_NULL, C_NULL, 0, n_out), sys.ptr)
    r = Vector{Float64}(undef, n_out[])
    cn = Vector{Float64}(undef, n_out[])
    _seams_check(ccall(_seams_sym(:seams_running_cn), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Csize_t, Ptr{Csize_t}),
                       sys.ptr, rmax, bins, 0, 0, r, cn, n_out[], n_out), sys.ptr)
    return cn
end

"""
    seams_coordination(sys; rmax=3.5, bins=8)

Coordination number integrated to `rmax` angstrom.
"""
function seams_coordination(sys::System; rmax::Real=3.5, bins::Integer=8)
    cn = Ref{Float64}(0.0)
    _seams_check(ccall(_seams_sym(:seams_coordination), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Float64}),
                       sys.ptr, rmax, bins, 0, 0, cn), sys.ptr)
    return cn[]
end

"""
    seams_hbonds(sys; cutoff=3.5, oxygen_type=1, hydrogen_type=2, dist=2.42, angle=30.0)

Hydrogen-bond count. Distances are angstrom. `angle` is degrees.
"""
function seams_hbonds(sys::System; cutoff::Real=3.5, oxygen_type::Integer=1, hydrogen_type::Integer=2,
                      dist::Real=2.42, angle::Real=30.0)
    n_bonds = Ref{Cint}(0)
    _seams_check(ccall(_seams_sym(:seams_hbonds), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cdouble, Cdouble, Ptr{Cint}),
                       sys.ptr, cutoff, oxygen_type, hydrogen_type, dist, angle, n_bonds), sys.ptr)
    return Int(n_bonds[])
end

"""
    seams_density_axis(sys; type=0, bins=4, axis=2)

Number density along one cell axis. Returns the density column.
"""
function seams_density_axis(sys::System; type::Integer=0, bins::Integer=4, axis::Integer=2)
    n_out = Ref{Csize_t}(0)
    _seams_check(ccall(_seams_sym(:seams_density_axis), Cint,
                       (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Csize_t, Ptr{Csize_t}),
                       sys.ptr, type, bins, axis, C_NULL, C_NULL, 0, n_out), sys.ptr)
    z = Vector{Float64}(undef, n_out[])
    rho = Vector{Float64}(undef, n_out[])
    _seams_check(ccall(_seams_sym(:seams_density_axis), Cint,
                       (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Csize_t, Ptr{Csize_t}),
                       sys.ptr, type, bins, axis, z, rho, n_out[], n_out), sys.ptr)
    return rho
end

"""
    seams_fingerprint(sys; cutoff=3.5, hops=2)

Frame-key length and the sum of the stored ring census. The handle is freed before return.
"""
function seams_fingerprint(sys::System; cutoff::Real=3.5, hops::Integer=2)
    out = Ref{Ptr{Cvoid}}(C_NULL)
    _seams_check(ccall(_seams_sym(:seams_fingerprint), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Cint, Cint, Ptr{Ptr{Cvoid}}),
                       sys.ptr, cutoff, 0, 4, 0, hops, 0, out), sys.ptr)
    n_out = Ref{Csize_t}(0)
    counts = zeros(Int32, 13)
    try
        _seams_check(ccall(_seams_sym(:seams_fingerprint_key), Cint,
                           (Ptr{Cvoid}, Ptr{Cchar}, Csize_t, Ptr{Csize_t}),
                           out[], C_NULL, 0, n_out), out[])
        _seams_check(ccall(_seams_sym(:seams_fingerprint_rings), Cint,
                           (Ptr{Cvoid}, Ptr{Int32}, Csize_t),
                           out[], counts, 13), out[])
    finally
        ccall(_seams_sym(:seams_fingerprint_destroy), Cvoid, (Ptr{Cvoid},), out[])
    end
    return Int(n_out[]), Int(sum(counts))
end

"""
    seams_contact_pairs(sys, spec)

Mutual nearest unlike contacts. Returns `(n_pairs, n_cation, n_anion)`.
"""
function seams_contact_pairs(sys::System, spec::AbstractString)
    n_pairs = Ref{Csize_t}(0)
    n_cation = Ref{Cint}(0)
    n_anion = Ref{Cint}(0)
    _seams_check(ccall(_seams_sym(:seams_contact_pairs), Cint,
                       (Ptr{Cvoid}, Cstring, Ptr{Csize_t}, Ptr{Cint}, Ptr{Cint}),
                       sys.ptr, spec, n_pairs, n_cation, n_anion), sys.ptr)
    return Int(n_pairs[]), Int(n_cation[]), Int(n_anion[])
end

"""
    seams_domains(sys, spec; cutoff=3.5, subset=0)

Largest domain. `subset` 0 is apolar and 1 is polar.
Returns `(n_subset, largest, percolation)`.
"""
function seams_domains(sys::System, spec::AbstractString; cutoff::Real=3.5, subset::Integer=0)
    n_subset = Ref{Cint}(0)
    largest = Ref{Cint}(0)
    percolation = Ref{Float64}(0.0)
    _seams_check(ccall(_seams_sym(:seams_domains), Cint,
                       (Ptr{Cvoid}, Cdouble, Cstring, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Float64}),
                       sys.ptr, cutoff, spec, subset, n_subset, largest, percolation), sys.ptr)
    return Int(n_subset[]), Int(largest[]), percolation[]
end

"""
    seams_prisms(sys; cutoff=3.5)

Perfect prism counts by ring size. Index 7 (1-based) is size 6.
"""
function seams_prisms(sys::System; cutoff::Real=3.5)
    perfect = zeros(Int32, 7)
    imperfect = zeros(Int32, 7)
    _seams_check(ccall(_seams_sym(:seams_prisms), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Ptr{Int32}, Ptr{Int32}, Csize_t),
                       sys.ptr, cutoff, 0, perfect, imperfect, 7), sys.ptr)
    return perfect
end

"""
    seams_monolayer(sys; cutoff=3.5, max_depth=6, sheet_area=100.0)

Primitive-ring histogram of a monolayer. `sheet_area` is square angstrom.
"""
function seams_monolayer(sys::System; cutoff::Real=3.5, max_depth::Integer=6, sheet_area::Real=100.0)
    counts = zeros(Int32, max_depth + 1)
    _seams_check(ccall(_seams_sym(:seams_monolayer), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cdouble, Ptr{Int32}, Csize_t),
                       sys.ptr, cutoff, 0, max_depth, sheet_area, counts, length(counts)), sys.ptr)
    return counts
end

"""
    seams_rings_vertices(sys; cutoff=3.5, max_depth=6, index=0)

Atom indices of one primitive ring. The ring handle is freed before return.
"""
function seams_rings_vertices(sys::System; cutoff::Real=3.5, max_depth::Integer=6, index::Integer=0)
    out = Ref{Ptr{Cvoid}}(C_NULL)
    _seams_check(ccall(_seams_sym(:seams_rings_find), Cint,
                       (Ptr{Cvoid}, Cdouble, Cint, Cint, Cint, Ptr{Ptr{Cvoid}}),
                       sys.ptr, cutoff, 0, 0, max_depth, out), sys.ptr)
    n_out = Ref{Csize_t}(0)
    try
        _seams_check(ccall(_seams_sym(:seams_rings_vertices), Cint,
                           (Ptr{Cvoid}, Cint, Ptr{Cint}, Csize_t, Ptr{Csize_t}),
                           out[], index, C_NULL, 0, n_out), out[])
        atoms = Vector{Int32}(undef, n_out[])
        _seams_check(ccall(_seams_sym(:seams_rings_vertices), Cint,
                           (Ptr{Cvoid}, Cint, Ptr{Int32}, Csize_t, Ptr{Csize_t}),
                           out[], index, atoms, n_out[], n_out), out[])
        return atoms
    finally
        ccall(_seams_sym(:seams_rings_destroy), Cvoid, (Ptr{Cvoid},), out[])
    end
end

"""
    seams_site_table(sys, spec)

Per-atom site kind from a site-table spec such as `"1=waterO"`.
"""
function seams_site_table(sys::System, spec::AbstractString)
    n = _count(sys)
    kinds = Vector{Int32}(undef, n)
    _seams_check(ccall(_seams_sym(:seams_site_table), Cint,
                       (Ptr{Cvoid}, Cstring, Ptr{Int32}, Csize_t),
                       sys.ptr, spec, kinds, n), sys.ptr)
    return kinds
end
