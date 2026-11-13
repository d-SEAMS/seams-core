# Copyright (c) d-SEAMS developers
"""
Julia front for the d-SEAMS C API. Set `SEAMS_LIB` to `libseams.so`
when it is not on the loader path.
"""
module Seams

export System, version, set_num_threads, num_threads, set_ortho!, chill_plus, q4q6, six_rings
export seams_system_count, seams_system_read, seams_chill, seams_steinhardt, seams_cages, seams_signature_cages
export seams_neighbours_knn, seams_rdf, seams_running_cn, seams_coordination, seams_hbonds
export seams_density_axis, seams_fingerprint, seams_contact_pairs, seams_domains, seams_prisms
export seams_monolayer, seams_rings_vertices, seams_site_table

include(joinpath(@__DIR__, "..", "..", "capi.jl"))

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

version() = unsafe_string(ccall((:seams_version_string, _lib()), Ptr{Cchar}, ()))

set_num_threads(n::Integer) = _check(ccall((:seams_set_num_threads, _lib()), Cint, (Cint,), Cint(n)))

num_threads() = Int(ccall((:seams_get_num_threads, _lib()), Cint, ()))

function set_ortho!(sys::System, xyz::AbstractMatrix{<:Real}, lengths::NTuple{3,Real})
    size(xyz, 1) == 3 || error("xyz is 3 by n")
    # A 3-by-n Float64 matrix is already the interleaved layout the C API copies.
    flat = xyz isa StridedMatrix{Float64} && stride(xyz, 1) == 1 && stride(xyz, 2) == 3 ? xyz :
           collect(Float64, xyz)
    _check(ccall((:seams_system_set_ortho, _lib()), Cint,
                 (Ptr{Cvoid}, Ptr{Float64}, Cint, Float64, Float64, Float64, Float64, Float64, Float64, Ptr{Cint}, Ptr{Cint}),
                 sys.ptr, flat, Cint(size(xyz, 2)),
                 0.0, 0.0, 0.0, Float64(lengths[1]), Float64(lengths[2]), Float64(lengths[3]),
                 C_NULL, C_NULL))
    return sys
end

function chill_plus(sys::System; cutoff=3.5, graph=CUTOFF, rule=PUBLISHED)
    n = Int(ccall((:seams_system_count, _lib()), Cint, (Ptr{Cvoid},), sys.ptr))
    labels = Vector{Int32}(undef, n)
    counts = Vector{Int32}(undef, 9)
    _check(ccall((:seams_chill_plus, _lib()), Cint,
                 (Ptr{Cvoid}, Float64, Cint, Cint, Cint, Ptr{Int32}, Ptr{Int32}),
                 sys.ptr, Float64(cutoff), Cint(0), Cint(graph), Cint(rule), labels, counts))
    return labels, counts
end

function q4q6(sys::System; cutoff=3.5)
    n = Int(ccall((:seams_system_count, _lib()), Cint, (Ptr{Cvoid},), sys.ptr))
    q4 = Vector{Float64}(undef, n)
    q6 = Vector{Float64}(undef, n)
    _check(ccall((:seams_q4q6, _lib()), Cint,
                 (Ptr{Cvoid}, Float64, Cint, Ptr{Float64}, Ptr{Float64}),
                 sys.ptr, Float64(cutoff), Cint(0), q4, q6))
    return q4, q6
end

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
