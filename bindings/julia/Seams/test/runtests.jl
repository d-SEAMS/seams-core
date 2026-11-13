# Copyright (c) d-SEAMS developers
using Test
using Seams

function diamond(nrep)
    a = 2.75 * 4 / sqrt(3)
    fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
    basis = ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
    cols = NTuple{3,Float64}[]
    for i in 0:nrep-1, j in 0:nrep-1, k in 0:nrep-1, shift in fcc, site in basis
        push!(cols, ((i + shift[1] + site[1]) * a,
                     (j + shift[2] + site[2]) * a,
                     (k + shift[3] + site[3]) * a))
    end
    xyz = reduce(hcat, collect.(cols))
    return xyz, nrep * a
end

@testset "cubic ice" begin
    set_num_threads(1)
    @test num_threads() == 1
    @test !isempty(version())
    xyz, edge = diamond(2)
    sys = System()
    set_ortho!(sys, xyz, (edge, edge, edge))
    labels, counts = chill_plus(sys; rule=Seams.PUBLISHED)
    @test all(labels .== 0)
    @test counts[1] == length(labels)
    labels, counts = chill_plus(sys; graph=Seams.MUTUAL, rule=Seams.OVITO)
    @test counts[1] == length(labels)
    q4, q6 = q4q6(sys)
    @test sum(q6) / length(q6) > 0.3
    @test all(isfinite, q4)
    six, fraction = six_rings(sys)
    @test six > 0
    @test fraction > 0.5
end
