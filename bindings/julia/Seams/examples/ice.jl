# Copyright (c) d-SEAMS developers
using Seams

a = 2.75 * 4 / sqrt(3)
nrep = 2
fcc = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5))
basis = ((0.0, 0.0, 0.0), (0.25, 0.25, 0.25))
cols = NTuple{3,Float64}[]
for i in 0:nrep-1, j in 0:nrep-1, k in 0:nrep-1, shift in fcc, site in basis
    push!(cols, ((i + shift[1] + site[1]) * a, (j + shift[2] + site[2]) * a, (k + shift[3] + site[3]) * a))
end
xyz = reduce(hcat, collect.(cols))
sys = System()
set_ortho!(sys, xyz, (nrep * a, nrep * a, nrep * a))
labels, counts = chill_plus(sys)
println("dseams ", version(), " cubic ", counts[1], " of ", length(labels))
