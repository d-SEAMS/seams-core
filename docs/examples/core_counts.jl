# Ice Ih and sII census through Dseams.
# Ice Ih: 128 hexagonal, 256 six-rings.
# sII: 136 clathrate, 144 five-rings, 16 six-rings.
using Dseams

function census(path)
    frame = read_frame(path)
    ice = chill_plus(frame; cutoff = 3.5, graph = "cutoff")
    rings = primitive_rings(frame; cutoff = 3.5, graph = "cutoff")
    ice, rings
end

ice, ice_rings = census("input/traj/ice1h.lammpstrj")
@assert ice.counts["hexagonal"] == 128
@assert ice_rings[6] == 256
sii, sii_rings = census("input/traj/genice_sII.lammpstrj")
@assert sii.counts["clathrate"] == 136
@assert sii_rings[5] == 144
@assert sii_rings[6] == 16
println("julia hexagonal 128 rings6 256 clathrate 136 rings5 144")
