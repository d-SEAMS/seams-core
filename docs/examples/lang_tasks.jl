# Neighbours, rings, and CHILL+ from Julia.
using Dseams

const labels = (
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

frame = read_frame(joinpath(ENV["DSEAMS_ROOT"], "input", "traj", "ice1h.lammpstrj"); type = -1)
pairs = neighbour_pairs(frame; cutoff = 3.5, graph = "cutoff")
rings = primitive_rings(frame; cutoff = 3.5, graph = "cutoff", max_depth = 6)
histogram = [get(rings, size, 0) for size in 0:6]
ice = chill_plus(frame; cutoff = 3.5, graph = "cutoff")
counts = [ice.counts[name] for name in labels]
println(
    "pairs ",
    pairs,
    " rings ",
    join(histogram, " "),
    " chill ",
    join(counts, " "),
)
