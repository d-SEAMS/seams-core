using Dseams

root = get(ENV, "DSEAMS_ROOT", "")
@assert !isempty(root)
frame = read_frame(joinpath(root, "input", "traj", "exampleTraj.lammpstrj"); type = 2)
@assert nop(frame) == 250
ice = chill_plus(frame)
@assert ice.counts["interClathrate"] == 9
@assert ice.counts["water"] == 241
order = q4_q6(frame)
@assert length(order.q6) == nop(frame)
println("julia api ok")
